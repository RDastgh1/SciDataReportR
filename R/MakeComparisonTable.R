#' Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts
#'
#' Create a label-aware comparison table using `gtsummary::tbl_summary()` with
#' optional global hypothesis tests, covariate-adjusted tests, effect sizes, and
#' pairwise comparisons.
#'
#' @details
#' Continuous variables are numeric variables with more than two unique
#' non-missing values. Numeric variables with exactly two unique values are
#' treated as dichotomous categorical variables.
#'
#' When covariates are supplied, continuous outcomes are tested using ANCOVA
#' with Type II tests. If `Parametric = FALSE`, robust HC3 covariance is used
#' for the group-level Wald test and adjusted pairwise comparisons.
#'
#' Binary categorical outcomes with covariates are tested using logistic
#' regression likelihood-ratio tests. Multicategory categorical outcomes with
#' covariates are tested using multinomial likelihood-ratio tests.
#'
#' Pairwise comparisons preserve non-standard group labels and variable names.
#'
#' @param DataFrame A data frame.
#' @param CompVariable Character scalar naming the grouping variable.
#' @param Variables Character vector of variables to summarize.
#' @param ... Optional additional variable names supplied individually.
#' @param Covariates Optional character vector of covariates for adjusted models.
#' @param ValueDigits Number of digits for descriptive statistics.
#' @param pDigits Number of digits for p-values.
#' @param AddEffectSize Logical; add effect-size columns.
#' @param EffectSizeDigits Number of digits for effect sizes.
#' @param AddPairwise Logical; add pairwise comparison columns.
#' @param PairwiseMethod P-value adjustment method. Use `"none"` for no adjustment.
#' @param Parametric Logical; use parametric tests for continuous outcomes.
#' @param ParametricDisplay Logical; display continuous summaries as mean (SD).
#' If `FALSE`, display median [IQR]. Defaults to `Parametric`.
#' @param IncludeOverallN Logical; add N column.
#' @param IncludeMissing Logical; include missing rows in summaries.
#' @param suppress_warnings Logical; suppress selected gtsummary warnings.
#' @param Referent Optional reference group for pairwise comparisons.
#' @param IncludeOverallStats Logical; add overall summary column.
#' @param ShowPositiveBinaryOnLabel Logical; for binary variables, show only the
#' positive level where identifiable.
#' @param CatMethod Categorical test method. One of `"auto"`, `"chisq"`, `"fisher"`.
#' @param MultiCatAdjusted Adjusted multicategory method. Currently
#' `"multinomial_LR"` or `"none"`.
#' @param ShowNotes Whether to show Notes column. One of `"auto"`, `"always"`,
#' or `"never"`.
#' @param NotesPosition Notes column position. One of `"last"`, `"after_test"`,
#' or `"before_pairwise"`.
#'
#' @return A `gtsummary` object.
#'
#' @examples
#' \dontrun{
#' MakeComparisonTable(
#'   DataFrame = mtcars,
#'   CompVariable = "am",
#'   Variables = c("mpg", "hp", "wt"),
#'   AddEffectSize = TRUE,
#'   AddPairwise = TRUE
#' )
#' }
#'
#' @export
MakeComparisonTable <- function(
    DataFrame,
    CompVariable = NULL,
    Variables,
    ...,
    Covariates = NULL,
    ValueDigits = 2,
    pDigits = 3,
    AddEffectSize = FALSE,
    EffectSizeDigits = 2,
    AddPairwise = FALSE,
    PairwiseMethod = "bonferroni",
    Parametric = TRUE,
    ParametricDisplay = NULL,
    IncludeOverallN = FALSE,
    IncludeMissing = FALSE,
    suppress_warnings = FALSE,
    Referent = NULL,
    IncludeOverallStats = FALSE,
    ShowPositiveBinaryOnLabel = TRUE,
    CatMethod = c("auto", "chisq", "fisher"),
    MultiCatAdjusted = c("multinomial_LR", "none"),
    ShowNotes = c("auto", "always", "never"),
    NotesPosition = c("last", "after_test", "before_pairwise")
) {

  if (is.null(ParametricDisplay)) {
    ParametricDisplay <- Parametric
  }

  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)
  ShowNotes <- match.arg(ShowNotes)
  NotesPosition <- match.arg(NotesPosition)

  extra_vars <- list(...)
  if (length(extra_vars) > 0) {
    extra_vars <- unlist(extra_vars, use.names = FALSE)
    if (!is.character(extra_vars)) {
      stop("Unnamed arguments in ... must be character variable names.")
    }
    Variables <- c(Variables, extra_vars)
  }

  if (!is.character(Variables)) {
    stop("Variables must be a character vector.")
  }

  Variables <- unique(Variables)

  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang",
    "tidyselect", "car", "emmeans", "effectsize", "sandwich"
  )

  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) {
    stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))
  }

  if (is.null(PairwiseMethod)) {
    PairwiseMethod <- "none"
  }

  valid_methods <- c("none", stats::p.adjust.methods)
  if (!is.character(PairwiseMethod) ||
      length(PairwiseMethod) != 1 ||
      !(PairwiseMethod %in% valid_methods)) {
    stop(
      "PairwiseMethod must be one of: ",
      paste(valid_methods, collapse = ", "),
      ". Got: ",
      paste(PairwiseMethod, collapse = ", ")
    )
  }

  comp_present <- !is.null(CompVariable) &&
    is.character(CompVariable) &&
    length(CompVariable) == 1 &&
    CompVariable %in% names(DataFrame)

  overall_mode <- !comp_present

  if (!all(Variables %in% names(DataFrame))) {
    stop(
      "Variable(s) not found: ",
      paste(setdiff(Variables, names(DataFrame)), collapse = ", ")
    )
  }

  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) {
    stop(
      "Covariate(s) not found: ",
      paste(setdiff(Covariates, names(DataFrame)), collapse = ", ")
    )
  }

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  if (!is.null(Covariates)) {
    drop_cov_from_vars <- intersect(Variables, Covariates)
    if (length(drop_cov_from_vars)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov_from_vars, collapse = ", "))
      Variables <- setdiff(Variables, Covariates)
    }
  }

  if (!length(Variables)) {
    stop("No variables left to summarise after removing covariates from Variables.")
  }

  as_factor_drop <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  btick <- function(x) {
    paste0("`", gsub("`", "\\\\`", as.character(x)), "`")
  }

  fmla <- function(lhs, rhs_terms) {
    rhs_terms <- as.character(rhs_terms)
    rhs <- paste(btick(rhs_terms), collapse = " + ")
    stats::as.formula(paste(btick(lhs), "~", rhs))
  }

  norm <- function(x) {
    gsub("\\s+", " ", trimws(as.character(x)))
  }

  pair_label <- function(a, b) {
    paste(a, b, sep = " - ")
  }

  pair_key <- function(a, b) {
    paste(sort(c(norm(a), norm(b))), collapse = "||")
  }

  fmt_p <- function(p, digits = 3) {
    if (is.na(p)) return(NA_character_)
    p <- max(p, .Machine$double.xmin)
    format.pval(p, digits = digits, eps = 10^-digits)
  }

  clean_tab <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    tab
  }

  cat_global_test <- function(tab, method = c("auto", "chisq", "fisher")) {
    method <- match.arg(method)
    tab <- clean_tab(tab)

    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(list(p = NA_real_, label = "Insufficient data"))
    }

    chi <- tryCatch(
      suppressWarnings(stats::chisq.test(tab, correct = FALSE)),
      error = function(e) NULL
    )

    if (method == "chisq") {
      if (!is.null(chi) && !is.null(chi$p.value)) {
        return(list(p = as.numeric(chi$p.value), label = "Pearson chi-squared"))
      }
      return(list(p = NA_real_, label = "Chi-squared failed"))
    }

    use_fisher <- method == "fisher" ||
      is.null(chi) ||
      any(chi$expected < 5)

    if (use_fisher) {
      fish <- tryCatch(
        if (nrow(tab) == 2 && ncol(tab) == 2) {
          stats::fisher.test(tab)
        } else {
          stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)
        },
        error = function(e) NULL
      )

      if (!is.null(fish) && !is.null(fish$p.value)) {
        return(list(
          p = as.numeric(fish$p.value),
          label = if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"
        ))
      }

      return(list(p = NA_real_, label = "Fisher failed"))
    }

    list(p = as.numeric(chi$p.value), label = "Pearson chi-squared")
  }

  cramers_v <- function(tab) {
    tab <- clean_tab(tab)

    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(NA_real_)
    }

    if (requireNamespace("DescTools", quietly = TRUE)) {
      return(as.numeric(DescTools::CramerV(tab, method = "bias.corrected")))
    }

    chi <- tryCatch(
      suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic),
      error = function(e) NA_real_
    )

    n <- sum(tab)
    m <- min(nrow(tab), ncol(tab)) - 1

    if (!is.finite(chi) || m <= 0 || n <= 0) {
      return(NA_real_)
    }

    sqrt(as.numeric(chi) / (n * m))
  }

  robust_group_p <- function(fit_lm, group_term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)

    mm <- tryCatch(stats::model.matrix(fit_lm), error = function(e) NULL)
    if (is.null(mm)) return(NA_real_)

    term_labels <- attr(stats::terms(fit_lm), "term.labels")
    group_pos <- match(group_term, term_labels)

    if (is.na(group_pos)) return(NA_real_)

    assign_vec <- attr(mm, "assign")
    idx <- which(assign_vec == group_pos)
    idx <- idx[idx != 1]

    if (!length(idx)) return(NA_real_)

    cn <- colnames(mm)

    L <- matrix(0, nrow = length(idx), ncol = length(cn))
    L[cbind(seq_along(idx), idx)] <- 1
    colnames(L) <- cn

    lh <- tryCatch(
      car::linearHypothesis(fit_lm, hypothesis.matrix = L, vcov. = V, test = "F"),
      error = function(e) NULL
    )

    if (is.null(lh) || !"Pr(>F)" %in% colnames(lh)) {
      return(NA_real_)
    }

    as.numeric(lh[2, "Pr(>F)"])
  }

  multinom_lr_p <- function(df_cc, y, group, covs) {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      return(list(p = NA_real_, label = "Multinomial LR (nnet missing)"))
    }

    df_cc[[y]] <- droplevels(as_factor_drop(df_cc[[y]]))

    if (nlevels(df_cc[[y]]) < 3) {
      return(list(p = NA_real_, label = "Multinomial LR (collapsed)"))
    }

    f_full <- fmla(y, c(group, covs))
    f_red <- if (length(covs)) fmla(y, covs) else stats::as.formula(paste(btick(y), "~ 1"))

    m_full <- tryCatch(nnet::multinom(f_full, data = df_cc, trace = FALSE), error = function(e) NULL)
    m_red <- tryCatch(nnet::multinom(f_red, data = df_cc, trace = FALSE), error = function(e) NULL)

    if (is.null(m_full) || is.null(m_red)) {
      return(list(p = NA_real_, label = "Multinomial LR (failed)"))
    }

    lr <- 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))
    df_lr <- attr(stats::logLik(m_full), "df") - attr(stats::logLik(m_red), "df")

    if (!is.finite(lr) || !is.finite(df_lr) || df_lr <= 0) {
      return(list(p = NA_real_, label = "Multinomial LR (invalid)"))
    }

    list(
      p = stats::pchisq(lr, df = df_lr, lower.tail = FALSE),
      label = "Multinomial LR"
    )
  }

  pairwise_cat_lr <- function(df_cc, y, group, covs, combos, PairwiseMethod) {
    out <- purrr::map_dfr(combos, function(cp) {
      a <- cp[1]
      b <- cp[2]

      sub <- df_cc[df_cc[[group]] %in% c(a, b), , drop = FALSE]
      sub[[group]] <- droplevels(as_factor_drop(sub[[group]]))
      sub[[y]] <- droplevels(as_factor_drop(sub[[y]]))

      if (nlevels(sub[[group]]) < 2 || nrow(sub) < 5) {
        return(tibble::tibble(
          key = pair_key(a, b),
          contrast_label = pair_label(a, b),
          p_val = NA_real_
        ))
      }

      k <- nlevels(sub[[y]])

      if (k < 2) {
        return(tibble::tibble(
          key = pair_key(a, b),
          contrast_label = pair_label(a, b),
          p_val = NA_real_
        ))
      }

      f_full <- fmla(y, c(group, covs))
      f_red <- if (length(covs)) fmla(y, covs) else stats::as.formula(paste(btick(y), "~ 1"))

      if (k == 2) {
        m_full <- tryCatch(stats::glm(f_full, data = sub, family = stats::binomial()), error = function(e) NULL)
        m_red <- tryCatch(stats::glm(f_red, data = sub, family = stats::binomial()), error = function(e) NULL)

        if (is.null(m_full) || is.null(m_red)) {
          return(tibble::tibble(
            key = pair_key(a, b),
            contrast_label = pair_label(a, b),
            p_val = NA_real_
          ))
        }

        an <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
        p <- if (!is.null(an) && "Pr(>Chi)" %in% names(an)) as.numeric(an$`Pr(>Chi)`[2]) else NA_real_

        return(tibble::tibble(
          key = pair_key(a, b),
          contrast_label = pair_label(a, b),
          p_val = p
        ))
      }

      if (!requireNamespace("nnet", quietly = TRUE)) {
        return(tibble::tibble(
          key = pair_key(a, b),
          contrast_label = pair_label(a, b),
          p_val = NA_real_
        ))
      }

      m_full <- tryCatch(nnet::multinom(f_full, data = sub, trace = FALSE), error = function(e) NULL)
      m_red <- tryCatch(nnet::multinom(f_red, data = sub, trace = FALSE), error = function(e) NULL)

      if (is.null(m_full) || is.null(m_red)) {
        return(tibble::tibble(
          key = pair_key(a, b),
          contrast_label = pair_label(a, b),
          p_val = NA_real_
        ))
      }

      lr <- 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))
      df_lr <- attr(stats::logLik(m_full), "df") - attr(stats::logLik(m_red), "df")

      p <- if (is.finite(lr) && is.finite(df_lr) && df_lr > 0) {
        stats::pchisq(lr, df = df_lr, lower.tail = FALSE)
      } else {
        NA_real_
      }

      tibble::tibble(
        key = pair_key(a, b),
        contrast_label = pair_label(a, b),
        p_val = p
      )
    })

    if (!identical(PairwiseMethod, "none")) {
      out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    }

    out
  }

  pairwise_cont_emmeans <- function(df_cc, var, group, covs, combos, PairwiseMethod, robust = FALSE) {
    fit <- tryCatch(
      stats::lm(fmla(var, c(group, covs)), data = df_cc),
      error = function(e) NULL
    )

    out_template <- purrr::map_dfr(combos, function(cp) {
      tibble::tibble(
        variable = var,
        contrast_label = pair_label(cp[1], cp[2]),
        key = pair_key(cp[1], cp[2]),
        p_val = NA_real_
      )
    })

    if (is.null(fit)) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    V <- NULL
    if (robust) {
      V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
    }

    emm <- tryCatch(
      if (is.null(V)) {
        emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", btick(group))))
      } else {
        emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", btick(group))), vcov. = V)
      },
      error = function(e) NULL
    )

    if (is.null(emm)) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    emm_df <- tryCatch(as.data.frame(emm), error = function(e) NULL)
    if (is.null(emm_df) || !(group %in% names(emm_df))) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    emm_levels <- as.character(emm_df[[group]])

    contrast_list <- list()

    for (cp in combos) {
      a <- cp[1]
      b <- cp[2]

      if (!(a %in% emm_levels) || !(b %in% emm_levels)) next

      v <- rep(0, length(emm_levels))
      v[match(a, emm_levels)] <- 1
      v[match(b, emm_levels)] <- -1

      contrast_list[[pair_key(a, b)]] <- v
    }

    if (!length(contrast_list)) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    ctr <- tryCatch(
      emmeans::contrast(emm, method = contrast_list),
      error = function(e) NULL
    )

    if (is.null(ctr)) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    adj_method <- if (identical(PairwiseMethod, "none")) "none" else PairwiseMethod

    res <- tryCatch(
      as.data.frame(summary(ctr, adjust = adj_method)),
      error = function(e) NULL
    )

    if (is.null(res) || !("p.value" %in% names(res))) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    got <- tibble::tibble(
      key = as.character(res$contrast),
      p_val = as.numeric(res$p.value)
    )

    out_template %>%
      dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
      dplyr::select("variable", "contrast_label", "p_val")
  }

  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)

  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))

  if (!overall_mode) {
    df[[CompVariable]] <- as_factor_drop(df[[CompVariable]])
  }

  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]

  drop <- setdiff(Variables, keep)
  if (length(drop)) {
    warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  }

  Variables <- keep

  if (!length(Variables)) {
    stop("No variables left to summarise after dropping constants.")
  }

  n_unique <- vapply(Variables, function(v) length(unique(df[[v]][!is.na(df[[v]])])), integer(1))
  is_num <- vapply(Variables, function(v) is.numeric(df[[v]]), logical(1))
  is_dich <- n_unique == 2
  treat_as_continuous <- is_num & !is_dich

  numeric_cont <- Variables[treat_as_continuous]
  dichotomous_numeric <- Variables[is_num & is_dich]

  type_list <- NULL
  if (length(numeric_cont) || length(dichotomous_numeric)) {
    type_list <- c(
      rlang::set_names(as.list(rep("continuous", length(numeric_cont))), numeric_cont),
      rlang::set_names(as.list(rep("dichotomous", length(dichotomous_numeric))), dichotomous_numeric)
    )
  }

  value_list <- NULL

  if (isTRUE(ShowPositiveBinaryOnLabel)) {
    pos_tokens <- c("TRUE", "True", "true", "1", "YES", "Yes", "yes")
    vmap <- list()

    for (v in Variables) {
      x <- df[[v]]
      ux <- unique(x[!is.na(x)])

      if (length(ux) != 2) next

      if (is.logical(x)) {
        vmap[[v]] <- TRUE
      } else if (is.numeric(x)) {
        vals <- sort(unique(as.numeric(ux)))
        vmap[[v]] <- if (any(vals == 1)) 1 else max(vals)
      } else if (is.factor(x)) {
        levs <- levels(droplevels(x))
        hit <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      } else {
        levs <- sort(unique(as.character(ux)))
        hit <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      }
    }

    if (length(vmap)) value_list <- vmap
  }

  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat <- "{n} ({p}%)"

  build_tbl_summary <- function(by_var = NULL) {
    if (is.null(by_var)) {
      return(gtsummary::tbl_summary(
        df,
        include = tidyselect::all_of(Variables),
        missing = if (IncludeMissing) "ifany" else "no",
        statistic = list(
          gtsummary::all_continuous() ~ stat_cont,
          gtsummary::all_categorical() ~ stat_cat
        ),
        digits = list(gtsummary::all_continuous() ~ ValueDigits),
        type = type_list,
        value = value_list
      ))
    }

    gtsummary::tbl_summary(
      df,
      by = tidyselect::all_of(by_var),
      include = tidyselect::all_of(Variables),
      missing = if (IncludeMissing) "ifany" else "no",
      statistic = list(
        gtsummary::all_continuous() ~ stat_cont,
        gtsummary::all_categorical() ~ stat_cat
      ),
      digits = list(gtsummary::all_continuous() ~ ValueDigits),
      type = type_list,
      value = value_list
    )
  }

  if (overall_mode) {
    tbl <- if (suppress_warnings) suppressWarnings(build_tbl_summary()) else build_tbl_summary()

    if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()

    cap <- sprintf(
      "Overall summary (display: %s; categorical global test: %s).",
      if (ParametricDisplay) "mean (SD)" else "median [IQR]",
      CatMethod
    )

    return(tbl %>% gtsummary::modify_caption(cap))
  }

  tbl <- if (suppress_warnings) {
    suppressWarnings(build_tbl_summary(CompVariable))
  } else {
    build_tbl_summary(CompVariable)
  }

  if (isTRUE(IncludeOverallStats)) {
    tbl <- tbl %>% gtsummary::add_overall(last = FALSE)
  }

  if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()

  if (nlevels(df[[CompVariable]]) < 2) {
    cap <- sprintf("Overall summary (only one level present in '%s').", CompVariable)
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  lvls_all <- levels(droplevels(as_factor_drop(df[[CompVariable]])))

  combos_all <- if (!is.null(Referent)) {
    if (!Referent %in% lvls_all) stop("Referent level not found: ", Referent)
    lapply(setdiff(lvls_all, Referent), function(x) c(Referent, x))
  } else {
    utils::combn(lvls_all, 2, simplify = FALSE)
  }

  pdat <- purrr::map_dfr(Variables, function(var) {
    is_cont <- isTRUE(treat_as_continuous[var])
    notes <- NA_character_

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))

    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(
        variable = var,
        p_unadj = NA_real_,
        p_adj = NA_real_,
        test_label = "Insufficient groups",
        Notes = "Insufficient group levels after filtering"
      ))
    }

    if (is_cont) {
      k <- nlevels(df_vg[[CompVariable]])

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))

        dropped_levels <- setdiff(lvls_all, levels(df_cc[[CompVariable]]))

        if (length(dropped_levels)) {
          notes <- paste0(
            "Adjusted analyses dropped group level(s): ",
            paste(dropped_levels, collapse = ", "),
            " (complete-case filtering on outcome + covariates)."
          )
        }

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(
            variable = var,
            p_unadj = NA_real_,
            p_adj = NA_real_,
            test_label = "Insufficient data (adjusted)",
            Notes = notes
          ))
        }

        fit <- tryCatch(stats::lm(fmla(var, c(CompVariable, Covariates)), data = df_cc), error = function(e) NULL)

        if (is.null(fit)) {
          return(tibble::tibble(
            variable = var,
            p_unadj = NA_real_,
            p_adj = NA_real_,
            test_label = "Model failed (adjusted)",
            Notes = notes
          ))
        }

        p_un <- tryCatch(
          summary(stats::aov(fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)

          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) {
            as.numeric(a2[CompVariable, "Pr(>F)"])
          } else {
            NA_real_
          }

          return(tibble::tibble(
            variable = var,
            p_unadj = p_un,
            p_adj = p_adj,
            test_label = "ANCOVA (Type II)",
            Notes = notes
          ))
        }

        p_rb <- robust_group_p(fit, CompVariable)

        return(tibble::tibble(
          variable = var,
          p_unadj = p_un,
          p_adj = p_rb,
          test_label = "Robust ANCOVA (HC3 Wald)",
          Notes = notes
        ))
      }

      if (Parametric) {
        if (k == 2) {
          p_un <- tryCatch(
            stats::t.test(fmla(var, CompVariable), data = df_vg, var.equal = FALSE)$p.value,
            error = function(e) NA_real_
          )

          return(tibble::tibble(
            variable = var,
            p_unadj = p_un,
            p_adj = NA_real_,
            test_label = "Welch t-test",
            Notes = notes
          ))
        }

        p_un <- tryCatch(
          summary(stats::aov(fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )

        return(tibble::tibble(
          variable = var,
          p_unadj = p_un,
          p_adj = NA_real_,
          test_label = "ANOVA",
          Notes = notes
        ))
      }

      p_un <- tryCatch(
        if (k == 2) {
          stats::wilcox.test(fmla(var, CompVariable), data = df_vg)$p.value
        } else {
          stats::kruskal.test(fmla(var, CompVariable), data = df_vg)$p.value
        },
        error = function(e) NA_real_
      )

      return(tibble::tibble(
        variable = var,
        p_unadj = p_un,
        p_adj = NA_real_,
        test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis",
        Notes = notes
      ))
    }

    x <- as_factor_drop(df_vg[[var]])
    g <- as_factor_drop(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- cat_global_test(tab, method = CatMethod)

    p_un <- tst$p
    test_label <- tst$label
    p_adj <- NA_real_

    if (!is.null(Covariates)) {
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(as_factor_drop(df_cc[[var]]))

      dropped_levels <- setdiff(lvls_all, levels(df_cc[[CompVariable]]))

      if (length(dropped_levels)) {
        notes <- paste0(
          "Adjusted analyses dropped group level(s): ",
          paste(dropped_levels, collapse = ", "),
          " (complete-case filtering on outcome + covariates)."
        )
      }

      if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
        df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
      }

      if (nlevels(df_cc[[var]]) == 2 &&
          nlevels(df_cc[[CompVariable]]) >= 2 &&
          nrow(df_cc) >= 5) {
        gm_full <- tryCatch(
          stats::glm(fmla(var, c(CompVariable, Covariates)), data = df_cc, family = stats::binomial()),
          error = function(e) NULL
        )

        gm_red <- tryCatch(
          if (length(Covariates)) {
            stats::glm(fmla(var, Covariates), data = df_cc, family = stats::binomial())
          } else {
            stats::glm(stats::as.formula(paste(btick(var), "~ 1")), data = df_cc, family = stats::binomial())
          },
          error = function(e) NULL
        )

        if (!is.null(gm_full) && !is.null(gm_red)) {
          an <- tryCatch(stats::anova(gm_red, gm_full, test = "Chisq"), error = function(e) NULL)
          if (!is.null(an) && "Pr(>Chi)" %in% names(an)) {
            p_adj <- as.numeric(an$`Pr(>Chi)`[2])
            test_label <- "Logistic regression (LR)"
          }
        }
      }

      if (is.na(p_adj) &&
          nlevels(df_cc[[var]]) >= 3 &&
          MultiCatAdjusted == "multinomial_LR") {
        mlr <- multinom_lr_p(df_cc, y = var, group = CompVariable, covs = Covariates)
        p_adj <- mlr$p
        test_label <- mlr$label
      }
    }

    tibble::tibble(
      variable = var,
      p_unadj = p_un,
      p_adj = p_adj,
      test_label = test_label,
      Notes = notes
    )
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(
      ~ .x %>%
        dplyr::left_join(pdat, by = "variable") %>%
        dplyr::group_by(.data$variable) %>%
        dplyr::mutate(.is_main_row = dplyr::row_number() == 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          p.value = dplyr::coalesce(.data$p_adj, .data$p_unadj),
          p.value_fmt = dplyr::if_else(
            .data$.is_main_row,
            vapply(.data$p.value, fmt_p, character(1), digits = pDigits),
            NA_character_
          ),
          Test = dplyr::if_else(.data$.is_main_row, .data$test_label, NA_character_),
          Notes = dplyr::if_else(.data$.is_main_row, .data$Notes, NA_character_)
        )
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**") %>%
    gtsummary::modify_column_unhide(columns = c("p.value_fmt", "Test"))

  tbl <- tbl %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows = .data$.is_main_row & !is.na(.data$p.value) & .data$p.value <= 0.05,
      text_format = "bold"
    )

  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {
      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))

      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(tibble::tibble(
          variable = var,
          effect_size = NA_real_,
          es_method = "Insufficient groups"
        ))
      }

      if (isTRUE(treat_as_continuous[var])) {
        if (!is.null(Covariates)) {
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))

          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            return(tibble::tibble(
              variable = var,
              effect_size = NA_real_,
              es_method = "partial eta-squared"
            ))
          }

          fit <- tryCatch(stats::lm(fmla(var, c(CompVariable, Covariates)), data = df_cc), error = function(e) NULL)
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          et <- tryCatch(effectsize::eta_squared(a2, partial = TRUE), error = function(e) NULL)

          val <- NA_real_

          if (!is.null(et)) {
            idx <- if ("Parameter" %in% names(et)) {
              et$Parameter == CompVariable
            } else {
              rownames(et) == CompVariable
            }

            val <- suppressWarnings(et$Eta2_partial[idx][1])
          }

          return(tibble::tibble(
            variable = var,
            effect_size = val,
            es_method = "partial eta-squared"
          ))
        }

        k <- nlevels(df_vg[[CompVariable]])

        if (Parametric) {
          if (k == 2) {
            val <- tryCatch(
              abs(effectsize::cohens_d(fmla(var, CompVariable), data = df_vg)$Cohens_d),
              error = function(e) NA_real_
            )

            return(tibble::tibble(
              variable = var,
              effect_size = val,
              es_method = "|d|"
            ))
          }

          val <- tryCatch(
            effectsize::eta_squared(
              stats::aov(fmla(var, CompVariable), data = df_vg),
              partial = FALSE
            )$Eta2[1],
            error = function(e) NA_real_
          )

          return(tibble::tibble(
            variable = var,
            effect_size = val,
            es_method = "eta-squared"
          ))
        }

        n <- nrow(df_vg)

        H <- tryCatch(
          stats::kruskal.test(fmla(var, CompVariable), data = df_vg)$statistic,
          error = function(e) NA_real_
        )

        eps2 <- suppressWarnings(as.numeric((H - k + 1) / (n - k)))

        return(tibble::tibble(
          variable = var,
          effect_size = eps2,
          es_method = "epsilon-squared"
        ))
      }

      tab <- table(
        as_factor_drop(df_vg[[var]]),
        as_factor_drop(df_vg[[CompVariable]])
      )

      tibble::tibble(
        variable = var,
        effect_size = cramers_v(tab),
        es_method = "Cramer's V"
      )
    })

    tbl <- tbl %>%
      gtsummary::modify_table_body(
        ~ .x %>%
          dplyr::left_join(es_df, by = "variable") %>%
          dplyr::mutate(
            effect_size = dplyr::if_else(.data$.is_main_row, .data$effect_size, NA_real_),
            ES_Method = dplyr::if_else(.data$.is_main_row, .data$es_method, NA_character_)
          )
      ) %>%
      gtsummary::modify_fmt_fun(
        effect_size ~ function(x) {
          ifelse(is.na(x), NA_character_, formatC(x, digits = EffectSizeDigits, format = "f"))
        }
      ) %>%
      gtsummary::modify_header(effect_size ~ "**Effect size**") %>%
      gtsummary::modify_header(ES_Method ~ "**ES method**") %>%
      gtsummary::modify_column_unhide(columns = c("effect_size", "ES_Method"))
  }

  if (AddPairwise && length(combos_all) > 0) {
    pw_long <- purrr::map_dfr(Variables, function(var) {
      is_cont <- isTRUE(treat_as_continuous[var])

      out_template <- purrr::map_dfr(combos_all, function(cp) {
        tibble::tibble(
          variable = var,
          contrast_label = pair_label(cp[1], cp[2]),
          key = pair_key(cp[1], cp[2]),
          p_val = NA_real_
        )
      })

      if (is_cont && !is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        return(pairwise_cont_emmeans(
          df_cc = df_cc,
          var = var,
          group = CompVariable,
          covs = Covariates,
          combos = combos_all,
          PairwiseMethod = PairwiseMethod,
          robust = !Parametric
        ))
      }

      if (!is_cont && !is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))
        df_cc[[var]] <- droplevels(as_factor_drop(df_cc[[var]]))

        got <- pairwise_cat_lr(
          df_cc,
          y = var,
          group = CompVariable,
          covs = Covariates,
          combos = combos_all,
          PairwiseMethod = PairwiseMethod
        )

        return(
          out_template %>%
            dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
            dplyr::select("variable", "contrast_label", "p_val")
        )
      }

      if (is_cont && is.null(Covariates)) {
        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))

        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        out <- purrr::map_dfr(combos_all, function(cp) {
          sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(as_factor_drop(sub[[CompVariable]]))

          p <- tryCatch(
            if (Parametric) {
              stats::t.test(fmla(var, CompVariable), data = sub, var.equal = FALSE)$p.value
            } else {
              stats::wilcox.test(fmla(var, CompVariable), data = sub)$p.value
            },
            error = function(e) NA_real_
          )

          tibble::tibble(
            variable = var,
            contrast_label = pair_label(cp[1], cp[2]),
            p_val = as.numeric(p)
          )
        })

        if (!identical(PairwiseMethod, "none")) {
          out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
        }

        return(out)
      }

      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))

      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
      }

      out <- purrr::map_dfr(combos_all, function(cp) {
        sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
        sub[[CompVariable]] <- droplevels(as_factor_drop(sub[[CompVariable]]))

        tab <- table(
          as_factor_drop(sub[[var]]),
          as_factor_drop(sub[[CompVariable]])
        )

        tst <- cat_global_test(tab, method = CatMethod)

        tibble::tibble(
          variable = var,
          contrast_label = pair_label(cp[1], cp[2]),
          p_val = as.numeric(tst$p)
        )
      })

      if (!identical(PairwiseMethod, "none")) {
        out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
      }

      out
    })

    if (nrow(pw_long) > 0) {
      contrast_levels <- unique(pw_long$contrast_label)
      safe_names <- make.unique(paste0("pw_", make.names(contrast_levels)))

      map_cols <- tibble::tibble(
        contrast_label = contrast_levels,
        col_safe = safe_names
      )

      pw_long <- dplyr::left_join(pw_long, map_cols, by = "contrast_label")

      pw_wide <- pw_long %>%
        dplyr::select("variable", "col_safe", "p_val") %>%
        tidyr::pivot_wider(
          id_cols = "variable",
          names_from = "col_safe",
          values_from = "p_val",
          values_fn = list(p_val = function(x) x[1])
        )

      tbl <- tbl %>%
        gtsummary::modify_table_body(
          ~ .x %>%
            dplyr::left_join(pw_wide, by = "variable") %>%
            dplyr::mutate(
              dplyr::across(
                tidyselect::all_of(setdiff(names(pw_wide), "variable")),
                ~ dplyr::if_else(.data$.is_main_row, as.numeric(.), NA_real_)
              )
            )
        )

      for (col in setdiff(names(pw_wide), "variable")) {
        lab <- map_cols$contrast_label[match(col, map_cols$col_safe)]

        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(
            !!rlang::sym(col) ~ function(x) vapply(x, fmt_p, character(1), digits = pDigits)
          ) %>%
          gtsummary::modify_table_styling(
            columns = col,
            rows = .data$.is_main_row & !is.na(.data[[col]]) & .data[[col]] <= 0.05,
            text_format = "bold"
          ) %>%
          gtsummary::modify_header(!!rlang::sym(col) := paste0("**", lab, "**")) %>%
          gtsummary::modify_column_unhide(columns = tidyselect::all_of(col))
      }
    }
  }

  if ("Notes" %in% names(tbl$table_body)) {
    has_notes <- any(!is.na(tbl$table_body$Notes) & nzchar(tbl$table_body$Notes))

    if (ShowNotes == "never" || (ShowNotes == "auto" && !has_notes)) {
      tbl <- tbl %>%
        gtsummary::modify_column_hide(columns = "Notes")
    } else {
      tbl <- tbl %>%
        gtsummary::modify_header(Notes ~ "**Notes**") %>%
        gtsummary::modify_column_unhide(columns = "Notes")

      if (NotesPosition == "last") {
        tbl <- tbl %>%
          gtsummary::modify_table_body(
            ~ .x %>% dplyr::relocate("Notes", .after = dplyr::last_col())
          )
      }

      if (NotesPosition == "after_test") {
        tbl <- tbl %>%
          gtsummary::modify_table_body(
            ~ .x %>% dplyr::relocate("Notes", .after = "Test")
          )
      }

      if (NotesPosition == "before_pairwise") {
        first_pw <- names(tbl$table_body)[grepl("^pw_", names(tbl$table_body))][1]

        if (!is.na(first_pw)) {
          tbl <- tbl %>%
            gtsummary::modify_table_body(
              ~ .x %>% dplyr::relocate("Notes", .before = tidyselect::all_of(first_pw))
            )
        } else {
          tbl <- tbl %>%
            gtsummary::modify_table_body(
              ~ .x %>% dplyr::relocate("Notes", .after = "Test")
            )
        }
      }
    }
  }

  cap <- sprintf(
    "Comparison table (display: %s). Global p-values: %s. Categorical global test: %s; adjusted multi-category: %s. Pairwise: %s (p-adjust: %s).",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    if (is.null(Covariates)) {
      "unadjusted (no covariates)"
    } else if (Parametric) {
      "adjusted (ANCOVA Type II / LR)"
    } else {
      "adjusted (robust ANCOVA HC3 / LR)"
    },
    CatMethod,
    MultiCatAdjusted,
    if (AddPairwise) "included" else "not included",
    PairwiseMethod
  )

  tbl %>% gtsummary::modify_caption(cap)
}
