#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, Pairwise Contrasts, and Optional Overall Summary
#'
#' @description
#' Creates a comparison table (via **gtsummary**) that
#'   • summarises variables by a grouping factor,
#'   • optionally adjusts continuous outcomes for covariates (ANCOVA),
#'   • adds p-values, effect sizes, and pair-wise contrasts,
#'   • automatically drops constant variables,
#'   • can suppress the "Unknown" row or add an overall *N* column, and
#'   • can return an overall-only summary when requested or when no valid group exists.
#'
#' @param DataFrame      Data frame with the raw data.
#' @param Variables      Character vector of columns to compare.
#' @param CompVariable   Grouping / comparison variable (character scalar). Optional if `IncludeOverallStats = TRUE`.
#' @param Covariates     Optional covariates (default `NULL`).
#' @param ValueDigits    Digits in summary stats (default `2`).
#' @param pDigits        Digits in p-values (default `3`).
#' @param AddEffectSize  Logical; add effect sizes? (default `FALSE`).
#' @param EffectSizeDigits Digits in effect sizes (default `2`).
#' @param AddPairwise    Logical; add pair-wise contrasts? (default `FALSE`).
#' @param PairwiseMethod P-adjustment method (default `"bonferroni"`).
#' @param Parametric     `TRUE` = parametric tests, `FALSE` = non-parametric tests (default `TRUE`).
#' @param ParametricDisplay `TRUE` = show mean/SD, `FALSE` = show median/IQR (default matches `Parametric`).
#' @param IncludeOverallN Logical; add overall *N*? (default `FALSE`).
#' @param IncludeMissing Logical; include "Unknown" row? (default `FALSE`).
#' @param suppress_warnings Suppress gtsummary warnings? (default `FALSE`).
#' @param Referent       Optional reference level for contrasts.
#' @param IncludeOverallStats Logical; if `TRUE` (default `FALSE`) or if `CompVariable`
#'                            is not present in `DataFrame`, return overall-only stats.
#' @param ShowPositiveBinaryOnLabel Logical; if `TRUE` (default) show only the
#'                            "positive" level (TRUE/1/YES/Yes) of binary
#'                            categorical/factor/logical variables on the label row.
#'
#' @return A **gtsummary::tbl_summary** object.
#' @export
MakeComparisonTable <- function(
    DataFrame,
    CompVariable         = NULL,
    Variables,
    Covariates           = NULL,
    ValueDigits          = 2,
    pDigits              = 3,
    AddEffectSize        = FALSE,
    EffectSizeDigits     = 2,
    AddPairwise          = FALSE,
    PairwiseMethod       = "bonferroni",
    Parametric           = TRUE,
    ParametricDisplay    = NULL,
    IncludeOverallN      = FALSE,
    IncludeMissing       = FALSE,
    suppress_warnings    = FALSE,
    Referent             = NULL,
    IncludeOverallStats  = FALSE,
    ShowPositiveBinaryOnLabel = TRUE
) {

  ## ── Set ParametricDisplay default if not specified ─────────────────────
  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  ## ── required packages ─────────────────────────────────────────────────
  req_pkgs <- c(
    "gtsummary", "dplyr", "car", "emmeans", "broom", "effectsize",
    "purrr", "tidyr", "rlang", "tibble", "tidyselect"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  ## ── mode detection ─────────────────────────────────────────────────────
  comp_present <- !is.null(CompVariable) && is.character(CompVariable) &&
    length(CompVariable) == 1 && CompVariable %in% names(DataFrame)
  overall_mode <- isTRUE(IncludeOverallStats) || !comp_present

  ## ── sanity checks ──────────────────────────────────────────────────────
  if (!overall_mode && !CompVariable %in% names(DataFrame))
    stop("Grouping variable not found: ", CompVariable)
  if (!all(Variables %in% names(DataFrame)))
    stop("Variable(s) not found: ",
         paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame)))
    stop("Covariate(s) not found: ",
         paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))

  ## ── helper: safe numeric coercion ──────────────────────────────────────
  .as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

  ## ── helper: render-safe eligibility checks ─────────────────────────────
  .eligible_numeric_by_group <- function(df_var, var, grp, min_per_group = 2) {
    x <- df_var[[var]]
    g <- df_var[[grp]]
    x <- .as_num(x)

    if (all(is.na(x))) return(FALSE)
    if (length(unique(x[!is.na(x)])) < 2) return(FALSE)

    # group-wise checks
    lv <- levels(droplevels(g))
    if (length(lv) < 2) return(FALSE)

    for (L in lv) {
      xs <- x[g == L]
      if (sum(!is.na(xs)) < min_per_group) return(FALSE)
      sdx <- suppressWarnings(stats::sd(xs, na.rm = TRUE))
      if (is.na(sdx) || !is.finite(sdx) || sdx == 0) return(FALSE)
    }

    # global Inf/NaN check
    if (any(!is.finite(x), na.rm = TRUE)) return(FALSE)

    TRUE
  }

  .eligible_categ_by_group <- function(df_var, var, grp) {
    x <- df_var[[var]]
    g <- droplevels(df_var[[grp]])

    # if all missing, not eligible
    if (all(is.na(x))) return(FALSE)

    tab <- table(droplevels(as.factor(x)), g)
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]

    # need at least 2 categories and 2 groups with data
    if (nrow(tab) < 2 || ncol(tab) < 2) return(FALSE)

    TRUE
  }

  ## ── data prep ──────────────────────────────────────────────────────────
  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- factor(df[[CompVariable]])

  ## drop constants (global) ------------------------------------------------
  keep <- Variables[sapply(Variables, function(v)
    length(unique(df[[v]][!is.na(df[[v]])])) >= 2)]
  drop <- setdiff(Variables, keep)
  if (length(drop))
    warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  Variables <- keep
  if (!length(Variables))
    stop("No variables left to summarize after dropping constants.")

  ## types for continuous vars ---------------------------------------------
  numeric_vars <- Variables[sapply(df[Variables], is.numeric)]
  type_list <- if (length(numeric_vars) > 0) {
    rlang::set_names(as.list(rep("continuous", length(numeric_vars))), numeric_vars)
  } else NULL

  ## ── binary-positive display map for tbl_summary(value=) -----------------
  value_list <- NULL
  if (isTRUE(ShowPositiveBinaryOnLabel)) {
    pos_tokens <- c("TRUE", "1", "YES", "Yes")
    vmap <- list()
    for (v in Variables) {
      x <- df[[v]]
      if (is.numeric(x)) next
      ux <- unique(x[!is.na(x)])
      if (length(ux) != 2) next

      if (is.logical(x)) {
        vmap[[v]] <- TRUE
      } else if (is.factor(x)) {
        levs <- levels(droplevels(x))
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      } else {
        levs <- sort(unique(as.character(ux)))
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      }
    }
    if (length(vmap)) value_list <- vmap
  }

  ## ── statistic templates ────────────────────────────────────────────────
  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  ## ── OVERALL-ONLY MODE ──────────────────────────────────────────────────
  if (overall_mode) {
    tbl <- gtsummary::tbl_summary(
      df,
      include   = tidyselect::all_of(Variables),
      missing   = if (IncludeMissing) "ifany" else "no",
      statistic = list(
        gtsummary::all_continuous()  ~ stat_cont,
        gtsummary::all_categorical() ~ stat_cat
      ),
      digits    = list(gtsummary::all_continuous() ~ ValueDigits),
      type      = type_list,
      value     = value_list
    )
    if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
    if (suppress_warnings) tbl <- suppressWarnings(tbl)

    cap <- sprintf("Overall Summary (display: %s, tests: %s)",
                   if (ParametricDisplay) "parametric" else "non-parametric",
                   if (Parametric) "parametric" else "non-parametric")
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  ## ── GROUPED SUMMARY ────────────────────────────────────────────────────
  tbl <- gtsummary::tbl_summary(
    df,
    by        = CompVariable,
    include   = tidyselect::all_of(Variables),
    missing   = if (IncludeMissing) "ifany" else "no",
    statistic = list(
      gtsummary::all_continuous()  ~ stat_cont,
      gtsummary::all_categorical() ~ stat_cat
    ),
    digits    = list(gtsummary::all_continuous() ~ ValueDigits),
    type      = type_list,
    value     = value_list
  )
  if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if (suppress_warnings) tbl <- suppressWarnings(tbl)

  if (nlevels(df[[CompVariable]]) < 2) {
    cap <- sprintf("Overall Summary (only one level in '%s')", CompVariable)
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  ## ── helpers for downstream steps ───────────────────────────────────────
  .cramers_v_rc <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
    if (requireNamespace("DescTools", quietly = TRUE))
      return(as.numeric(DescTools::CramerV(tab, method = "bias.corrected")))
    chi <- suppressWarnings(stats::chisq.test(tab)$statistic)
    nt  <- sum(tab); m <- min(nrow(tab), ncol(tab)) - 1
    sqrt(as.numeric(chi) / (nt * m))
  }
  .pair_label      <- function(a, b) paste(a, b, sep = " - ")
  .canonical_label <- function(lab) {
    p <- trimws(strsplit(lab, " - ", fixed = TRUE)[[1]])
    .pair_label(sort(p)[1], sort(p)[2])
  }

  ## ── p-values & test labels (render-safe) ───────────────────────────────
  pdat <- purrr::map_dfr(Variables, function(var) {

    # start with rows where group and var are available for that var's test
    df_var <- df[!is.na(df[[CompVariable]]) & !is.na(df[[var]]), , drop = FALSE]
    df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
    k <- nlevels(df_var[[CompVariable]])

    if (k < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = "Insufficient groups"))
    }

    pu  <- NA_real_
    pad <- NA_real_
    tl  <- NA_character_

    if (is.numeric(df[[var]])) {
      # For numeric, enforce per-group nonmissing + variance checks before any test call
      if (!.eligible_numeric_by_group(df_var, var, CompVariable, min_per_group = 2)) {
        return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                              test_label = "Insufficient numeric data"))
      }

      fmla_unadj <- reformulate(CompVariable, response = as.name(var))

      if (Parametric) {
        if (!is.null(Covariates)) {
          # complete-case for covariates too (prevents NA-induced weirdness)
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df_var[stats::complete.cases(df_var[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(df_cc[[CompVariable]])
          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            pu <- pad <- NA_real_; tl <- "Insufficient data (ANCOVA)"
          } else {
            fmla_adj <- reformulate(c(CompVariable, Covariates), response = as.name(var))
            m1 <- tryCatch(stats::lm(fmla_adj, data = df_cc), error = function(e) NULL)
            pu <- tryCatch(
              summary(stats::aov(fmla_unadj, data = df_cc))[[1]][CompVariable, "Pr(>F)"],
              error = function(e) NA_real_
            )
            pad <- if (is.null(m1)) NA_real_ else tryCatch(
              car::Anova(m1, type = 2)[CompVariable, "Pr(>F)"],
              error = function(e) NA_real_
            )
            tl <- "ANCOVA (Type II)"
          }

        } else if (k == 2) {
          grp_n <- table(df_var[[CompVariable]])
          if (min(grp_n) < 2) {
            pu <- pad <- NA_real_; tl <- "Insufficient data (t-test)"
          } else {
            pu <- tryCatch(stats::t.test(fmla_unadj, data = df_var, var.equal = FALSE)$p.value,
                           error = function(e) NA_real_)
            pad <- NA_real_
            tl  <- "Welch t-test"
          }

        } else {
          pu  <- tryCatch(summary(stats::aov(fmla_unadj, data = df_var))[[1]][CompVariable, "Pr(>F)"],
                          error = function(e) NA_real_)
          pad <- NA_real_
          tl  <- "ANOVA"
        }

      } else {
        pu <- tryCatch(
          if (k == 2) stats::wilcox.test(fmla_unadj, data = df_var)$p.value
          else stats::kruskal.test(fmla_unadj, data = df_var)$p.value,
          error = function(e) NA_real_
        )
        pad <- NA_real_
        tl  <- if (k == 2) "Wilcoxon rank-sum" else "Kruskal–Wallis"
      }

    } else {
      # categorical
      df_cat <- df[!is.na(df[[CompVariable]]) & !is.na(df[[var]]), , drop = FALSE]
      df_cat[[CompVariable]] <- droplevels(df_cat[[CompVariable]])

      if (!.eligible_categ_by_group(df_cat, var, CompVariable)) {
        return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                              test_label = "Insufficient categorical data"))
      }

      tbl0 <- table(droplevels(as.factor(df_cat[[var]])), df_cat[[CompVariable]])
      tbl0 <- tbl0[rowSums(tbl0) > 0, , drop = FALSE]
      tbl0 <- tbl0[, colSums(tbl0) > 0, drop = FALSE]

      if (Parametric && !is.null(Covariates) &&
          nlevels(droplevels(as.factor(df_cat[[var]]))) == 2) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df_cat[stats::complete.cases(df_cat[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(df_cc[[CompVariable]])

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 5) {
          pu <- pad <- NA_real_; tl <- "Insufficient data (logistic)"
        } else {
          fmla_adj <- reformulate(c(CompVariable, Covariates), response = as.name(var))
          gm <- tryCatch(stats::glm(fmla_adj, df_cc, family = stats::binomial()), error = function(e) NULL)
          pad <- if (is.null(gm)) NA_real_ else tryCatch(
            stats::drop1(gm, test = "Chisq")[CompVariable, "Pr(>Chi)"],
            error = function(e) NA_real_
          )
          pu <- tryCatch(suppressWarnings(stats::chisq.test(tbl0)$p.value), error = function(e) NA_real_)
          tl <- "Logistic regression (LR)"
        }

      } else {
        expct <- tryCatch(suppressWarnings(stats::chisq.test(tbl0)$expected), error = function(e) NULL)
        if (is.null(expct)) {
          pu <- NA_real_; tl <- "Chi-squared (failed)"
        } else if (any(expct < 5)) {
          pu <- tryCatch(stats::fisher.test(tbl0, simulate.p.value = TRUE, B = 1e4)$p.value,
                         error = function(e) NA_real_)
          tl <- "Fisher (sim.)"
        } else {
          pu <- tryCatch(stats::chisq.test(tbl0)$p.value, error = function(e) NA_real_)
          tl <- "Chi-squared"
        }
        pad <- NA_real_
      }
    }

    tibble::tibble(variable = var, p_unadj = pu, p_adj = pad, test_label = tl)
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::mutate(
                                     p.value     = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value_fmt = gtsummary::style_pvalue(.data$p.value, digits = pDigits),
                                     p.value_fmt = ifelse(.data$row_type == "label", .data$p.value_fmt, NA_character_),
                                     Test        = ifelse(.data$row_type == "label", .data$test_label, NA_character_)
                                   )) %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows    = .data$p.value <= 0.05 & .data$row_type == "label",
      text_format = "bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  ## ── effect sizes (render-safe) ─────────────────────────────────────────
  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {

      df_var <- df[!is.na(df[[CompVariable]]) & !is.na(df[[var]]), , drop = FALSE]
      df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
      k <- nlevels(df_var[[CompVariable]])
      n <- nrow(df_var)

      if (n < 2 || k < 2) {
        return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = NA_character_))
      }

      es_val <- NA_real_
      method <- NA_character_

      if (is.numeric(df[[var]])) {

        if (!.eligible_numeric_by_group(df_var, var, CompVariable, min_per_group = 2)) {
          return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "Insufficient numeric data"))
        }

        if (Parametric && k == 2 && is.null(Covariates)) {
          es_val <- tryCatch({
            # cohens_d sometimes trips internal stderr logic if data are degenerate, so we prechecked + still tryCatch
            abs(effectsize::cohens_d(
              reformulate(CompVariable, response = as.name(var)),
              data = df_var
            )$Cohens_d)
          }, error = function(e) NA_real_)
          method <- "|d|"

        } else if (Parametric && !is.null(Covariates)) {
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df_var[stats::complete.cases(df_var[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(df_cc[[CompVariable]])

          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            es_val <- NA_real_; method <- "partial η²"
          } else {
            fmla_es <- reformulate(c(CompVariable, Covariates), response = as.name(var))
            lm_es <- tryCatch(stats::lm(fmla_es, df_cc), error = function(e) NULL)
            if (!is.null(lm_es)) {
              a2 <- tryCatch(car::Anova(lm_es, type = 2), error = function(e) NULL)
              if (!is.null(a2)) {
                et <- tryCatch(effectsize::eta_squared(a2, partial = TRUE), error = function(e) NULL)
                if (!is.null(et)) {
                  idx <- if ("Parameter" %in% names(et)) et$Parameter == CompVariable else rownames(et) == CompVariable
                  es_val <- suppressWarnings(et$Eta2_partial[idx][1])
                }
              }
            }
            method <- "partial η²"
          }

        } else if (Parametric) {
          fmla_es <- reformulate(CompVariable, response = as.name(var))
          es_val <- tryCatch(
            effectsize::eta_squared(stats::aov(fmla_es, df_var), partial = FALSE)$Eta2[1],
            error = function(e) NA_real_
          )
          method <- "η²"

        } else {
          H <- tryCatch(stats::kruskal.test(reformulate(CompVariable, response = as.name(var)), data = df_var)$statistic,
                        error = function(e) NA_real_)
          es_val <- suppressWarnings(as.numeric((H - k + 1) / (n - k)))
          method <- "ε²"
        }

      } else {
        df_cat <- df[!is.na(df[[CompVariable]]) & !is.na(df[[var]]), , drop = FALSE]
        df_cat[[CompVariable]] <- droplevels(df_cat[[CompVariable]])
        tab <- table(droplevels(as.factor(df_cat[[var]])), df_cat[[CompVariable]])
        es_val <- .cramers_v_rc(tab)
        method <- "Cramer's V"
      }

      tibble::tibble(
        variable    = var,
        effect_size = ifelse(is.finite(es_val),
                             round(es_val, EffectSizeDigits), NA_real_),
        es_method   = method
      )
    })

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     dplyr::left_join(es_df, by = "variable") %>%
                                     dplyr::mutate(
                                       effect_size = dplyr::if_else(.data$row_type == "label", .data$effect_size, NA_real_),
                                       ES_Method   = dplyr::if_else(.data$row_type == "label", .data$es_method, NA_character_)
                                     )) %>%
      gtsummary::modify_fmt_fun(effect_size ~ function(x)
        ifelse(is.na(x), NA_character_, formatC(x, digits = EffectSizeDigits, format = "f"))) %>%
      gtsummary::modify_header(effect_size ~ "**Effect Size**") %>%
      gtsummary::modify_header(ES_Method   ~ "**ES Method**")
  }

  ## ── pair-wise contrasts (leave as-is, but wrapped more defensively) ─────
  if (AddPairwise && nlevels(df[[CompVariable]]) > 1) {

    lvls   <- levels(df[[CompVariable]])
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% lvls) stop("Referent level not found: ", Referent)
      lapply(setdiff(lvls, Referent), function(x) c(Referent, x))
    } else utils::combn(lvls, 2, simplify = FALSE)

    pw <- purrr::map_dfr(Variables, function(var) {

      empty_pw <- tibble::tibble(variable = character(), contrast = character(), p_val = numeric())

      if (length(combos) == 0) return(empty_pw)

      if (is.numeric(df[[var]])) {

        # if numeric is not eligible for group comparisons, skip pairwise safely
        df_var <- df[!is.na(df[[CompVariable]]) & !is.na(df[[var]]), , drop = FALSE]
        df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
        if (nlevels(df_var[[CompVariable]]) < 2) return(empty_pw)
        if (!.eligible_numeric_by_group(df_var, var, CompVariable, min_per_group = 2)) return(empty_pw)

        if (Parametric) {
          res <- tryCatch(
            stats::pairwise.t.test(
              x = df_var[[var]],
              g = df_var[[CompVariable]],
              p.adjust.method = PairwiseMethod,
              pool.sd = FALSE
            ),
            error = function(e) NULL
          )
          if (is.null(res) || is.null(res$p.value)) return(empty_pw)
          r <- as.data.frame(as.table(res$p.value))
          if (!all(c("Var1","Var2","Freq") %in% names(r))) return(empty_pw)
          r <- r %>%
            dplyr::filter(!is.na(Freq)) %>%
            dplyr::transmute(
              variable = var,
              g1 = as.character(Var1),
              g2 = as.character(Var2),
              contrast = .canonical_label(.pair_label(g1, g2)),
              p_val    = Freq
            )
          if (!is.null(Referent))
            r <- dplyr::filter(r, g1 == Referent | g2 == Referent)
          return(dplyr::select(r, variable, contrast, p_val))
        } else {
          purrr::map_dfr(combos, function(cp) {
            sub <- df_var[df_var[[CompVariable]] %in% cp, , drop = FALSE]
            sub[[CompVariable]] <- droplevels(sub[[CompVariable]])
            if (nlevels(sub[[CompVariable]]) < 2) return(empty_pw)
            out <- tryCatch(
              stats::pairwise.wilcox.test(
                x = sub[[var]],
                g = sub[[CompVariable]],
                p.adjust.method = PairwiseMethod
              )$p.value,
              error = function(e) NULL
            )
            if (is.null(out)) return(empty_pw)
            rr <- as.data.frame(as.table(out))
            if (!all(c("Var1","Var2","Freq") %in% names(rr))) return(empty_pw)
            rr %>%
              dplyr::transmute(
                variable = var,
                contrast = .canonical_label(.pair_label(Var1, Var2)),
                p_val    = Freq
              )
          })
        }

      } else {
        purrr::map_dfr(combos, function(cp) {
          sub <- df[df[[CompVariable]] %in% cp & !is.na(df[[var]]), , drop = FALSE]
          sub[[CompVariable]] <- droplevels(sub[[CompVariable]])
          if (nlevels(sub[[CompVariable]]) < 2) return(empty_pw)

          tbl0 <- table(droplevels(as.factor(sub[[var]])), sub[[CompVariable]])
          tbl0 <- tbl0[rowSums(tbl0) > 0, , drop = FALSE]
          tbl0 <- tbl0[, colSums(tbl0) > 0, drop = FALSE]
          if (nrow(tbl0) < 2 || ncol(tbl0) < 2) {
            pv <- NA_real_
          } else {
            expct <- tryCatch(suppressWarnings(stats::chisq.test(tbl0)$expected), error = function(e) NULL)
            pv <- if (is.null(expct)) NA_real_
            else if (any(expct < 5)) tryCatch(stats::fisher.test(tbl0, simulate.p.value = TRUE, B = 1e4)$p.value,
                                              error = function(e) NA_real_)
            else tryCatch(stats::chisq.test(tbl0)$p.value, error = function(e) NA_real_)
          }
          tibble::tibble(
            variable = var,
            contrast = .canonical_label(.pair_label(cp[1], cp[2])),
            p_val    = pv
          )
        }) %>%
          dplyr::mutate(
            p_val = if (PairwiseMethod != "none") stats::p.adjust(p_val, PairwiseMethod) else p_val
          )
      }
    })

    if (nrow(pw) > 0) {
      pw_wide <- tidyr::pivot_wider(
        pw, id_cols = "variable",
        names_from = "contrast",
        values_from = "p_val",
        values_fn = list(p_val = ~ mean(.x, na.rm = TRUE))
      )

      tbl <- tbl %>%
        gtsummary::modify_table_body(~ .x %>%
                                       dplyr::left_join(pw_wide, by = "variable") %>%
                                       dplyr::mutate(
                                         dplyr::across(
                                           tidyselect::all_of(setdiff(names(pw_wide), "variable")),
                                           ~ ifelse(.data$row_type == "label", ., NA_real_)
                                         )
                                       ))

      for (col in setdiff(names(pw_wide), "variable")) {
        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(
            !!rlang::sym(col) ~ function(x) gtsummary::style_pvalue(x, digits = pDigits)
          ) %>%
          gtsummary::modify_table_styling(
            columns     = col,
            rows        = .data[[col]] <= 0.05 & .data$row_type == "label",
            text_format = "bold"
          ) %>%
          gtsummary::modify_header(
            !!rlang::sym(col) := paste0("**", col, "**")
          ) %>%
          gtsummary::modify_footnote(
            !!rlang::sym(col) ~ paste0("Pair-wise p-value (", PairwiseMethod,
                                       if (!is.null(Referent)) paste0("; vs ", Referent), ")")
          )
      }
    }
  }

  ## ── caption & return ────────────────────────────────────────────────────
  cap <- sprintf(
    "Comparison Table (display: %s, tests: %s%s) — see **Test** column for per-variable method",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    if (Parametric) "parametric" else "non-parametric",
    if (!is.null(Covariates) && Parametric)
      paste0("; adjusted for ", paste(Covariates, collapse = ", "))
    else ""
  )
  tbl %>% gtsummary::modify_caption(cap)
}
