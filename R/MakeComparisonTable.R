#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, and Pairwise Contrasts
#'
#' @description
#' Creates a comparison table using {gtsummary} summarizing continuous and categorical variables
#' by a grouping factor, with optional covariate adjustment for continuous outcomes only (ANCOVA),
#' effect-size calculations, pairwise contrasts with p-value adjustments (or logistic regression for binary outcomes),
#' and optional suppression of the "Unknown" (missing) row. Variables with fewer than
#' two unique non-NA values are dropped automatically before any tests.
#'
#' @param DataFrame A `data.frame` containing the raw data.
#' @param Variables Character vector of column names to include in the comparison.
#' @param CompVariable Character string specifying the grouping (comparison) variable.
#' @param Covariates Optional character vector of covariate column names for adjustment (default `NULL`).
#' @param ValueDigits Integer; number of digits for summary statistics (default `2`).
#' @param pDigits Integer; number of digits for formatted p-values (default `3`).
#' @param AddEffectSize Logical; include an effect-size column? (default `FALSE`).
#' @param EffectSizeDigits Integer; number of digits for effect sizes (default `2`).
#' @param AddPairwise Logical; include pairwise contrast columns? (default `FALSE`).
#' @param PairwiseMethod Character; p-value adjustment method for contrasts (default `"bonferroni"`; `"none"` for unadjusted).
#' @param Parametric Logical; use parametric tests if `TRUE` (ANOVA/ANCOVA or Welch t-test), or nonparametric if `FALSE` (Kruskal–Wallis / Wilcoxon) (default `TRUE`).
#' @param IncludeOverallN Logical; include a column with the overall N in the table? (default `FALSE`).
#' @param IncludeMissing Logical; include a row summarizing missing data ("Unknown")? (default `FALSE`).
#' @param suppress_warnings Logical; suppress intermediate warnings from gtsummary (default `FALSE`).
#' @param Referent Optional character string specifying the reference level for pairwise contrasts.
#'
#' @return A `gtsummary::tbl_summary` object with p-values, a per-variable **Test** label,
#'         effect sizes (incl. Cramér's V for categorical), and optional pairwise contrasts.
#' @export

MakeComparisonTable <- function(
    DataFrame,
    Variables,
    CompVariable,
    Covariates        = NULL,
    ValueDigits       = 2,
    pDigits           = 3,
    AddEffectSize     = FALSE,
    EffectSizeDigits  = 2,
    AddPairwise       = FALSE,
    PairwiseMethod    = "bonferroni",
    Parametric        = TRUE,
    IncludeOverallN   = FALSE,
    IncludeMissing    = FALSE,
    suppress_warnings = FALSE,
    Referent          = NULL
) {

  ## ── packages ────────────────────────────────────────────────────────────
  req_pkgs <- c(
    "gtsummary", "dplyr", "car", "emmeans", "broom", "effectsize",
    "purrr", "tidyr", "rlang", "tibble", "tidyselect"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok))
    stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  ## ── helpers ─────────────────────────────────────────────────────────────
  .cramers_v_rc <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
    if (requireNamespace("DescTools", quietly = TRUE))
      return(as.numeric(DescTools::CramerV(tab, method = "bias.corrected")))
    chi <- suppressWarnings(stats::chisq.test(tab)$statistic)
    n_tot <- sum(tab); m <- min(nrow(tab), ncol(tab)) - 1
    sqrt(as.numeric(chi) / (n_tot * m))
  }
  .pair_label      <- function(a, b) paste(a, b, sep = " - ")
  .canonical_label <- function(lab) {
    parts <- trimws(strsplit(lab, " - ", fixed = TRUE)[[1]])
    .pair_label(sort(parts)[1], sort(parts)[2])
  }

  ## ── validation ──────────────────────────────────────────────────────────
  if (!CompVariable %in% names(DataFrame))
    stop("Grouping var not found: ", CompVariable)
  if (!all(Variables %in% names(DataFrame)))
    stop("Variable(s) not found: ",
         paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame)))
    stop("Covariate(s) not found: ",
         paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))

  ## ── data prep ───────────────────────────────────────────────────────────
  df <- DataFrame %>%
    dplyr::select(tidyselect::all_of(c(CompVariable, Variables, Covariates)))
  df[[CompVariable]] <- factor(df[[CompVariable]])

  safe_comp   <- paste0("`", CompVariable, "`")
  numeric_vars <- Variables[sapply(df[Variables], is.numeric)]
  type_list <- if (length(numeric_vars)) {
    rlang::set_names(as.list(rep("continuous", length(numeric_vars))),
                     numeric_vars)
  } else list()

  ## drop constants ----------------------------------------------------------
  keep <- Variables[sapply(Variables, function(v)
    length(unique(df[[v]][!is.na(df[[v]])])) >= 2)]

  drop <- setdiff(Variables, keep)
  if (length(drop))
    warning("Dropping: ", paste(drop, collapse = ", "))

  Variables <- keep
  if (!length(Variables))
    stop("No variables to compare after dropping constants.")

  ## ── base gtsummary table ────────────────────────────────────────────────
  tbl <- gtsummary::tbl_summary(
    df,
    by      = CompVariable,
    include = tidyselect::all_of(Variables),
    missing = if (IncludeMissing) "ifany" else "no",
    statistic = list(
      gtsummary::all_continuous()  ~
        if (Parametric) "{mean} ({sd})" else "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(gtsummary::all_continuous() ~ ValueDigits),
    type   = type_list
  )
  if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if (suppress_warnings) tbl <- suppressWarnings(tbl)
  if (nlevels(df[[CompVariable]]) < 2) return(tbl)

  ## ── p-values & test labels ───────────────────────────────────────────────
  pdat <- purrr::map_dfr(Variables, function(var) {
    df_var <- df[!is.na(df[[var]]), , drop = FALSE]
    df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
    safe_var <- paste0("`", var, "`")
    k <- nlevels(df_var[[CompVariable]])

    if (k < 2)
      return(tibble::tibble(variable = var,
                            p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = NA_character_))

    if (is.numeric(df_var[[var]])) {                         # continuous
      fmla_unadj <- stats::as.formula(paste0(safe_var, " ~ ", safe_comp))

      if (Parametric) {
        if (!is.null(Covariates)) {
          cov_terms <- paste(c(safe_comp, Covariates), collapse = " + ")
          fmla_adj  <- stats::as.formula(paste0(safe_var, " ~ ", cov_terms))
          m1   <- stats::lm(fmla_adj, data = df_var)
          pu   <- summary(stats::aov(fmla_unadj, data = df_var))[[1]][
            CompVariable, "Pr(>F)"]
          pad  <- car::Anova(m1, type = 2)[CompVariable, "Pr(>F)"]
          tl   <- "ANCOVA (Type II)"

        } else if (k == 2) {                                  # Welch t-test
          grp_n <- table(df_var[[CompVariable]])
          if (min(grp_n) < 2) {
            pu  <- NA_real_; pad <- NA_real_; tl <- "Insufficient data"
          } else {
            pu  <- stats::t.test(fmla_unadj, data = df_var,
                                 var.equal = FALSE)$p.value
            pad <- NA_real_; tl <- "Welch t-test"
          }

        } else {
          pu  <- summary(stats::aov(fmla_unadj, data = df_var))[[1]][
            CompVariable, "Pr(>F)"]
          pad <- NA_real_; tl <- "ANOVA"
        }

      } else {                                              # non-parametric
        if (k == 2) {
          pu <- stats::wilcox.test(fmla_unadj, data = df_var)$p.value
          tl <- "Wilcoxon rank-sum"
        } else {
          pu <- stats::kruskal.test(fmla_unadj, data = df_var)$p.value
          tl <- "Kruskal–Wallis"
        }
        pad <- NA_real_
      }

    } else {                                                # categorical
      tbl0 <- table(df_var[[var]], df_var[[CompVariable]])
      tbl0 <- tbl0[rowSums(tbl0) > 0, , drop = FALSE]
      tbl0 <- tbl0[, colSums(tbl0) > 0, drop = FALSE]

      if (nrow(tbl0) < 2 || ncol(tbl0) < 2) {
        pu <- pad <- NA_real_; tl <- "Insufficient categories"

      } else if (Parametric && !is.null(Covariates) &&
                 nlevels(droplevels(df_var[[var]])) == 2) {
        cov_terms <- paste(c(safe_comp, Covariates), collapse = " + ")
        fmla_adj  <- stats::as.formula(paste0(safe_var, " ~ ", cov_terms))
        gm   <- stats::glm(fmla_adj, df_var, family = stats::binomial())
        pad  <- stats::drop1(gm, test = "Chisq")[CompVariable, "Pr(>Chi)"]
        pu   <- suppressWarnings(stats::chisq.test(tbl0)$p.value)
        tl   <- "Logistic regression (LR)"

      } else {
        expct <- suppressWarnings(stats::chisq.test(tbl0)$expected)
        if (any(expct < 5)) {
          pu <- stats::fisher.test(tbl0, simulate.p.value = TRUE,
                                   B = 1e4)$p.value
          tl <- "Fisher (sim.)"
        } else {
          pu <- stats::chisq.test(tbl0)$p.value
          tl <- "Chi-squared"
        }
        pad <- NA_real_
      }
    }

    tibble::tibble(variable = var,
                   p_unadj = pu, p_adj = pad, test_label = tl)
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::mutate(
                                     p.value     = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value_fmt = gtsummary::style_pvalue(.data$p.value, digits = pDigits),
                                     p.value_fmt = ifelse(.data$row_type == "label", .data$p.value_fmt,
                                                          NA_character_),
                                     Test        = ifelse(.data$row_type == "label", .data$test_label,
                                                          NA_character_)
                                   )) %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows    = .data$p.value <= 0.05 & .data$row_type == "label",
      text_format = "bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  ## ── effect sizes ─────────────────────────────────────────────────────────
  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {
      df_var <- df[!is.na(df[[var]]), , drop = FALSE]
      df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
      k <- nlevels(df_var[[CompVariable]])
      n <- nrow(df_var)
      safe_var <- paste0("`", var, "`")

      if (n < 2 || k < 2)
        return(tibble::tibble(variable = var,
                              effect_size = NA_real_, es_method = NA_character_))

      if (is.numeric(df_var[[var]])) {

        if (Parametric && k == 2 && is.null(Covariates)) {
          fmla_es <- stats::as.formula(paste0(safe_var, " ~ ", safe_comp))
          es_val  <- abs(effectsize::cohens_d(fmla_es, df_var)$Cohens_d)
          method  <- "|d|"

        } else if (Parametric && !is.null(Covariates)) {
          cov_terms <- paste(c(safe_comp, Covariates), collapse = " + ")
          fmla_es   <- stats::as.formula(paste0(safe_var, " ~ ", cov_terms))
          aov_tab   <- car::Anova(stats::lm(fmla_es, df_var), type = 2)
          et        <- effectsize::eta_squared(aov_tab, partial = TRUE)
          idx       <- match(CompVariable, rownames(aov_tab))
          es_val    <- et$Eta2_partial[idx]
          method    <- "partial η²"

        } else if (Parametric) {
          fmla_es <- stats::as.formula(paste0(safe_var, " ~ ", safe_comp))
          et      <- effectsize::eta_squared(
            stats::aov(fmla_es, df_var), partial = FALSE)
          es_val  <- et$Eta2[1]; method <- "η²"

        } else {
          fmla_es <- stats::as.formula(paste0(safe_var, " ~ ", safe_comp))
          H       <- stats::kruskal.test(fmla_es, df_var)$statistic
          es_val  <- as.numeric((H - k + 1) / (n - k))
          method  <- "ε²"
        }

      } else {
        tab    <- table(df_var[[var]], df_var[[CompVariable]])
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
                                       effect_size = dplyr::if_else(.data$row_type == "label",
                                                                    .data$effect_size, NA_real_),
                                       ES_Method   = dplyr::if_else(.data$row_type == "label",
                                                                    .data$es_method,   NA_character_)
                                     )) %>%
      gtsummary::modify_fmt_fun(effect_size ~ function(x)
        ifelse(is.na(x), NA_character_,
               formatC(x, digits = EffectSizeDigits, format = "f"))) %>%
      gtsummary::modify_header(effect_size ~ "**Effect Size**") %>%
      gtsummary::modify_header(ES_Method   ~ "**ES Method**")
  }

  ## ── pairwise contrasts ──────────────────────────────────────────────────
  if (AddPairwise && nlevels(df[[CompVariable]]) > 1) {

    lvls   <- levels(df[[CompVariable]])
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% lvls)
        stop("Referent level not found: ", Referent)
      lapply(setdiff(lvls, Referent), function(x) c(Referent, x))
    } else utils::combn(lvls, 2, simplify = FALSE)

    pw <- purrr::map_dfr(Variables, function(var) {
      safe_var <- paste0("`", var, "`")

      if (is.numeric(df[[var]])) {                         # continuous
        if (Parametric) {
          fit <- stats::aov(
            stats::as.formula(paste0(safe_var, " ~ ", safe_comp)), df)
          em  <- emmeans::emmeans(fit, CompVariable)
          ct  <- if (!is.null(Referent))
            emmeans::contrast(em, "trt.vs.ctrl",
                              ref = match(Referent, lvls),
                              adjust = PairwiseMethod)
          else emmeans::contrast(em, "pairwise",
                                 adjust = PairwiseMethod)

          r <- broom::tidy(ct)
          if (!"adj.p.value" %in% names(r))           # ensure column exists
            r$adj.p.value <- NA_real_
          r <- r %>%
            dplyr::transmute(variable = var,
                             contrast = .canonical_label(contrast),
                             p_val    = dplyr::coalesce(adj.p.value,
                                                        p.value))

        } else {                                          # Wilcoxon
          r <- purrr::map_dfr(combos, function(cp) {
            idx <- df[[CompVariable]] %in% cp
            res <- stats::pairwise.wilcox.test(
              x = df[[var]][idx],
              g = droplevels(df[[CompVariable]][idx]),
              p.adjust.method = PairwiseMethod
            )$p.value
            as.data.frame(as.table(res)) %>%
              dplyr::transmute(
                variable = var,
                contrast = .canonical_label(.pair_label(Var1, Var2)),
                p_val    = Freq
              )
          })
        }

      } else {                                            # categorical
        r <- purrr::map_dfr(combos, function(cp) {
          sub   <- df[df[[CompVariable]] %in% cp, ]
          tbl0  <- table(droplevels(as.factor(sub[[var]])),
                         droplevels(sub[[CompVariable]]))
          pv    <- if (nrow(tbl0) < 2 || ncol(tbl0) < 2) NA_real_
          else if (any(
            suppressWarnings(stats::chisq.test(tbl0)$expected) < 5))
            stats::fisher.test(tbl0, simulate.p.value = TRUE,
                               B = 1e4)$p.value
          else stats::chisq.test(tbl0)$p.value
          tibble::tibble(
            variable = var,
            contrast = .canonical_label(.pair_label(cp[1], cp[2])),
            p_val    = pv
          )
        }) %>%
          dplyr::mutate(
            p_val = if (PairwiseMethod != "none")
              p.adjust(p_val, PairwiseMethod) else p_val)
      }
    })

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
                                       dplyr::across(tidyselect::all_of(
                                         setdiff(names(pw_wide), "variable")),
                                         ~ ifelse(.data$row_type == "label", ., NA_real_))
                                     ))

    for (col in setdiff(names(pw_wide), "variable")) {
      tbl <- tbl %>%
        gtsummary::modify_fmt_fun(!!rlang::sym(col) ~
                                    function(x) gtsummary::style_pvalue(x, digits = pDigits)) %>%
        gtsummary::modify_table_styling(
          columns     = col,
          rows        = .data[[col]] <= 0.05 & .data$row_type == "label",
          text_format = "bold"
        ) %>%
        gtsummary::modify_header(
          !!rlang::sym(col) := paste0("**", col, "**")) %>%
        gtsummary::modify_footnote(
          !!rlang::sym(col) ~ paste0("Pairwise p-value (", PairwiseMethod,
                                     if (!is.null(Referent)) paste0("; vs ", Referent), ")"))
    }
  }

  ## ── caption & return ────────────────────────────────────────────────────
  cap <- sprintf(
    "Comparison Table (%s analysis%s) — see **Test** column for method used per variable",
    if (Parametric) "parametric" else "non-parametric",
    if (!is.null(Covariates) && Parametric)
      paste0("; adjusted for ", paste(Covariates, collapse = ", "))
    else "")
  tbl %>% gtsummary::modify_caption(cap)
}


