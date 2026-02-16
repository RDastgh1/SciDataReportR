#' Make a publication-ready comparison table using gtsummary
#'
#' `MakeComparisonTable()` is a label-friendly wrapper around **gtsummary** that produces a
#' publication-ready comparison table (`tbl_summary`) and optionally adds:
#' covariate-adjusted global tests, effect sizes, and pairwise contrasts (with p-value adjustment).
#'
#' The function is designed to be **predictable** and **deterministic**:
#' outputs are stable, p-values are clamped to avoid numeric underflow to 0,
#' and complete-case filtering is explicit and can drop group levels (reported via Notes).
#'
#' ## Quick start (defaults)
#' Only the data and grouping variable are required:
#' \preformatted{
#' MakeComparisonTable(df, "Cluster")
#' }
#'
#' - If `Variables` is `NULL` (default), the function uses **all columns** in `DataFrame`
#'   except the grouping variable (`CompVariable`) and any covariates (`Covariates`).
#' - All other arguments have defaults and can be overridden as needed.
#'
#' ## Defaults and when to change them
#' - `Variables = NULL`: include all non-group, non-covariate columns. Provide `Variables`
#'   when you want a curated set or a specific order.
#' - `Covariates = NULL`: unadjusted tests only. Provide covariates to run adjusted global tests and
#'   adjusted pairwise contrasts.
#' - `Parametric = TRUE`: use parametric tests (Welch t-test or ANOVA; ANCOVA with Type II test when adjusted).
#'   Set `Parametric = FALSE` to use nonparametric tests when unadjusted (Wilcoxon or Kruskal),
#'   and **robust ANCOVA** (HC3) when adjusted.
#' - `CatMethod = "auto"`: use chi-squared when assumptions are acceptable; fall back to Fisher when they are not.
#'   Set `"chisq"` or `"fisher"` to force a method.
#' - `MultiCatAdjusted = "multinomial_LR"`: required default for adjusted multi-category outcomes
#'   (prevents “silent breakage” when users forget to set it).
#' - `AddPairwise = FALSE`: set `TRUE` to add pairwise p-value columns.
#' - `PairwiseMethod = "bonferroni"`: change to `"holm"`, `"BH"`/`"fdr"`, etc, or `"none"`.
#' - `Referent = NULL`: when set (e.g., `"0"`), pairwise contrasts become treatment-vs-control against the referent.
#' - `AddEffectSize = FALSE`: set `TRUE` to add an `effect_size` column using **effectsize**.
#' - `IncludeNotes = TRUE`: set `FALSE` to drop Notes entirely. When included, Notes is always the **last** column.
#'
#' ## Statistical methods used
#'
#' ### Continuous outcomes (unadjusted)
#' If `Covariates` is `NULL`:
#' - `Parametric = TRUE`:
#'   - 2 groups: Welch t-test (`stats::t.test`, var.equal = FALSE)
#'   - 3+ groups: one-way ANOVA (`stats::aov`)
#' - `Parametric = FALSE`:
#'   - 2 groups: Wilcoxon rank-sum (`stats::wilcox.test`)
#'   - 3+ groups: Kruskal-Wallis (`stats::kruskal.test`)
#'
#' Pairwise (if `AddPairwise = TRUE`):
#' - `Parametric = TRUE`: `stats::pairwise.t.test(pool.sd = FALSE)`
#' - `Parametric = FALSE`: `stats::pairwise.wilcox.test`
#' - p-adjust via `PairwiseMethod` (or none when `"none"`).
#'
#' ### Continuous outcomes (adjusted)
#' If `Covariates` is provided, complete-case filtering is performed on outcome + group + covariates:
#' - Global adjusted test:
#'   - `Parametric = TRUE`: ANCOVA via `lm`, Type II test of the group term via `car::Anova(type = 2)`
#'   - `Parametric = FALSE`: Robust ANCOVA via `lm` + HC3 covariance (`sandwich::vcovHC(type="HC3")`)
#'     passed to `car::Anova(vcov.=..., test.statistic="F")`
#' - Adjusted pairwise:
#'   - Uses `emmeans::emmeans()` on the adjusted model
#'   - In robust mode, HC3 covariance is passed through `emmeans(..., vcov. = V)` so pairwise does not disappear.
#'
#' ### Categorical outcomes (unadjusted global)
#' Uses a contingency table between outcome and group:
#' - `CatMethod = "chisq"`: Pearson chi-squared (`stats::chisq.test`, no Yates correction)
#' - `CatMethod = "fisher"`: Fisher exact for 2x2; Fisher with simulation for larger tables
#' - `CatMethod = "auto"` (default): attempt chi-squared; if expected counts are too small, fall back to Fisher.
#'
#' ### Categorical outcomes (adjusted global)
#' Complete-case filtering is performed on outcome + group + covariates:
#' - Binary outcome: logistic regression likelihood ratio (LR) test of the group term
#'   (`glm(family=binomial)` + `drop1(test="Chisq")`)
#' - Multi-category outcome (3+ levels): multinomial LR test by comparing reduced vs full model
#'   (`nnet::multinom` + `anova(test="Chisq")` with a logLik LR fallback)
#'
#' ### Categorical outcomes (adjusted pairwise)
#' For each pairwise subset (or referent-vs-each):
#' - If the outcome collapses to 2 levels in that subset: logistic LR
#' - If still 3+ levels: multinomial LR
#'
#' This per-pair collapse detection prevents missing pairwise results for variables like Race.
#'
#' ## Notes column and dropped group levels
#' Adjusted analyses use complete-case filtering on outcome + group + covariates.
#' This can drop entire group levels for a given outcome (e.g. Cluster "3" has zero complete cases).
#' Pairwise contrasts involving dropped levels are returned as `NA`, and (if `IncludeNotes = TRUE`)
#' the Notes column explains which levels were dropped.
#'
#' ## P-values
#' Extremely small p-values can underflow to numeric 0. All numeric p-values are clamped to at least
#' `.Machine$double.xmin`, while formatted p-values still display nicely (e.g. "<0.001").
#'
#' ## References (packages and methods)
#' - gtsummary: https://www.danieldsjoberg.com/gtsummary/
#' - car (Type II tests): https://cran.r-project.org/package=car
#' - emmeans: https://cran.r-project.org/package=emmeans
#' - sandwich (HC3 covariance): https://cran.r-project.org/package=sandwich
#' - nnet (multinomial): https://cran.r-project.org/package=nnet
#' - effectsize: https://cran.r-project.org/package=effectsize
#'
#' @param DataFrame A data frame containing the grouping variable and variables to summarize.
#' @param CompVariable A single character string giving the grouping variable name (required for group comparisons).
#' @param Variables Character vector of variables to include. Default `NULL` uses all columns except the grouping
#'   variable and covariates.
#' @param ... Unnamed additional variable names (character). These are appended to `Variables`.
#' @param Covariates Optional character vector of covariate names for adjusted models (default `NULL`).
#' @param ValueDigits Digits for continuous summary values (default `2`).
#' @param pDigits Digits for p-value formatting (default `3`).
#' @param AddEffectSize If `TRUE`, adds an `effect_size` column (default `FALSE`).
#' @param EffectSizeDigits Digits for effect size formatting (default `2`).
#' @param AddPairwise If `TRUE`, adds pairwise contrast p-value columns (default `FALSE`).
#' @param PairwiseMethod P-value adjustment method for pairwise tests. One of `"none"` or `stats::p.adjust.methods`
#'   (default `"bonferroni"`).
#' @param Parametric If `TRUE`, use parametric tests; if `FALSE`, use nonparametric (unadjusted) and robust (adjusted)
#'   tests (default `TRUE`).
#' @param ParametricDisplay If `NULL` (default), display matches `Parametric`. If `TRUE`, show mean (SD). If `FALSE`,
#'   show median [IQR].
#' @param IncludeOverallN If `TRUE`, adds N column via `gtsummary::add_n()` (default `FALSE`).
#' @param IncludeMissing If `TRUE`, include missing rows via `missing="ifany"` (default `FALSE`).
#' @param suppress_warnings If `TRUE`, suppress warnings from `gtsummary::tbl_summary()` (default `FALSE`).
#' @param Referent Optional referent level for treatment-vs-control contrasts (default `NULL`).
#' @param IncludeOverallStats If `TRUE`, produce an overall-only summary (no grouping) and return early (default `FALSE`).
#' @param ShowPositiveBinaryOnLabel If `TRUE`, shows the "positive" level for dichotomous variables when possible (default `TRUE`).
#' @param BinaryPairwiseScale Scale hint for binary effect sizes: `"logit"` or `"prob"` (default `"logit"`).
#' @param CatMethod Categorical global test method: `"auto"` (default), `"chisq"`, or `"fisher"`.
#' @param MultiCatAdjusted Adjusted method for 3+ level outcomes. Default `"multinomial_LR"`.
#' @param IncludeNotes If `TRUE`, include Notes column (default `TRUE`). If included, Notes is always the last column.
#'
#' @return A `gtsummary::tbl_summary` object with added columns in `table_body`:
#' - `p.value` (numeric clamped) and `p.value_fmt` (formatted) on the label row
#' - `Test` (test label) on the label row
#' - optional `effect_size` on the label row
#' - optional `pw_*` pairwise columns on the label row
#' - optional `Notes` on the label row (always last when included)
#'
#' @examples
#' \dontrun{
#' # Minimal: just data + grouping variable (Variables inferred)
#' tbl <- MakeComparisonTable(df, "Cluster")
#'
#' # Adjusted robust with pairwise and effect sizes
#' tbl2 <- MakeComparisonTable(
#'   df, "Cluster",
#'   Covariates = c("Age", "Sex", "wrat3"),
#'   Parametric = FALSE,
#'   AddPairwise = TRUE,
#'   PairwiseMethod = "holm",
#'   AddEffectSize = TRUE
#' )
#'
#' # Referent contrasts
#' tbl3 <- MakeComparisonTable(
#'   df, "Cluster",
#'   Covariates = c("Age", "Sex", "wrat3"),
#'   AddPairwise = TRUE,
#'   Referent = "0"
#' )
#' }
#' @export
MakeComparisonTable <- function(
    DataFrame,
    CompVariable,
    Variables = NULL,
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
    BinaryPairwiseScale = c("logit", "prob"),
    CatMethod = c("auto", "chisq", "fisher"),
    MultiCatAdjusted = c("multinomial_LR"),
    IncludeNotes = TRUE
) {

  # Validate inputs -----------------------------------------------------------
  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "sandwich"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  if (isTRUE(AddEffectSize) && !requireNamespace("effectsize", quietly = TRUE)) {
    stop("Please install: effectsize (required when AddEffectSize = TRUE).")
  }
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Please install: nnet (required for multinomial adjusted tests).")
  }

  if (!is.data.frame(DataFrame)) stop("DataFrame must be a data frame.")
  if (!is.character(CompVariable) || length(CompVariable) != 1) {
    stop("CompVariable must be a single character string.")
  }
  if (!(CompVariable %in% names(DataFrame))) stop("Grouping variable not found: ", CompVariable)

  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale)
  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)

  valid_methods <- c("none", stats::p.adjust.methods)
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1 || !(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "),
         ". Got: ", paste(PairwiseMethod, collapse = ", "))
  }

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  if (!is.null(Covariates)) {
    if (!is.character(Covariates)) stop("Covariates must be a character vector or NULL.")
    if (!all(Covariates %in% names(DataFrame))) {
      stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))
    }
  }

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  # Prepare variables (defaults + ...) --------------------------------------
  extra_vars <- list(...)
  if (length(extra_vars) > 0) {
    extra_vars <- unlist(extra_vars, use.names = FALSE)
    if (!is.character(extra_vars)) stop("Unnamed arguments in ... must be character variable names.")
  } else {
    extra_vars <- character(0)
  }

  if (is.null(Variables)) {
    Variables <- setdiff(names(DataFrame), c(CompVariable, Covariates))
    Variables <- c(Variables, extra_vars)
  } else {
    if (!is.character(Variables)) stop("Variables must be a character vector or NULL.")
    Variables <- c(Variables, extra_vars)
  }
  Variables <- unique(Variables)

  if (!length(Variables)) stop("No Variables to summarise (after defaults and exclusions).")
  if (!all(Variables %in% names(DataFrame))) {
    stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  }

  # Helpers -----------------------------------------------------------------
  .as_factor <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  .norm <- function(x) gsub("\\s+", " ", trimws(as.character(x)))

  .pair_label <- function(a, b) paste(a, b, sep = " - ")

  .pair_key <- function(a, b) {
    a <- .norm(a); b <- .norm(b)
    paste(sort(c(a, b)), collapse = "||")
  }

  .btick <- function(x) {
    x <- gsub("`", "\\\\`", as.character(x))
    paste0("`", x, "`")
  }

  .fmla <- function(lhs, rhs_terms) {
    rhs_terms <- as.character(rhs_terms)
    rhs <- paste(.btick(rhs_terms), collapse = " + ")
    stats::as.formula(paste(.btick(lhs), "~", rhs))
  }

  # Vectorized clamp to avoid numeric 0 and scalar-if bugs
  .clamp_p <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    out[ok] <- pmax(p[ok], .Machine$double.xmin)
    out
  }

  .fmt_p <- function(p, digits = 3) {
    p <- suppressWarnings(as.numeric(p))
    out <- rep(NA_character_, length(p))
    ok <- is.finite(p)
    out[ok] <- format.pval(p[ok], digits = digits, eps = 10^-digits)
    out
  }

  .clean_tab <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    tab
  }

  .cat_global_test <- function(tab, method = c("auto", "chisq", "fisher")) {
    method <- match.arg(method)
    tab <- .clean_tab(tab)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(list(p = NA_real_, label = "Insufficient data"))

    if (method == "chisq") {
      chi <- tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)), error = function(e) NULL)
      if (!is.null(chi) && !is.null(chi$p.value)) return(list(p = .clamp_p(chi$p.value), label = "Pearson chi-squared"))
      return(list(p = NA_real_, label = "Chi-squared failed"))
    }

    if (method == "fisher") {
      fish <- tryCatch(
        if (nrow(tab) == 2 && ncol(tab) == 2) stats::fisher.test(tab)
        else stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e4),
        error = function(e) NULL
      )
      if (!is.null(fish) && !is.null(fish$p.value)) {
        return(list(
          p = .clamp_p(fish$p.value),
          label = if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"
        ))
      }
      return(list(p = NA_real_, label = "Fisher failed"))
    }

    # auto
    chi_obj <- tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)), error = function(e) NULL)
    if (is.null(chi_obj) || is.null(chi_obj$expected)) return(.cat_global_test(tab, method = "fisher"))
    if (any(chi_obj$expected < 5)) return(.cat_global_test(tab, method = "fisher"))
    list(p = .clamp_p(chi_obj$p.value), label = "Pearson chi-squared")
  }

  .robust_type2_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)
    a2 <- tryCatch(car::Anova(fit_lm, type = 2, vcov. = V, test.statistic = "F"),
                   error = function(e) NULL)
    if (is.null(a2) || !(term %in% rownames(a2))) return(NA_real_)
    .clamp_p(a2[term, "Pr(>F)"])
  }

  .cc_levels_note <- function(full_levels, kept_levels) {
    dropped <- setdiff(full_levels, kept_levels)
    if (!length(dropped)) return(NA_character_)
    paste0(
      "Adjusted pairwise unavailable for contrasts involving dropped group level(s): ",
      paste(dropped, collapse = ", "),
      " (complete-case filtering on outcome and covariates)."
    )
  }

  .strip_group_prefix <- function(x, compvar) {
    x <- .norm(x)
    x <- gsub(paste0("^", compvar, "\\s*=?\\s*"), "", x)
    x <- gsub(paste0("^", compvar), "", x)
    .norm(x)
  }

  .extract_anova_lr_p <- function(a) {
    if (is.null(a) || nrow(a) < 2) return(NA_real_)
    cn <- colnames(a)
    pcol <- intersect(c("Pr(>Chi)", "Pr(Chi)", "Pr(>Chisq)", "Pr(Chisq)"), cn)
    if (!length(pcol)) return(NA_real_)
    as.numeric(a[2, pcol[1]])
  }

  .lrt_loglik_p <- function(m_red, m_full) {
    ll0 <- tryCatch(as.numeric(stats::logLik(m_red)),  error = function(e) NA_real_)
    ll1 <- tryCatch(as.numeric(stats::logLik(m_full)), error = function(e) NA_real_)
    df0 <- tryCatch(attr(stats::logLik(m_red),  "df"), error = function(e) NA_integer_)
    df1 <- tryCatch(attr(stats::logLik(m_full), "df"), error = function(e) NA_integer_)
    if (anyNA(c(ll0, ll1, df0, df1))) return(NA_real_)
    LR <- 2 * (ll1 - ll0)
    ddf <- df1 - df0
    if (is.na(LR) || is.na(ddf) || ddf <= 0) return(NA_real_)
    .clamp_p(stats::pchisq(LR, df = ddf, lower.tail = FALSE))
  }

  .multinom_global_lr <- function(df_cc, var, compvar, covs) {
    df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))
    df_cc[[compvar]] <- droplevels(.as_factor(df_cc[[compvar]]))
    if (nlevels(df_cc[[var]]) < 3 || nlevels(df_cc[[compvar]]) < 2) {
      return(list(p = NA_real_, label = "Insufficient data (multinomial)"))
    }

    f_full <- .fmla(var, c(compvar, covs))
    f_red  <- .fmla(var, covs)

    m_full <- tryCatch(nnet::multinom(f_full, data = df_cc, trace = FALSE), error = function(e) NULL)
    m_red  <- tryCatch(nnet::multinom(f_red,  data = df_cc, trace = FALSE), error = function(e) NULL)
    if (is.null(m_full) || is.null(m_red)) return(list(p = NA_real_, label = "Multinomial failed"))

    a <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
    p <- .extract_anova_lr_p(a)
    if (is.na(p)) p <- .lrt_loglik_p(m_red, m_full)
    if (is.na(p)) return(list(p = NA_real_, label = "Multinomial LR failed"))

    list(p = .clamp_p(p), label = "Multinomial LR")
  }

  .cat_pairwise_lr <- function(df_cc, var, compvar, covs, combos, PairwiseMethod) {

    out <- purrr::map_dfr(combos, function(cp) {

      sub <- df_cc[df_cc[[compvar]] %in% cp, , drop = FALSE]
      sub[[compvar]] <- droplevels(.as_factor(sub[[compvar]]))
      sub[[var]] <- droplevels(.as_factor(sub[[var]]))

      if (nlevels(sub[[compvar]]) < 2 || nrow(sub) < 5 || nlevels(sub[[var]]) < 2) {
        return(tibble::tibble(
          contrast_label = .pair_label(cp[1], cp[2]),
          key = .pair_key(cp[1], cp[2]),
          p_val = NA_real_
        ))
      }

      # binary logistic LR
      if (nlevels(sub[[var]]) == 2) {
        f_full <- .fmla(var, c(compvar, covs))
        gm <- tryCatch(stats::glm(f_full, data = sub, family = stats::binomial()), error = function(e) NULL)
        if (is.null(gm)) {
          return(tibble::tibble(
            contrast_label = .pair_label(cp[1], cp[2]),
            key = .pair_key(cp[1], cp[2]),
            p_val = NA_real_
          ))
        }
        d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
        p <- NA_real_
        if (!is.null(d1) && compvar %in% rownames(d1) && "Pr(>Chi)" %in% colnames(d1)) {
          p <- as.numeric(d1[compvar, "Pr(>Chi)"])
        }
        return(tibble::tibble(
          contrast_label = .pair_label(cp[1], cp[2]),
          key = .pair_key(cp[1], cp[2]),
          p_val = .clamp_p(p)
        ))
      }

      # multinomial LR
      if (nlevels(sub[[var]]) >= 3) {
        f_full <- .fmla(var, c(compvar, covs))
        f_red  <- .fmla(var, covs)

        m_full <- tryCatch(nnet::multinom(f_full, data = sub, trace = FALSE), error = function(e) NULL)
        m_red  <- tryCatch(nnet::multinom(f_red,  data = sub, trace = FALSE), error = function(e) NULL)

        if (is.null(m_full) || is.null(m_red)) {
          return(tibble::tibble(
            contrast_label = .pair_label(cp[1], cp[2]),
            key = .pair_key(cp[1], cp[2]),
            p_val = NA_real_
          ))
        }

        a <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
        p <- .extract_anova_lr_p(a)
        if (is.na(p)) p <- .lrt_loglik_p(m_red, m_full)

        return(tibble::tibble(
          contrast_label = .pair_label(cp[1], cp[2]),
          key = .pair_key(cp[1], cp[2]),
          p_val = .clamp_p(p)
        ))
      }

      tibble::tibble(
        contrast_label = .pair_label(cp[1], cp[2]),
        key = .pair_key(cp[1], cp[2]),
        p_val = NA_real_
      )
    })

    if (!identical(PairwiseMethod, "none")) {
      out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    }
    out
  }

  .effect_size_str <- function(var, df_vg, is_cont, full_group_levels) {
    if (!isTRUE(AddEffectSize)) return(NA_character_)

    # Continuous effect sizes
    if (isTRUE(is_cont)) {
      k <- nlevels(df_vg[[CompVariable]])
      if (k < 2) return(NA_character_)

      # Adjusted (covariates) -> partial eta squared for group term
      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df_vg[stats::complete.cases(df_vg[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) return(NA_character_)

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) return(NA_character_)

        es <- tryCatch(effectsize::eta_squared(fit, partial = TRUE), error = function(e) NULL)
        if (is.null(es) || !("Parameter" %in% names(es)) || !("Eta2_partial" %in% names(es))) return(NA_character_)
        row <- es[es$Parameter == CompVariable, , drop = FALSE]
        if (!nrow(row)) return(NA_character_)
        val <- row$Eta2_partial[1]
        if (!is.finite(val)) return(NA_character_)
        return(paste0("η²p=", formatC(val, digits = EffectSizeDigits, format = "f")))
      }

      # Unadjusted
      if (Parametric) {
        if (k == 2) {
          d <- tryCatch(effectsize::cohens_d(.fmla(var, CompVariable), data = df_vg, hedges.correction = TRUE),
                        error = function(e) NULL)
          if (is.null(d) || !("Cohens_d" %in% names(d))) return(NA_character_)
          val <- d$Cohens_d[1]
          if (!is.finite(val)) return(NA_character_)
          return(paste0("d=", formatC(val, digits = EffectSizeDigits, format = "f")))
        }

        aov_obj <- tryCatch(stats::aov(.fmla(var, CompVariable), data = df_vg), error = function(e) NULL)
        if (is.null(aov_obj)) return(NA_character_)
        es <- tryCatch(effectsize::eta_squared(aov_obj, partial = FALSE), error = function(e) NULL)
        if (is.null(es) || !("Parameter" %in% names(es)) || !("Eta2" %in% names(es))) return(NA_character_)
        row <- es[es$Parameter == CompVariable, , drop = FALSE]
        if (!nrow(row)) return(NA_character_)
        val <- row$Eta2[1]
        if (!is.finite(val)) return(NA_character_)
        return(paste0("η²=", formatC(val, digits = EffectSizeDigits, format = "f")))
      }

      # Nonparametric unadjusted
      if (k == 2) {
        rb <- tryCatch(effectsize::rank_biserial(.fmla(var, CompVariable), data = df_vg),
                       error = function(e) NULL)
        if (is.null(rb) || !("r_rank_biserial" %in% names(rb))) return(NA_character_)
        val <- rb$r_rank_biserial[1]
        if (!is.finite(val)) return(NA_character_)
        return(paste0("rrb=", formatC(val, digits = EffectSizeDigits, format = "f")))
      } else {
        kw <- tryCatch(stats::kruskal.test(.fmla(var, CompVariable), data = df_vg),
                       error = function(e) NULL)
        if (is.null(kw)) return(NA_character_)
        es <- tryCatch(effectsize::epsilon_squared(kw), error = function(e) NULL)
        if (is.null(es) || !("Epsilon2" %in% names(es))) return(NA_character_)
        val <- es$Epsilon2[1]
        if (!is.finite(val)) return(NA_character_)
        return(paste0("ε²=", formatC(val, digits = EffectSizeDigits, format = "f")))
      }
    }

    # Categorical effect sizes -> Cramer's V from var x group contingency
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)
    tab <- .clean_tab(tab)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_character_)

    v <- tryCatch(effectsize::cramers_v(tab, ci = FALSE), error = function(e) NULL)
    if (is.null(v)) return(NA_character_)
    # effectsize output can vary; handle common names
    val <- NA_real_
    if ("Cramers_v" %in% names(v)) val <- v$Cramers_v[1]
    if (!is.finite(val) && "Cramers_V" %in% names(v)) val <- v$Cramers_V[1]
    if (!is.finite(val)) return(NA_character_)

    paste0("V=", formatC(val, digits = EffectSizeDigits, format = "f"))
  }

  # Exclude covariates from Variables ---------------------------------------
  if (!is.null(Covariates)) {
    drop_cov <- intersect(Variables, Covariates)
    if (length(drop_cov)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov, collapse = ", "))
      Variables <- setdiff(Variables, drop_cov)
    }
  }

  # Subset and coerce group --------------------------------------------------
  cols <- c(CompVariable, Variables, Covariates)
  cols <- unique(cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  df[[CompVariable]] <- .as_factor(df[[CompVariable]])

  # Drop constants -----------------------------------------------------------
  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]
  drop <- setdiff(Variables, keep)
  if (length(drop)) warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  Variables <- keep
  if (!length(Variables)) stop("No variables left to summarise after dropping constants.")

  # Type inference (gtsummary) ----------------------------------------------
  n_unique <- vapply(Variables, function(v) length(unique(df[[v]][!is.na(df[[v]])])), integer(1))
  is_num   <- vapply(Variables, function(v) is.numeric(df[[v]]), logical(1))
  is_dich  <- n_unique == 2
  treat_as_continuous <- is_num & !is_dich

  type_list <- NULL
  numeric_cont <- Variables[treat_as_continuous]
  dichotomous_numeric <- Variables[is_num & is_dich]
  if (length(numeric_cont) || length(dichotomous_numeric)) {
    type_list <- c(
      rlang::set_names(as.list(rep("continuous",  length(numeric_cont))),        numeric_cont),
      rlang::set_names(as.list(rep("dichotomous", length(dichotomous_numeric))), dichotomous_numeric)
    )
  }

  # gtsummary display control -----------------------------------------------
  value_list <- NULL
  if (isTRUE(ShowPositiveBinaryOnLabel)) {
    pos_tokens <- c("TRUE", "1", "YES", "Yes")
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
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit)) vmap[[v]] <- hit[1]
      } else {
        levs <- sort(unique(as.character(ux)))
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit)) vmap[[v]] <- hit[1]
      }
    }
    if (length(vmap)) value_list <- vmap
  }

  stat_cont <- if (isTRUE(ParametricDisplay)) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  # Caption builder ----------------------------------------------------------
  .caption_text <- function() {
    cont_disp <- if (isTRUE(ParametricDisplay)) "Continuous summaries: mean (SD)" else "Continuous summaries: median [IQR]"
    test_mode <- if (isTRUE(Parametric)) "Parametric tests" else "Robust nonparametric mode (HC3 for adjusted continuous; Wilcoxon/Kruskal unadjusted)"
    cov_txt <- if (!is.null(Covariates) && length(Covariates)) paste0("Adjusted models include covariates: ", paste(Covariates, collapse = ", "), ".") else "No covariate adjustment."
    cat_txt <- paste0("Categorical global test: ", CatMethod, " (auto uses Fisher only when chi-squared assumptions fail).")
    multi_txt <- paste0("Multi-category adjusted: ", MultiCatAdjusted, ".")
    pw_txt <- if (isTRUE(AddPairwise)) paste0("Pairwise contrasts included; p-adjust: ", PairwiseMethod, ".") else "No pairwise contrasts."
    paste(cont_disp, test_mode, cov_txt, cat_txt, multi_txt, pw_txt, sep = " ")
  }

  # Overall-only mode --------------------------------------------------------
  if (isTRUE(IncludeOverallStats)) {
    tbl0 <- gtsummary::tbl_summary(
      df,
      include   = tidyselect::all_of(Variables),
      missing   = if (isTRUE(IncludeMissing)) "ifany" else "no",
      statistic = list(
        gtsummary::all_continuous()  ~ stat_cont,
        gtsummary::all_categorical() ~ stat_cat
      ),
      digits    = list(gtsummary::all_continuous() ~ ValueDigits),
      type      = type_list,
      value     = value_list
    )
    if (isTRUE(IncludeOverallN)) tbl0 <- tbl0 %>% gtsummary::add_n()
    if (isTRUE(suppress_warnings)) tbl0 <- suppressWarnings(tbl0)
    tbl0 <- tbl0 %>% gtsummary::modify_caption(.caption_text())
    return(tbl0)
  }

  # Grouped tbl_summary ------------------------------------------------------
  tbl <- gtsummary::tbl_summary(
    df,
    by        = CompVariable,
    include   = tidyselect::all_of(Variables),
    missing   = if (isTRUE(IncludeMissing)) "ifany" else "no",
    statistic = list(
      gtsummary::all_continuous()  ~ stat_cont,
      gtsummary::all_categorical() ~ stat_cat
    ),
    digits    = list(gtsummary::all_continuous() ~ ValueDigits),
    type      = type_list,
    value     = value_list
  )
  if (isTRUE(IncludeOverallN)) tbl <- tbl %>% gtsummary::add_n()
  if (isTRUE(suppress_warnings)) tbl <- suppressWarnings(tbl)
  if (nlevels(df[[CompVariable]]) < 2) {
    tbl <- tbl %>% gtsummary::modify_caption(.caption_text())
    return(tbl)
  }

  full_group_levels <- levels(droplevels(.as_factor(df[[CompVariable]])))

  # Global p-values (+ adjusted) + Notes + effect sizes ----------------------
  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    # complete cases for unadjusted global tests
    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))

    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(
        variable = var, p_unadj = NA_real_, p_adj = NA_real_,
        test_label = "Insufficient groups",
        Notes = NA_character_,
        effect_size = .effect_size_str(var, df_vg, is_cont, full_group_levels)
      ))
    }

    note <- NA_character_
    p_adj <- NA_real_
    test_label <- NA_character_
    p_un <- NA_real_

    # Continuous -------------------------------------------------------------
    if (isTRUE(is_cont)) {

      # Adjusted continuous
      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))

        note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Insufficient data (adjusted)",
            Notes = note,
            effect_size = .effect_size_str(var, df_cc, is_cont, full_group_levels)
          ))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Model failed (adjusted)",
            Notes = note,
            effect_size = .effect_size_str(var, df_cc, is_cont, full_group_levels)
          ))
        }

        # “Unadjusted” p within same CC set for reporting consistency
        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        p_un <- .clamp_p(p_un)

        if (isTRUE(Parametric)) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) .clamp_p(a2[CompVariable, "Pr(>F)"]) else NA_real_
          test_label <- "ANCOVA (Type II)"
        } else {
          p_adj <- .robust_type2_p(fit, CompVariable)
          test_label <- "Robust ANCOVA (HC3 Type II)"
        }

        return(tibble::tibble(
          variable = var, p_unadj = p_un, p_adj = p_adj,
          test_label = test_label,
          Notes = note,
          effect_size = .effect_size_str(var, df_cc, is_cont, full_group_levels)
        ))
      }

      # Unadjusted continuous
      k <- nlevels(df_vg[[CompVariable]])
      if (isTRUE(Parametric)) {
        if (k == 2) {
          p_un <- tryCatch(stats::t.test(.fmla(var, CompVariable), data = df_vg)$p.value,
                           error = function(e) NA_real_)
          test_label <- "Welch t-test"
        } else {
          p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
                           error = function(e) NA_real_)
          test_label <- "ANOVA"
        }
      } else {
        p_un <- tryCatch(
          if (k == 2) stats::wilcox.test(.fmla(var, CompVariable), data = df_vg)$p.value
          else stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$p.value,
          error = function(e) NA_real_
        )
        test_label <- if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"
      }

      return(tibble::tibble(
        variable = var, p_unadj = .clamp_p(p_un), p_adj = NA_real_,
        test_label = test_label,
        Notes = NA_character_,
        effect_size = .effect_size_str(var, df_vg, is_cont, full_group_levels)
      ))
    }

    # Categorical ------------------------------------------------------------
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- .cat_global_test(tab, method = CatMethod)
    p_un <- tst$p
    test_label <- tst$label

    # Adjusted categorical if covariates provided
    if (!is.null(Covariates)) {
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))
      note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

      if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
        df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
      }

      # binary: logistic LR
      if (nlevels(df_cc[[var]]) == 2 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5) {
        gm <- tryCatch(stats::glm(.fmla(var, c(CompVariable, Covariates)),
                                  data = df_cc, family = stats::binomial()),
                       error = function(e) NULL)
        if (!is.null(gm)) {
          d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
          if (!is.null(d1) && CompVariable %in% rownames(d1) && "Pr(>Chi)" %in% colnames(d1)) {
            p_adj <- .clamp_p(d1[CompVariable, "Pr(>Chi)"])
            test_label <- "Logistic regression (LR)"
          }
        }
      }

      # multi-category: multinomial LR
      if (nlevels(df_cc[[var]]) >= 3 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 10) {
        if (MultiCatAdjusted == "multinomial_LR") {
          g_lr <- .multinom_global_lr(df_cc, var, CompVariable, Covariates)
          if (!is.na(g_lr$p)) {
            p_adj <- .clamp_p(g_lr$p)
            test_label <- g_lr$label
          }
        }
      }
    }

    tibble::tibble(
      variable = var,
      p_unadj = p_un,
      p_adj = p_adj,
      test_label = test_label,
      Notes = note,
      effect_size = .effect_size_str(var, df_vg, is_cont, full_group_levels)
    )
  })

  # Merge into gtsummary table_body -----------------------------------------
  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::group_by(.data$variable) %>%
                                   dplyr::mutate(.is_main_row = dplyr::row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(
                                     p.value = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value = dplyr::if_else(.data$.is_main_row, .clamp_p(.data$p.value), NA_real_),
                                     p.value_fmt = dplyr::if_else(.data$.is_main_row, .fmt_p(.data$p.value, digits = pDigits), NA_character_),
                                     Test = dplyr::if_else(.data$.is_main_row, .data$test_label, NA_character_),
                                     effect_size = dplyr::if_else(.data$.is_main_row, .data$effect_size, NA_character_),
                                     Notes = dplyr::if_else(.data$.is_main_row, .data$Notes, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**") %>%
    gtsummary::modify_caption(.caption_text())

  if (isTRUE(AddEffectSize)) {
    tbl <- tbl %>% gtsummary::modify_header(effect_size ~ "**Effect size**")
  } else {
    tbl <- tbl %>% gtsummary::modify_table_body(~ dplyr::select(.x, -dplyr::any_of("effect_size")))
  }

  if (isTRUE(IncludeNotes)) {
    tbl <- tbl %>% gtsummary::modify_header(Notes ~ "**Notes**")
  } else {
    tbl <- tbl %>% gtsummary::modify_table_body(~ dplyr::select(.x, -dplyr::any_of("Notes")))
  }

  # Pairwise contrasts -------------------------------------------------------
  if (isTRUE(AddPairwise) && nlevels(df[[CompVariable]]) > 1) {

    lvls <- levels(droplevels(.as_factor(df[[CompVariable]])))
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% lvls) stop("Referent level not found: ", Referent)
      lapply(setdiff(lvls, Referent), function(x) c(Referent, x))
    } else {
      utils::combn(lvls, 2, simplify = FALSE)
    }

    pw_long <- purrr::map_dfr(Variables, function(var) {

      is_cont <- isTRUE(treat_as_continuous[var])

      out_template <- purrr::map_dfr(combos, function(cp) {
        tibble::tibble(
          variable = var,
          contrast_label = .pair_label(cp[1], cp[2]),
          key = .pair_key(cp[1], cp[2]),
          p_val = NA_real_
        )
      })

      # Continuous -----------------------------------------------------------
      if (isTRUE(is_cont)) {

        # Adjusted continuous: emmeans on lm (robust vcov passed when Parametric=FALSE)
        if (!is.null(Covariates)) {
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))

          if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
            df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
          }

          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
          }

          fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                          error = function(e) NULL)
          if (is.null(fit)) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

          V <- NULL
          if (!isTRUE(Parametric)) V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)

          emm <- tryCatch(
            if (!is.null(V)) emmeans::emmeans(fit, specs = CompVariable, vcov. = V)
            else emmeans::emmeans(fit, specs = CompVariable),
            error = function(e) NULL
          )
          if (is.null(emm)) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

          ctr <- tryCatch(
            if (!is.null(Referent)) emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)
            else emmeans::contrast(emm, method = "pairwise"),
            error = function(e) NULL
          )
          if (is.null(ctr)) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

          res <- tryCatch(as.data.frame(summary(ctr, adjust = ifelse(PairwiseMethod == "none", "none", PairwiseMethod))),
                          error = function(e) NULL)
          if (is.null(res) || !("p.value" %in% names(res))) {
            return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
          }

          got <- tibble::tibble(
            key = vapply(as.character(res$contrast), function(s) {
              parts <- strsplit(.norm(s), "\\s*-\\s*")[[1]]
              if (length(parts) != 2) return(NA_character_)
              a <- .strip_group_prefix(parts[1], CompVariable)
              b <- .strip_group_prefix(parts[2], CompVariable)
              .pair_key(a, b)
            }, character(1)),
            p_val = .clamp_p(res$p.value)
          )

          return(out_template %>%
                   dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                   dplyr::select("variable", "contrast_label", "p_val"))
        }

        # Unadjusted continuous pairwise ------------------------------------
        df_pw <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_pw[[CompVariable]] <- droplevels(.as_factor(df_pw[[CompVariable]]))
        if (nlevels(df_pw[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        if (!is.null(Referent) && Referent %in% levels(df_pw[[CompVariable]])) {
          df_pw[[CompVariable]] <- stats::relevel(df_pw[[CompVariable]], ref = Referent)
        }

        pw <- NULL
        if (isTRUE(Parametric)) {
          pw <- tryCatch(stats::pairwise.t.test(df_pw[[var]], df_pw[[CompVariable]],
                                                p.adjust.method = ifelse(PairwiseMethod == "none", "none", PairwiseMethod),
                                                pool.sd = FALSE),
                         error = function(e) NULL)
        } else {
          pw <- tryCatch(stats::pairwise.wilcox.test(df_pw[[var]], df_pw[[CompVariable]],
                                                     p.adjust.method = ifelse(PairwiseMethod == "none", "none", PairwiseMethod)),
                         error = function(e) NULL)
        }

        if (is.null(pw) || is.null(pw$p.value)) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        # Convert p.value matrix to keys
        mat <- pw$p.value
        got <- tibble::tibble(
          a = rep(rownames(mat), times = ncol(mat)),
          b = rep(colnames(mat), each = nrow(mat)),
          p_val = as.numeric(as.vector(mat))
        ) %>%
          dplyr::mutate(
            key = purrr::map2_chr(.data$a, .data$b, ~ .pair_key(.x, .y)),
            p_val = .clamp_p(.data$p_val)
          )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Categorical ----------------------------------------------------------
      # Unadjusted pairwise: run CatMethod (auto/chisq/fisher) on each pair
      if (is.null(Covariates)) {
        got <- purrr::map_dfr(combos, function(cp) {
          sub <- df[df[[CompVariable]] %in% cp, , drop = FALSE]
          sub <- sub[stats::complete.cases(sub[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
          sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
          sub[[var]] <- droplevels(.as_factor(sub[[var]]))
          if (nlevels(sub[[CompVariable]]) < 2 || nlevels(sub[[var]]) < 2) {
            return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
          }
          tab <- table(sub[[var]], sub[[CompVariable]])
          tst <- .cat_global_test(tab, method = CatMethod)
          tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = .clamp_p(tst$p))
        })

        if (!identical(PairwiseMethod, "none")) {
          got$p_val <- stats::p.adjust(got$p_val, method = PairwiseMethod)
        }

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Adjusted categorical pairwise: LR logic (binary vs multinom per pair)
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))

      if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
        df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
      }

      got <- .cat_pairwise_lr(df_cc, var, CompVariable, Covariates, combos, PairwiseMethod)

      out_template %>%
        dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
        dplyr::select("variable", "contrast_label", "p_val")
    })

    if (nrow(pw_long) > 0) {
      contrast_levels <- unique(pw_long$contrast_label)
      safe_names <- make.unique(paste0("pw_", make.names(contrast_levels)))
      map <- tibble::tibble(contrast_label = contrast_levels, col_safe = safe_names)

      pw_long <- dplyr::left_join(pw_long, map, by = "contrast_label")

      pw_wide <- pw_long %>%
        dplyr::select("variable", "col_safe", "p_val") %>%
        tidyr::pivot_wider(
          id_cols = "variable",
          names_from = "col_safe",
          values_from = "p_val",
          values_fn = list(p_val = function(x) x[1])
        )

      tbl <- tbl %>%
        gtsummary::modify_table_body(~ .x %>%
                                       dplyr::left_join(pw_wide, by = "variable") %>%
                                       dplyr::mutate(
                                         dplyr::across(
                                           tidyselect::all_of(setdiff(names(pw_wide), "variable")),
                                           ~ dplyr::if_else(.data$.is_main_row, .clamp_p(as.numeric(.)), NA_real_)
                                         )
                                       )
        )

      for (col in setdiff(names(pw_wide), "variable")) {
        lab <- map$contrast_label[match(col, map$col_safe)]
        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(!!rlang::sym(col) ~ function(x) .fmt_p(x, digits = pDigits)) %>%
          gtsummary::modify_header(!!rlang::sym(col) := paste0("**", lab, "**"))
      }
    }
  }

  # Ensure Notes is optional AND always last ---------------------------------
  if (isTRUE(IncludeNotes) && ("Notes" %in% names(tbl$table_body))) {
    tbl <- tbl %>%
      gtsummary::modify_table_body(~ dplyr::relocate(.x, .data$Notes, .after = dplyr::last_col()))
  }

  tbl
}
