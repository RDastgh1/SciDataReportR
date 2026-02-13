#' Make comparison table with optional covariate adjustment, effect sizes, and pairwise tests
#'
#' Create a label friendly comparison table (via **gtsummary**) that summarizes one or more variables
#' across levels of a grouping variable, optionally adjusts for covariates, adds global p-values,
#' effect sizes, and pairwise p-values (with optional multiplicity adjustment).
#'
#' @description
#' This function is designed for “Table 1 style” descriptive comparisons, with flexible choices
#' for statistical tests and a predictable output structure that stays stable across variable types.
#'
#' @details
#' ## How tests are chosen (the guide you wish every table function came with)
#'
#' ### Continuous variables (numeric with more than 2 distinct values)
#' **No covariates**
#' - `Parametric = TRUE`
#'   - 2 groups: Welch two sample t-test
#'   - 3+ groups: one-way ANOVA
#' - `Parametric = FALSE`
#'   - 2 groups: Wilcoxon rank-sum (Mann–Whitney)
#'   - 3+ groups: Kruskal–Wallis
#'
#' **With covariates**
#' - Continuous outcomes are modeled with `lm(outcome ~ group + covariates)`.
#' - Global group test depends on `ContinuousAdjustedMethod`:
#'   - `"ancova"`: Type II ANCOVA via `car::Anova(fit, type = 2)`.
#'   - `"robust_ancova"`: robust (heteroskedasticity consistent) Wald F test using HC3
#'     (`sandwich::vcovHC(type = "HC3")` + `car::linearHypothesis()`).
#'
#' **Pairwise p-values**
#' - No covariates, parametric: `pairwise.t.test(pool.sd = FALSE)` (Welch style).
#' - No covariates, nonparametric: per-pair `wilcox.test()`.
#' - With covariates: estimated marginal means and contrasts via **emmeans** on the fitted `lm`.
#'   If `ContinuousAdjustedMethod="robust_ancova"`, the same HC3 covariance is passed to emmeans
#'   so pairwise tests align with robust inference.
#'
#' **Important note about “nonparametric + covariates”**
#' - There is no single universally accepted “rank ANCOVA” that behaves like Kruskal–Wallis
#'   while adjusting for covariates and still meaningfully targets adjusted means.
#' - This function’s approach (robust ANCOVA) keeps the adjusted mean interpretation and uses
#'   robust standard errors for inference. It is not rank based.
#'
#' ---
#'
#' ### Categorical variables (including dichotomous)
#' **Key principle:** The `Parametric` argument does *not* change categorical tests.
#' Categorical tests are controlled by `CategoricalMethod` and covariate presence.
#'
#' **No covariates**
#' - `CategoricalMethod = "chisq"`: Pearson chi-squared without Yates correction
#'   (`chisq.test(correct = FALSE)`) to match common Stata defaults.
#' - `CategoricalMethod = "fisher"`:
#'   - 2x2: Fisher’s exact
#'   - larger RxC: Fisher with Monte Carlo simulation (`simulate.p.value=TRUE`) using `FisherB`.
#'
#' **With covariates**
#' - If the outcome has exactly 2 levels, the global adjusted test is:
#'   - Logistic regression likelihood ratio test for the grouping factor via `drop1(..., test="Chisq")`.
#' - If the outcome has 3+ levels, covariate-adjusted global tests are **not implemented**
#'   (multinomial models are possible but require additional design decisions). The function will
#'   fall back to the unadjusted categorical test for that variable.
#'
#' **Pairwise p-values**
#' - No covariates: pairwise uses the *same* method selected for that variable’s global test
#'   (chi-squared or Fisher), so the pairwise columns match the Test column.
#' - With covariates + binary outcome: pairwise is computed from the fitted logistic regression
#'   using Wald z-tests of linear predictor contrasts on the logit scale. This is the closest match
#'   to typical Stata pairwise coefficient contrast workflows.
#'
#' ---
#'
#' ### Effect sizes
#' - Continuous, 2 groups, unadjusted parametric: \|Cohen’s d\|.
#' - Continuous, 3+ groups, unadjusted parametric: eta squared (η²).
#' - Continuous, unadjusted nonparametric: epsilon squared (ε²) derived from Kruskal–Wallis.
#' - Continuous with covariates: partial eta squared (partial η²) for the grouping term (Type II).
#' - Categorical: Cramer’s V (unadjusted).
#'
#' **Note:** Adjusted odds ratios or standardized OR effect sizes for logistic regression are not
#' returned by this function (by design). If you want adjusted ORs, compute them from the fitted
#' model outside this table.
#'
#' ## References (methods and software)
#' - Pearson chi-squared test: Pearson K. (1900) \doi{10.1098/rspl.1900.0024}
#' - Fisher’s exact test: Fisher RA. (1922) \doi{10.2307/2340521}
#' - Welch t-test: Welch BL. (1947) \doi{10.2307/2332510}
#' - Kruskal–Wallis: Kruskal WH, Wallis WA. (1952) \doi{10.1080/01621459.1952.10483441}
#' - Robust (HC) covariance estimators: MacKinnon JG, White H. (1985) \doi{10.1016/0304-4076(85)90158-7}
#' - Type II/III ANOVA in practice: Fox J, Weisberg S. *An R Companion to Applied Regression* (car package ecosystem).
#' - Estimated marginal means: Lenth RV. (2021) emmeans package vignette and documentation.
#' - Effect sizes in R: Ben-Shachar MS, Lüdecke D, Makowski D. (2020) \doi{10.21105/joss.02815} (effectsize package).
#'
#' @param DataFrame Data frame containing variables to summarize.
#' @param CompVariable Grouping/comparison variable name (character scalar). If `NULL` and
#'   `IncludeOverallStats=FALSE`, overall-only mode is used.
#' @param Variables Character vector of variable names to include in the table.
#' @param Covariates Optional character vector of covariate names for adjusted analyses.
#' @param ValueDigits Digits for continuous summary statistics in `tbl_summary()`.
#' @param pDigits Digits shown for formatted p-values (uses `format.pval()`).
#' @param AddEffectSize Logical; add effect size columns.
#' @param EffectSizeDigits Digits shown for effect sizes.
#' @param AddPairwise Logical; add pairwise p-value columns.
#' @param PairwiseMethod Multiplicity correction passed to `p.adjust()` and emmeans summaries.
#'   Use `"none"` for no adjustment.
#' @param Parametric Logical controlling tests for *continuous* outcomes when there are no covariates.
#'   Categorical tests are controlled by `CategoricalMethod` and covariate presence.
#' @param ParametricDisplay Logical controlling whether continuous variables are displayed as mean (SD)
#'   or median [IQR]. Defaults to `Parametric`.
#' @param IncludeOverallN Logical; add an overall N column via `gtsummary::add_n()`.
#' @param IncludeMissing Logical; include missing row (Unknown) in gtsummary output.
#' @param suppress_warnings Logical; suppress gtsummary warnings.
#' @param Referent Optional reference level (character scalar) for contrasts.
#' @param IncludeOverallStats Logical; force overall-only summary even if `CompVariable` is provided.
#' @param ShowPositiveBinaryOnLabel Logical; for binary categorical variables, show only the “positive”
#'   level on the label row (TRUE/1/YES etc.).
#' @param BinaryPairwiseScale Currently `"logit"` only. Pairwise for adjusted binary outcomes is on logit scale.
#' @param ContinuousAdjustedMethod `"robust_ancova"` (default) or `"ancova"`. Controls inference for
#'   continuous outcomes with covariates.
#' @param CategoricalMethod `"chisq"` (default) or `"fisher"` for unadjusted categorical tests.
#' @param FisherB Number of Monte Carlo replicates for Fisher simulation when `CategoricalMethod="fisher"`
#'   and the table is larger than 2x2.
#' @param FisherSeed Optional seed for Fisher simulation reproducibility.
#'
#' @return A `gtsummary::tbl_summary` object with additional columns for p-values, test labels,
#'   optional effect sizes, and optional pairwise p-values.
#'
#' @examples
#' suppressPackageStartupMessages({
#'   library(dplyr)
#' })
#'
#' set.seed(1)
#' n <- 400
#' df_fake <- tibble::tibble(
#'   CTQ_sexual_abuse_3CAT = factor(sample(c("None or minimal", "Low to moderate", "Moderate to extreme"),
#'                                         n, replace = TRUE, prob = c(0.65, 0.18, 0.17))),
#'   age = round(rnorm(n, 45, 12)),
#'   sex = factor(sample(c("F", "M"), n, replace = TRUE)),
#'   CTQ_physical_abuse_di = rbinom(n, 1, prob = dplyr::case_when(
#'     CTQ_sexual_abuse_3CAT == "None or minimal" ~ 0.35,
#'     CTQ_sexual_abuse_3CAT == "Low to moderate" ~ 0.55,
#'     TRUE                                      ~ 0.60
#'   )),
#'   CTQ_physical_abuse = rnorm(n, mean = dplyr::case_when(
#'     CTQ_sexual_abuse_3CAT == "None or minimal" ~ 5.8,
#'     CTQ_sexual_abuse_3CAT == "Low to moderate" ~ 6.2,
#'     TRUE                                      ~ 7.3
#'   ), sd = 2.5)
#' )
#'
#' # Unadjusted categorical (Pearson chi-squared) with pairwise
#' t_cat_unadj <- MakeComparisonTable(
#'   df_fake,
#'   Variables = "CTQ_physical_abuse_di",
#'   CompVariable = "CTQ_sexual_abuse_3CAT",
#'   AddPairwise = TRUE,
#'   PairwiseMethod = "none",
#'   Parametric = FALSE,
#'   CategoricalMethod = "chisq"
#' )
#' t_cat_unadj$table_body %>% dplyr::select(variable, p.value, starts_with("pw_"))
#'
#' # Adjusted binary outcome (logistic LR global + Wald z pairwise)
#' t_cat_adj <- MakeComparisonTable(
#'   df_fake,
#'   Variables = "CTQ_physical_abuse_di",
#'   CompVariable = "CTQ_sexual_abuse_3CAT",
#'   Covariates = c("age", "sex"),
#'   AddPairwise = TRUE,
#'   PairwiseMethod = "none",
#'   Parametric = TRUE
#' )
#' t_cat_adj$table_body %>% dplyr::select(variable, p.value, starts_with("pw_"))
#'
#' # Continuous adjusted robust ANCOVA + emmeans pairwise
#' t_cont_adj <- MakeComparisonTable(
#'   df_fake,
#'   Variables = "CTQ_physical_abuse",
#'   CompVariable = "CTQ_sexual_abuse_3CAT",
#'   Covariates = c("age", "sex"),
#'   AddPairwise = TRUE,
#'   PairwiseMethod = "bonferroni",
#'   Parametric = FALSE, # affects display + unadjusted tests only
#'   ContinuousAdjustedMethod = "robust_ancova"
#' )
#' t_cont_adj$table_body %>% dplyr::select(variable, p.value, starts_with("pw_"))
#'
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
    ShowPositiveBinaryOnLabel = TRUE,
    BinaryPairwiseScale  = c("logit"),
    ContinuousAdjustedMethod = c("robust_ancova", "ancova"),
    CategoricalMethod    = c("chisq", "fisher"),
    FisherB              = 1e5,
    FisherSeed           = NULL
) {

  # Validate inputs ---------------------------------------------------------

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric
  if (is.null(PairwiseMethod)) PairwiseMethod <- "none"

  Variables <- as.character(Variables)
  if (length(Variables) < 1) stop("Variables must be a character vector with length >= 1.")

  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale)
  ContinuousAdjustedMethod <- match.arg(ContinuousAdjustedMethod)
  CategoricalMethod <- match.arg(CategoricalMethod)

  valid_methods <- c("none", stats::p.adjust.methods)
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1 || !(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "),
         ". Got: ", paste(PairwiseMethod, collapse = ", "))
  }

  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "effectsize", "sandwich"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  comp_present <- !is.null(CompVariable) && is.character(CompVariable) &&
    length(CompVariable) == 1 && CompVariable %in% names(DataFrame)
  overall_mode <- isTRUE(IncludeOverallStats) || !comp_present

  if (!overall_mode && !CompVariable %in% names(DataFrame))
    stop("Grouping variable not found: ", CompVariable)

  if (!all(Variables %in% names(DataFrame)))
    stop("Variable(s) not found: ",
         paste(setdiff(Variables, names(DataFrame)), collapse = ", "))

  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame)))
    stop("Covariate(s) not found: ",
         paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))

  # Prepare data ------------------------------------------------------------

  .as_factor <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  .norm <- function(x) gsub("\\s+", " ", trimws(as.character(x)))
  .pair_label <- function(a, b) paste(a, b, sep = " - ")
  .pair_key <- function(a, b) paste(sort(c(.norm(a), .norm(b))), collapse = "||")

  .clean_tab <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    tab
  }

  .fmt_p <- function(p, digits = 3) {
    ifelse(is.na(p), NA_character_, format.pval(p, digits = digits, eps = 10^(-digits)))
  }

  .cramers_v <- function(tab) {
    tab <- .clean_tab(tab)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)

    chi <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE)$statistic, error = function(e) NA_real_))
    if (is.na(chi)) return(NA_real_)
    n <- sum(tab)
    m <- min(nrow(tab), ncol(tab)) - 1
    if (m <= 0 || n <= 0) return(NA_real_)
    sqrt(as.numeric(chi) / (n * m))
  }

  .cat_test <- function(tab, method = c("chisq", "fisher"), B = 1e5, seed = NULL) {
    method <- match.arg(method)
    tab <- .clean_tab(tab)

    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(list(p = NA_real_, label = "Insufficient data", method = method))
    }

    if (method == "chisq") {
      chi <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL))
      p <- if (!is.null(chi) && !is.null(chi$p.value)) as.numeric(chi$p.value) else NA_real_
      return(list(p = p, label = "Chi-squared", method = "chisq"))
    }

    if (!is.null(seed)) set.seed(seed)
    fish <- tryCatch(
      if (nrow(tab) == 2 && ncol(tab) == 2) stats::fisher.test(tab)
      else stats::fisher.test(tab, simulate.p.value = TRUE, B = B),
      error = function(e) NULL
    )
    p <- if (!is.null(fish) && !is.null(fish$p.value)) as.numeric(fish$p.value) else NA_real_
    lab <- if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"
    list(p = p, label = lab, method = "fisher")
  }

  .robust_group_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)

    cn <- names(stats::coef(fit_lm))
    idx <- grep(paste0("^", term), cn)
    if (!length(idx)) return(NA_real_)

    L <- matrix(0, nrow = length(idx), ncol = length(cn))
    L[cbind(seq_along(idx), idx)] <- 1

    lh <- tryCatch(
      car::linearHypothesis(fit_lm, hypothesis.matrix = L, vcov. = V, test = "F"),
      error = function(e) NULL
    )
    if (is.null(lh) || !"Pr(>F)" %in% colnames(lh)) return(NA_real_)
    as.numeric(lh[2, "Pr(>F)"])
  }

  .pairwise_logit_xb <- function(gm, df_cc, compvar, PairwiseMethod, Referent = NULL) {

    comp_levels <- levels(df_cc[[compvar]])
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% comp_levels) stop("Referent level not found: ", Referent)
      lapply(setdiff(comp_levels, Referent), function(x) c(Referent, x))
    } else {
      utils::combn(comp_levels, 2, simplify = FALSE)
    }

    tt <- stats::delete.response(stats::terms(gm))
    beta <- stats::coef(gm)
    V <- tryCatch(stats::vcov(gm), error = function(e) NULL)
    if (is.null(V)) return(NULL)

    nd <- df_cc[rep(1, length(comp_levels)), , drop = FALSE]
    nd[[compvar]] <- factor(comp_levels, levels = comp_levels)

    X <- tryCatch(stats::model.matrix(tt, data = nd, xlev = gm$xlevels), error = function(e) NULL)
    if (is.null(X)) return(NULL)

    idx <- setNames(seq_along(comp_levels), comp_levels)

    out <- purrr::map_dfr(combos, function(cp) {
      i <- idx[[cp[1]]]
      j <- idx[[cp[2]]]
      d <- X[i, , drop = FALSE] - X[j, , drop = FALSE]
      est <- as.numeric(d %*% beta)
      se <- sqrt(as.numeric(d %*% V %*% t(d)))
      z <- est / se
      p <- 2 * stats::pnorm(-abs(z))

      tibble::tibble(
        contrast_label = .pair_label(cp[1], cp[2]),
        key = .pair_key(cp[1], cp[2]),
        p_val = p
      )
    })

    if (!identical(PairwiseMethod, "none")) {
      out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    }
    out
  }

  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)

  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- .as_factor(df[[CompVariable]])

  # Drop constants ----------------------------------------------------------

  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]

  drop <- setdiff(Variables, keep)
  if (length(drop)) warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))

  Variables <- keep
  if (!length(Variables)) stop("No variables left to summarise after dropping constants.")

  # Variable typing ---------------------------------------------------------

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

  # Positive-level display --------------------------------------------------

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
        if (any(vals == 1)) vmap[[v]] <- 1 else vmap[[v]] <- max(vals)
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

  # Build outputs -----------------------------------------------------------

  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

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
    cap <- sprintf("Overall Summary (display: %s)", if (ParametricDisplay) "mean (SD)" else "median [IQR]")
    return(tbl %>% gtsummary::modify_caption(cap))
  }

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

  # Global p-values ---------------------------------------------------------

  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(
        variable = var, p_unadj = NA_real_, p_adj = NA_real_,
        test_label = "Insufficient groups", cat_method = NA_character_
      ))
    }

    if (is_cont) {
      k <- nlevels(df_vg[[CompVariable]])

      if (!is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Insufficient data (adjusted)", cat_method = NA_character_
          ))
        }

        fit <- tryCatch(
          stats::lm(reformulate(c(CompVariable, Covariates), response = var), data = df_cc),
          error = function(e) NULL
        )
        if (is.null(fit)) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Model failed (adjusted)", cat_method = NA_character_
          ))
        }

        p_un <- tryCatch(
          summary(stats::aov(reformulate(CompVariable, response = var), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )

        if (ContinuousAdjustedMethod == "ancova") {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) as.numeric(a2[CompVariable, "Pr(>F)"]) else NA_real_
          return(tibble::tibble(
            variable = var, p_unadj = p_un, p_adj = p_adj,
            test_label = "ANCOVA (Type II)", cat_method = NA_character_
          ))
        }

        p_rb <- .robust_group_p(fit, CompVariable)
        return(tibble::tibble(
          variable = var, p_unadj = p_un, p_adj = p_rb,
          test_label = "Robust ANCOVA (HC3 Wald)", cat_method = NA_character_
        ))
      }

      if (Parametric) {
        if (k == 2) {
          p_un <- tryCatch(
            stats::t.test(reformulate(CompVariable, response = var), data = df_vg, var.equal = FALSE)$p.value,
            error = function(e) NA_real_
          )
          return(tibble::tibble(
            variable = var, p_unadj = p_un, p_adj = NA_real_,
            test_label = "Welch t-test", cat_method = NA_character_
          ))
        }
        p_un <- tryCatch(
          summary(stats::aov(reformulate(CompVariable, response = var), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )
        return(tibble::tibble(
          variable = var, p_unadj = p_un, p_adj = NA_real_,
          test_label = "ANOVA", cat_method = NA_character_
        ))
      }

      p_un <- tryCatch(
        if (k == 2) stats::wilcox.test(reformulate(CompVariable, response = var), data = df_vg)$p.value
        else stats::kruskal.test(reformulate(CompVariable, response = var), data = df_vg)$p.value,
        error = function(e) NA_real_
      )
      return(tibble::tibble(
        variable = var, p_unadj = p_un, p_adj = NA_real_,
        test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis",
        cat_method = NA_character_
      ))
    }

    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- .cat_test(tab, method = CategoricalMethod, B = FisherB, seed = FisherSeed)
    p_un <- tst$p
    test_label <- tst$label
    p_adj <- NA_real_
    cat_method <- tst$method

    if (!is.null(Covariates) && nlevels(droplevels(x)) == 2) {

      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- .as_factor(df_cc[[var]])
      if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

      if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5 && nlevels(df_cc[[var]]) == 2) {
        gm <- tryCatch(
          stats::glm(reformulate(c(CompVariable, Covariates), response = var),
                     data = df_cc, family = stats::binomial()),
          error = function(e) NULL
        )
        if (!is.null(gm)) {
          d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
          if (!is.null(d1) && CompVariable %in% rownames(d1)) {
            p_adj <- as.numeric(d1[CompVariable, "Pr(>Chi)"])
            test_label <- "Logistic regression (LR)"
          }
        }
      }
    }

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = test_label, cat_method = cat_method)
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::group_by(variable) %>%
                                   dplyr::mutate(.is_main_row = dplyr::row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(
                                     p.value = dplyr::coalesce(p_adj, p_unadj),
                                     p.value_fmt = dplyr::if_else(.is_main_row, .fmt_p(p.value, digits = pDigits), NA_character_),
                                     Test = dplyr::if_else(.is_main_row, test_label, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows    = .data$.is_main_row & !is.na(.data$p.value) & .data$p.value <= 0.05,
      text_format = "bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  # Effect sizes ------------------------------------------------------------

  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {

      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "Insufficient groups"))
      }

      if (isTRUE(treat_as_continuous[var])) {

        if (!is.null(Covariates)) {
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
          if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "partial η²"))
          }

          fit <- tryCatch(stats::lm(reformulate(c(CompVariable, Covariates), response = var), data = df_cc),
                          error = function(e) NULL)
          a2  <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          et  <- tryCatch(effectsize::eta_squared(a2, partial = TRUE), error = function(e) NULL)

          val <- NA_real_
          if (!is.null(et)) {
            idx <- if ("Parameter" %in% names(et)) et$Parameter == CompVariable else rownames(et) == CompVariable
            val <- suppressWarnings(et$Eta2_partial[idx][1])
          }
          return(tibble::tibble(variable = var, effect_size = val, es_method = "partial η²"))
        }

        k <- nlevels(df_vg[[CompVariable]])
        if (Parametric) {
          if (k == 2) {
            val <- tryCatch(
              abs(effectsize::cohens_d(reformulate(CompVariable, response = var), data = df_vg)$Cohens_d),
              error = function(e) NA_real_
            )
            return(tibble::tibble(variable = var, effect_size = val, es_method = "|d|"))
          }
          val <- tryCatch(
            effectsize::eta_squared(stats::aov(reformulate(CompVariable, response = var), data = df_vg),
                                    partial = FALSE)$Eta2[1],
            error = function(e) NA_real_
          )
          return(tibble::tibble(variable = var, effect_size = val, es_method = "η²"))
        }

        n <- nrow(df_vg)
        H <- tryCatch(stats::kruskal.test(reformulate(CompVariable, response = var), data = df_vg)$statistic,
                      error = function(e) NA_real_)
        eps2 <- suppressWarnings(as.numeric((H - k + 1) / (n - k)))
        return(tibble::tibble(variable = var, effect_size = eps2, es_method = "ε²"))
      }

      x <- .as_factor(df_vg[[var]])
      g <- .as_factor(df_vg[[CompVariable]])
      tab <- table(x, g)
      v <- .cramers_v(tab)
      tibble::tibble(variable = var, effect_size = v, es_method = "Cramer's V")
    })

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     dplyr::left_join(es_df, by = "variable") %>%
                                     dplyr::mutate(
                                       effect_size = dplyr::if_else(.data$.is_main_row, effect_size, NA_real_),
                                       ES_Method   = dplyr::if_else(.data$.is_main_row, es_method, NA_character_)
                                     )
      ) %>%
      gtsummary::modify_fmt_fun(effect_size ~ function(x)
        ifelse(is.na(x), NA_character_, formatC(x, digits = EffectSizeDigits, format = "f"))) %>%
      gtsummary::modify_header(effect_size ~ "**Effect size**") %>%
      gtsummary::modify_header(ES_Method   ~ "**ES method**")
  }

  # Pairwise contrasts ------------------------------------------------------

  if (AddPairwise && nlevels(df[[CompVariable]]) > 1) {

    lvls <- levels(droplevels(.as_factor(df[[CompVariable]])))
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% lvls) stop("Referent level not found: ", Referent)
      lapply(setdiff(lvls, Referent), function(x) c(Referent, x))
    } else {
      utils::combn(lvls, 2, simplify = FALSE)
    }

    cat_method_map <- pdat %>% dplyr::select(variable, cat_method, test_label)

    pw_long <- purrr::map_dfr(Variables, function(var) {

      if (length(combos) == 0) {
        return(tibble::tibble(variable = character(), contrast_label = character(), p_val = numeric()))
      }

      is_cont <- isTRUE(treat_as_continuous[var])

      out_template <- purrr::map_dfr(combos, function(cp) {
        tibble::tibble(
          variable = var,
          contrast_label = .pair_label(cp[1], cp[2]),
          key = .pair_key(cp[1], cp[2]),
          p_val = NA_real_
        )
      })

      if (is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        fit <- tryCatch(stats::lm(reformulate(c(CompVariable, Covariates), response = var), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        vc <- NULL
        if (ContinuousAdjustedMethod == "robust_ancova") {
          V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
          if (!is.null(V)) vc <- V
        }

        emm <- tryCatch(
          if (is.null(vc)) emmeans::emmeans(fit, specs = CompVariable)
          else emmeans::emmeans(fit, specs = CompVariable, vcov. = vc),
          error = function(e) NULL
        )
        if (is.null(emm)) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        ctr <- tryCatch(
          if (!is.null(Referent)) emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)
          else emmeans::contrast(emm, method = "pairwise"),
          error = function(e) NULL
        )
        if (is.null(ctr)) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        res <- tryCatch(as.data.frame(summary(ctr, adjust = PairwiseMethod)), error = function(e) NULL)
        if (is.null(res) || !("p.value" %in% names(res)) || !("contrast" %in% names(res))) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        got <- tibble::tibble(
          key = vapply(res$contrast, function(s) {
            parts <- strsplit(.norm(s), "\\s*-\\s*")[[1]]
            if (length(parts) != 2) return(NA_character_)
            .pair_key(parts[1], parts[2])
          }, character(1)),
          p_val = res$p.value
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(key, got$key)]) %>%
                 dplyr::select(variable, contrast_label, p_val))
      }

      if (!is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        df_cc[[var]] <- .as_factor(df_cc[[var]])
        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5 && nlevels(df_cc[[var]]) == 2) {
          gm <- tryCatch(
            stats::glm(reformulate(c(CompVariable, Covariates), response = var),
                       data = df_cc, family = stats::binomial()),
            error = function(e) NULL
          )
          if (!is.null(gm)) {
            got <- .pairwise_logit_xb(gm, df_cc, CompVariable, PairwiseMethod, Referent)
            if (!is.null(got) && nrow(got)) {
              return(out_template %>%
                       dplyr::mutate(p_val = got$p_val[match(key, got$key)]) %>%
                       dplyr::select(variable, contrast_label, p_val))
            }
          }
        }

        return(dplyr::select(out_template, variable, contrast_label, p_val))
      }

      if (is_cont && is.null(Covariates)) {

        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, variable, contrast_label, p_val))
        }

        if (Parametric) {
          res <- tryCatch(
            stats::pairwise.t.test(df_vg[[var]], df_vg[[CompVariable]],
                                   p.adjust.method = PairwiseMethod, pool.sd = FALSE),
            error = function(e) NULL
          )
          if (is.null(res) || is.null(res$p.value)) {
            return(dplyr::select(out_template, variable, contrast_label, p_val))
          }

          r <- as.data.frame(as.table(res$p.value))
          r <- r[!is.na(r$Freq), , drop = FALSE]

          out <- purrr::map_dfr(combos, function(cp) {
            pv <- NA_real_
            hit1 <- r$Freq[r$Var1 == cp[2] & r$Var2 == cp[1]]
            hit2 <- r$Freq[r$Var1 == cp[1] & r$Var2 == cp[2]]
            if (length(hit1)) pv <- hit1[1]
            if (is.na(pv) && length(hit2)) pv <- hit2[1]
            tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = pv)
          })
          return(out)
        }

        out <- purrr::map_dfr(combos, function(cp) {
          sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
          p <- tryCatch(stats::wilcox.test(reformulate(CompVariable, response = var), data = sub)$p.value,
                        error = function(e) NA_real_)
          tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = p)
        })
        if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
        return(out)
      }

      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(dplyr::select(out_template, variable, contrast_label, p_val))
      }

      this_method <- cat_method_map$cat_method[match(var, cat_method_map$variable)]
      if (is.na(this_method) || !this_method %in% c("chisq", "fisher")) this_method <- CategoricalMethod

      out <- purrr::map_dfr(combos, function(cp) {
        sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
        sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
        tab2 <- table(.as_factor(sub[[var]]), .as_factor(sub[[CompVariable]]))
        tst2 <- .cat_test(tab2, method = this_method, B = FisherB, seed = FisherSeed)
        tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = tst2$p)
      })

      if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
      out
    })

    if (nrow(pw_long) > 0) {
      contrast_levels <- unique(pw_long$contrast_label)
      safe_names <- make.unique(paste0("pw_", make.names(contrast_levels)))
      map <- tibble::tibble(contrast_label = contrast_levels, col_safe = safe_names)

      pw_long <- dplyr::left_join(pw_long, map, by = "contrast_label")

      pw_wide <- pw_long %>%
        dplyr::select(variable, col_safe, p_val) %>%
        tidyr::pivot_wider(id_cols = variable, names_from = col_safe, values_from = p_val)

      tbl <- tbl %>%
        gtsummary::modify_table_body(~ .x %>%
                                       dplyr::left_join(pw_wide, by = "variable") %>%
                                       dplyr::mutate(
                                         dplyr::across(
                                           tidyselect::all_of(setdiff(names(pw_wide), "variable")),
                                           ~ dplyr::if_else(.data$.is_main_row, ., NA_real_)
                                         )
                                       )
        )

      for (col in setdiff(names(pw_wide), "variable")) {
        lab <- map$contrast_label[match(col, map$col_safe)]

        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(!!rlang::sym(col) ~ function(x) .fmt_p(x, digits = pDigits)) %>%
          gtsummary::modify_table_styling(
            columns = col,
            rows    = .data$.is_main_row & !is.na(.data[[col]]) & .data[[col]] <= 0.05,
            text_format = "bold"
          ) %>%
          gtsummary::modify_header(!!rlang::sym(col) := paste0("**", lab, "**")) %>%
          gtsummary::modify_footnote(
            !!rlang::sym(col) ~ paste0(
              "Pairwise p-value (", PairwiseMethod, "). ",
              "Categorical pairwise uses the same method as **Test**. ",
              if (!is.null(Covariates)) "Adjusted pairwise uses lm-emmeans (continuous) or logistic Wald z (binary)." else ""
            )
          )
      }
    }
  }

  # Return result -----------------------------------------------------------

  cap <- sprintf(
    "Comparison Table (display: %s; Parametric=%s; Covariates=%s; Categorical=%s; Continuous adjusted=%s). Pairwise: %s (%s).",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    Parametric,
    if (is.null(Covariates)) "none" else paste(Covariates, collapse = ", "),
    CategoricalMethod,
    ContinuousAdjustedMethod,
    if (AddPairwise) "included" else "not included",
    PairwiseMethod
  )

  tbl %>% gtsummary::modify_caption(cap)
}
