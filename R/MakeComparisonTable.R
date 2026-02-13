#' Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts
#'
#' Create a label-friendly comparison table (via **gtsummary**) that summarizes variables
#' across groups, adds global tests, optional covariate adjustment, effect sizes, and
#' optional pairwise comparisons. Designed for “publication table” workflows where you
#' want the **Test** column to reflect the method that generated the p-values.
#'
#' @details
#' ## What gets tested (global p-value)
#' The global p-value shown for each variable depends on variable type and options:
#'
#' **Continuous outcomes**
#' - `Parametric = TRUE`, no covariates:
#'   - 2 groups: Welch two-sample t test (`stats::t.test`, unequal variances by default).
#'   - 3+ groups: one-way ANOVA (`stats::aov`).
#' - Covariates present:
#'   - `Parametric = TRUE`: ANCOVA via linear model (`stats::lm`) with Type II test for the grouping term
#'     (`car::Anova(type = 2)`).
#'   - `Parametric = FALSE`: Robust ANCOVA via `stats::lm` with **HC3** robust covariance
#'     (`sandwich::vcovHC(type = "HC3")`) and a Wald F test using `car::linearHypothesis`.
#'
#' **Categorical outcomes** (dichotomous or multi-category)
#' - Global test is Pearson chi-squared (`stats::chisq.test(correct = FALSE)`) by default.
#' - If `CatMethod = "auto"`, Fisher’s exact test is used when expected counts are small.
#' - If covariates are present and the outcome is **binary**, an adjusted global test is
#'   provided via logistic regression likelihood-ratio test (`stats::drop1(test="Chisq")`).
#'
#' ## What gets tested (pairwise p-values)
#' Pairwise p-values are intended to match the method implied by the **Test** column:
#'
#' **Continuous outcomes**
#' - No covariates, `Parametric = TRUE`: `stats::pairwise.t.test(pool.sd = FALSE)`.
#' - No covariates, `Parametric = FALSE`: `stats::pairwise.wilcox.test`.
#' - With covariates: pairwise comparisons of adjusted means using `emmeans::emmeans` +
#'   `emmeans::pairs` (or treatment-vs-control if `Referent` is set).
#'
#' **Important note about adjusted continuous pairwise**
#' This function uses complete-case rows for each outcome (outcome + group + covariates).
#' After complete-case filtering, some covariates can become constant (for example, all remaining
#' participants share the same Sex), or some group levels can disappear. That can make the model
#' rank-deficient and prevent `emmeans` from estimating certain contrasts.
#'
#' To avoid blank pairwise columns due to purely technical issues, this function:
#' - Removes the outcome variable from covariates per-row (prevents `y` on both sides),
#' - Drops covariates that are constant after complete-case filtering,
#' - Falls back to unadjusted pairwise if no usable covariates remain for that outcome.
#'
#' If pairwise values are still `NA`, that typically reflects truly non-estimable contrasts
#' (data limitation), not a join or formatting failure.
#'
#' **Categorical outcomes**
#' - No covariates: pairwise chi-squared or Fisher (per `CatMethod`) on 2x2 or Rx2 subtables.
#' - With covariates and binary outcome: pairwise Wald z tests of linear predictor contrasts
#'   from logistic regression (coefficient-contrast behavior).
#'
#' ## Effect sizes
#' - Continuous, 2 groups, unadjusted: |Cohen’s d| (`effectsize::cohens_d`).
#' - Continuous, 3+ groups, unadjusted parametric: eta-squared (`effectsize::eta_squared`).
#' - Continuous, adjusted: partial eta-squared from Type II ANOVA table.
#' - Continuous, nonparametric: epsilon-squared approximation from Kruskal-Wallis.
#' - Categorical: Cramer’s V.
#'
#' ## Why the backticks matter
#' This function supports non-syntactic column names (e.g., `"sCD163(ng/mL)"`, `"PSQI TOTAL"`).
#' Any internally constructed formulas always backtick variable names so model calls do not fail.
#'
#' ## Method references
#' - Pearson chi-squared and Fisher’s exact tests: Agresti A. *Categorical Data Analysis*.
#' - ANCOVA / linear models: see `stats::lm`, `stats::aov`.
#' - Type II tests: Fox J, Weisberg S. *car* package documentation for `Anova`.
#' - Estimated marginal means and contrasts: Lenth RV. *emmeans* package documentation.
#' - Robust (HC3) covariance: Long & Ervin (2000) and the `sandwich` package documentation.
#'
#' @param DataFrame A data frame.
#' @param CompVariable Grouping variable name (character scalar). If `NULL` or not present and
#'   `IncludeOverallStats = TRUE`, an overall-only table is returned.
#' @param Variables Character vector of variables to include. You may also pass additional
#'   variable names as unnamed arguments via `...` (convenience).
#' @param ... Optional additional variable names (character scalars).
#' @param Covariates Optional covariate names (character vector) used for adjusted models.
#'   Any variables listed here are automatically removed from `Variables` so they do not appear
#'   as rows in the output table.
#' @param ValueDigits Digits for summary statistics (default 2).
#' @param pDigits Digits to display p-values (default 3). Formatting uses `format.pval`.
#' @param AddEffectSize Add effect sizes column (default FALSE).
#' @param EffectSizeDigits Digits for effect sizes (default 2).
#' @param AddPairwise Add pairwise p-value columns (default FALSE).
#' @param PairwiseMethod P-adjustment method (default "bonferroni"). Use "none" for no adjustment.
#'   Must be one of `"none"` or `stats::p.adjust.methods`.
#' @param Parametric If TRUE, use parametric tests for continuous outcomes where relevant.
#'   If FALSE and covariates are present, robust ANCOVA is used for continuous outcomes.
#' @param ParametricDisplay If TRUE show mean (SD); if FALSE show median [IQR]. Defaults to `Parametric`.
#' @param IncludeOverallN If TRUE add overall N column.
#' @param IncludeMissing If TRUE include missing/unknown category rows in summaries.
#' @param suppress_warnings If TRUE suppress gtsummary warnings.
#' @param Referent Optional reference level for pairwise contrasts (character scalar). If set,
#'   pairwise comparisons are against this referent.
#' @param IncludeOverallStats If TRUE return overall-only summary (ignores grouping).
#' @param ShowPositiveBinaryOnLabel If TRUE, display only the “positive” level for binary categorical variables.
#' @param BinaryPairwiseScale For adjusted binary outcomes: "logit" (default) contrasts on the logit scale,
#'   or "prob" for probability-scale delta-method contrasts.
#' @param CatMethod Categorical test method: "auto" (default), "chisq", or "fisher".
#'
#' @return A `gtsummary::tbl_summary` object.
#'
#' @examples
#' \dontrun{
#' MakeComparisonTable(
#'   df_Clusters, "Cluster", vars_Predictors,
#'   Covariates = c("Age","Sex","wrat3"),
#'   AddEffectSize = TRUE, AddPairwise = TRUE, PairwiseMethod = "none"
#' )
#' }
#' @export
MakeComparisonTable <- function(
    DataFrame,
    CompVariable         = NULL,
    Variables,
    ...,
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
    BinaryPairwiseScale  = c("logit", "prob"),
    CatMethod            = "auto"
) {

  # Validate inputs ---------------------------------------------------------

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  extra_vars <- list(...)
  if (length(extra_vars) > 0) {
    extra_vars <- unlist(extra_vars, use.names = FALSE)
    if (!is.character(extra_vars)) stop("Unnamed arguments in ... must be character variable names.")
    Variables <- c(Variables, extra_vars)
  }

  if (!is.character(Variables)) stop("Variables must be a character vector.")
  Variables <- unique(Variables)

  if (is.null(PairwiseMethod)) PairwiseMethod <- "none"
  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale, choices = c("logit", "prob"))
  CatMethod <- match.arg(CatMethod, choices = c("auto", "chisq", "fisher"))

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

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  # Remove covariates from Variables (requested behavior) -------------------

  if (!is.null(Covariates) && length(Covariates)) {
    drop_from_vars <- intersect(Variables, Covariates)
    if (length(drop_from_vars)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_from_vars, collapse = ", "))
      Variables <- setdiff(Variables, Covariates)
    }
  }
  if (!length(Variables)) stop("No variables left in `Variables` after removing covariates.")

  # Helpers ----------------------------------------------------------------

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

  .fmt_p <- function(p, digits = 3) {
    ifelse(is.na(p), NA_character_, format.pval(p, digits = digits, eps = 10^-digits))
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
      chi <- tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL)
      if (!is.null(chi) && !is.null(chi$p.value)) return(list(p = as.numeric(chi$p.value), label = "Pearson chi-squared"))
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
          p = as.numeric(fish$p.value),
          label = if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"
        ))
      }
      return(list(p = NA_real_, label = "Fisher failed"))
    }

    chi_exp <- tryCatch(stats::chisq.test(tab, correct = FALSE)$expected, error = function(e) NULL)
    if (is.null(chi_exp)) return(.cat_global_test(tab, method = "fisher"))
    if (any(chi_exp < 5)) return(.cat_global_test(tab, method = "fisher"))
    .cat_global_test(tab, method = "chisq")
  }

  .cramers_v <- function(tab) {
    tab <- .clean_tab(tab)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)

    if (requireNamespace("DescTools", quietly = TRUE)) {
      return(as.numeric(DescTools::CramerV(tab, method = "bias.corrected")))
    }

    chi <- tryCatch(stats::chisq.test(tab, correct = FALSE)$statistic, error = function(e) NA_real_)
    if (is.na(chi)) return(NA_real_)
    n <- sum(tab)
    m <- min(nrow(tab), ncol(tab)) - 1
    if (m <= 0 || n <= 0) return(NA_real_)
    sqrt(as.numeric(chi) / (n * m))
  }

  .robust_group_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)

    cn <- names(stats::coef(fit_lm))
    idx <- grep(paste0("^", term), cn)
    if (!length(idx)) return(NA_real_)

    L <- matrix(0, nrow = length(idx), ncol = length(cn))
    L[cbind(seq_along(idx), idx)] <- 1

    lh <- tryCatch(car::linearHypothesis(fit_lm, hypothesis.matrix = L, vcov. = V, test = "F"),
                   error = function(e) NULL)
    if (is.null(lh) || !"Pr(>F)" %in% colnames(lh)) return(NA_real_)
    as.numeric(lh[2, "Pr(>F)"])
  }

  # Pick covariates per-variable: remove outcome, drop constants after CC ----
  .cov_use <- function(df_cc, covariates, outcome) {
    if (is.null(covariates) || !length(covariates)) return(character(0))
    covs <- setdiff(covariates, outcome)
    covs[vapply(covs, function(cv) {
      x <- df_cc[[cv]]
      length(unique(x[!is.na(x)])) >= 2
    }, logical(1))]
  }

  # Logistic pairwise contrasts (logit scale) --------------------------------
  .pairwise_logit_xb <- function(gm, df_cc, compvar, PairwiseMethod, Referent = NULL) {

    comp_levels <- levels(df_cc[[compvar]])
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% comp_levels) stop("Referent level not found: ", Referent)
      lapply(setdiff(comp_levels, Referent), function(x) c(Referent, x))
    } else {
      utils::combn(comp_levels, 2, simplify = FALSE)
    }

    beta <- stats::coef(gm)
    cn_beta <- names(beta)

    V0 <- tryCatch(stats::vcov(gm), error = function(e) NULL)
    if (is.null(V0)) return(NULL)

    tt <- stats::delete.response(stats::terms(gm))

    nd <- df_cc[rep(1, length(comp_levels)), , drop = FALSE]
    nd[[compvar]] <- factor(comp_levels, levels = comp_levels)

    X0 <- tryCatch(stats::model.matrix(tt, data = nd, xlev = gm$xlevels), error = function(e) NULL)
    if (is.null(X0)) return(NULL)

    X <- matrix(0, nrow = nrow(X0), ncol = length(beta), dimnames = list(NULL, cn_beta))
    commonX <- intersect(colnames(X0), cn_beta)
    X[, commonX] <- X0[, commonX, drop = FALSE]

    V <- matrix(0, nrow = length(beta), ncol = length(beta), dimnames = list(cn_beta, cn_beta))
    commonV <- intersect(rownames(V0), cn_beta)
    V[commonV, commonV] <- V0[commonV, commonV, drop = FALSE]

    idx <- setNames(seq_along(comp_levels), comp_levels)

    out <- purrr::map_dfr(combos, function(cp) {
      i <- idx[[cp[1]]]; j <- idx[[cp[2]]]
      d <- X[i, , drop = FALSE] - X[j, , drop = FALSE]
      est <- as.numeric(d %*% beta)
      se2 <- as.numeric(d %*% V %*% t(d))
      se <- sqrt(se2)

      p <- NA_real_
      if (is.finite(est) && is.finite(se) && se > 0) {
        z <- est / se
        p <- 2 * stats::pnorm(-abs(z))
      }

      tibble::tibble(
        contrast_label = .pair_label(cp[1], cp[2]),
        key = .pair_key(cp[1], cp[2]),
        p_val = p
      )
    })

    if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    out
  }

  # Prepare data ------------------------------------------------------------

  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- .as_factor(df[[CompVariable]])

  # Drop constant outcomes --------------------------------------------------

  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]

  drop <- setdiff(Variables, keep)
  if (length(drop)) warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))

  Variables <- keep
  if (!length(Variables)) stop("No variables left to summarise after dropping constants.")

  # Type detection ----------------------------------------------------------

  n_unique <- vapply(Variables, function(v) length(unique(df[[v]][!is.na(df[[v]])])), integer(1))
  is_num   <- vapply(Variables, function(v) is.numeric(df[[v]]), logical(1))
  is_dich  <- n_unique == 2

  treat_as_continuous <- is_num & !is_dich
  numeric_cont <- Variables[treat_as_continuous]
  dichotomous_numeric <- Variables[is_num & is_dich]

  type_list <- NULL
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
        vmap[[v]] <- if (any(vals == 1)) 1 else max(vals)
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

  # Statistic templates -----------------------------------------------------

  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  # Overall-only mode -------------------------------------------------------

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

    cap <- sprintf("Overall Summary (display: %s; continuous tests setting: %s)",
                   if (ParametricDisplay) "mean (SD)" else "median [IQR]",
                   if (Parametric) "parametric" else "non-parametric")
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  # Grouped summary ---------------------------------------------------------

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

  # P-values & test labels --------------------------------------------------

  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_, test_label = "Insufficient groups"))
    }

    if (is_cont) {
      k <- nlevels(df_vg[[CompVariable]])

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))

        cov_use <- .cov_use(df_cc, Covariates, outcome = var)

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3 || length(cov_use) == 0) {
          # If adjusted model isn't viable, fall back to unadjusted global tests
          if (Parametric) {
            if (k == 2) {
              p_un <- tryCatch(stats::t.test(.fmla(var, CompVariable), data = df_vg, var.equal = FALSE)$p.value,
                               error = function(e) NA_real_)
              return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_, test_label = "Welch t-test"))
            }
            p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
                             error = function(e) NA_real_)
            return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_, test_label = "ANOVA"))
          }

          p_un <- tryCatch(
            if (k == 2) stats::wilcox.test(.fmla(var, CompVariable), data = df_vg)$p.value
            else stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$p.value,
            error = function(e) NA_real_
          )
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                                test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, cov_use)), data = df_cc), error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_, test_label = "Model failed (adjusted)"))
        }

        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) as.numeric(a2[CompVariable, "Pr(>F)"]) else NA_real_
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = "ANCOVA (Type II)"))
        }

        p_rb <- .robust_group_p(fit, CompVariable)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_rb, test_label = "Robust ANCOVA (HC3 Wald)"))
      }

      # No covariates
      if (Parametric) {
        if (k == 2) {
          p_un <- tryCatch(stats::t.test(.fmla(var, CompVariable), data = df_vg, var.equal = FALSE)$p.value,
                           error = function(e) NA_real_)
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_, test_label = "Welch t-test"))
        }
        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_, test_label = "ANOVA"))
      }

      p_un <- tryCatch(
        if (k == 2) stats::wilcox.test(.fmla(var, CompVariable), data = df_vg)$p.value
        else stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$p.value,
        error = function(e) NA_real_
      )
      return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                            test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"))
    }

    # Categorical global tests
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)
    tst <- .cat_global_test(tab, method = CatMethod)

    tibble::tibble(variable = var, p_unadj = tst$p, p_adj = NA_real_, test_label = tst$label)
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::group_by(.data$variable) %>%
                                   dplyr::mutate(.is_main_row = dplyr::row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(
                                     p.value     = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value_fmt = dplyr::if_else(.data$.is_main_row, .fmt_p(.data$p.value, digits = pDigits), NA_character_),
                                     Test        = dplyr::if_else(.data$.is_main_row, .data$test_label, NA_character_)
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
          cov_use <- .cov_use(df_cc, Covariates, outcome = var)

          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3 || length(cov_use) == 0) {
            return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "partial η²"))
          }

          fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, cov_use)), data = df_cc), error = function(e) NULL)
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
            val <- tryCatch(abs(effectsize::cohens_d(.fmla(var, CompVariable), data = df_vg)$Cohens_d),
                            error = function(e) NA_real_)
            return(tibble::tibble(variable = var, effect_size = val, es_method = "|d|"))
          }
          val <- tryCatch(effectsize::eta_squared(stats::aov(.fmla(var, CompVariable), data = df_vg),
                                                  partial = FALSE)$Eta2[1],
                          error = function(e) NA_real_)
          return(tibble::tibble(variable = var, effect_size = val, es_method = "η²"))
        }

        n <- nrow(df_vg)
        H <- tryCatch(stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$statistic, error = function(e) NA_real_)
        eps2 <- suppressWarnings(as.numeric((H - k + 1) / (n - k)))
        return(tibble::tibble(variable = var, effect_size = eps2, es_method = "ε²"))
      }

      x <- .as_factor(df_vg[[var]])
      g <- .as_factor(df_vg[[CompVariable]])
      tab <- table(x, g)
      tibble::tibble(variable = var, effect_size = .cramers_v(tab), es_method = "Cramer's V")
    })

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     dplyr::left_join(es_df, by = "variable") %>%
                                     dplyr::mutate(
                                       effect_size = dplyr::if_else(.data$.is_main_row, .data$effect_size, NA_real_),
                                       ES_Method   = dplyr::if_else(.data$.is_main_row, .data$es_method, NA_character_)
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

      if (nrow(out_template) == 0) return(out_template[, c("variable", "contrast_label", "p_val")])

      # Continuous with covariates: emmeans, else fallback to unadjusted
      if (is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))

        cov_use <- .cov_use(df_cc, Covariates, outcome = var)

        # If adjusted model isn't viable, fall back to unadjusted pairwise (so you don't get blank columns)
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3 || length(cov_use) == 0) {

          df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
          df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))

          if (nlevels(df_vg[[CompVariable]]) < 2) {
            return(out_template[, c("variable", "contrast_label", "p_val")])
          }

          if (Parametric) {
            res <- tryCatch(
              stats::pairwise.t.test(
                x = df_vg[[var]],
                g = df_vg[[CompVariable]],
                p.adjust.method = PairwiseMethod,
                pool.sd = FALSE
              ),
              error = function(e) NULL
            )
            if (is.null(res) || is.null(res$p.value)) {
              return(out_template[, c("variable", "contrast_label", "p_val")])
            }

            r <- as.data.frame(as.table(res$p.value))
            r <- r[!is.na(r$Freq), , drop = FALSE]

            got <- purrr::map_dfr(combos, function(cp) {
              pv <- NA_real_
              hit1 <- r$Freq[r$Var1 == cp[2] & r$Var2 == cp[1]]
              hit2 <- r$Freq[r$Var1 == cp[1] & r$Var2 == cp[2]]
              if (length(hit1)) pv <- hit1[1]
              if (is.na(pv) && length(hit2)) pv <- hit2[1]
              tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = as.numeric(pv))
            })

            return(out_template %>%
                     dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                     dplyr::select("variable", "contrast_label", "p_val"))
          }

          out <- purrr::map_dfr(combos, function(cp) {
            sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
            sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
            p <- NA_real_
            if (nlevels(sub[[CompVariable]]) >= 2) {
              p <- tryCatch(stats::wilcox.test(.fmla(var, CompVariable), data = sub)$p.value, error = function(e) NA_real_)
            }
            tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = as.numeric(p))
          })

          if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)

          return(out_template %>%
                   dplyr::mutate(p_val = out$p_val[match(.data$key, out$key)]) %>%
                   dplyr::select("variable", "contrast_label", "p_val"))
        }

        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, cov_use)), data = df_cc), error = function(e) NULL)
        if (is.null(fit)) return(out_template[, c("variable", "contrast_label", "p_val")])

        emm <- tryCatch(emmeans::emmeans(fit, specs = CompVariable), error = function(e) NULL)
        if (is.null(emm)) return(out_template[, c("variable", "contrast_label", "p_val")])

        res <- NULL
        if (!is.null(Referent)) {
          ctr <- tryCatch(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1), error = function(e) NULL)
          if (!is.null(ctr)) {
            res <- tryCatch(as.data.frame(summary(ctr)), error = function(e) NULL)
            if (!is.null(res) && "p.value" %in% names(res) && !identical(PairwiseMethod, "none")) {
              res$p.value <- stats::p.adjust(res$p.value, method = PairwiseMethod)
            }
          }
        } else {
          ctr <- tryCatch(emmeans::pairs(emm, adjust = if (PairwiseMethod == "none") "none" else PairwiseMethod),
                          error = function(e) NULL)
          if (!is.null(ctr)) res <- tryCatch(as.data.frame(summary(ctr)), error = function(e) NULL)
        }

        if (is.null(res) || !("contrast" %in% names(res)) || !("p.value" %in% names(res))) {
          return(out_template[, c("variable", "contrast_label", "p_val")])
        }

        got <- tibble::tibble(
          key = vapply(as.character(res$contrast), function(s) {
            parts <- strsplit(.norm(s), "\\s*-\\s*")[[1]]
            if (length(parts) != 2) return(NA_character_)
            .pair_key(parts[1], parts[2])
          }, character(1)),
          p_val = as.numeric(res$p.value)
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Binary categorical with covariates: logistic pairwise (logit scale)
      if (!is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        df_cc[[var]] <- .as_factor(df_cc[[var]])

        cov_use <- .cov_use(df_cc, Covariates, outcome = var)

        if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5 && nlevels(df_cc[[var]]) == 2 && length(cov_use) > 0) {

          if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

          gm <- tryCatch(stats::glm(.fmla(var, c(CompVariable, cov_use)),
                                    data = df_cc, family = stats::binomial()),
                         error = function(e) NULL)

          if (!is.null(gm)) {
            got <- .pairwise_logit_xb(gm, df_cc, CompVariable, PairwiseMethod, Referent)
            if (!is.null(got) && nrow(got)) {
              return(out_template %>%
                       dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                       dplyr::select("variable", "contrast_label", "p_val"))
            }
          }
        }

        return(out_template[, c("variable", "contrast_label", "p_val")])
      }

      # Unadjusted continuous
      if (is_cont && is.null(Covariates)) {
        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) return(out_template[, c("variable", "contrast_label", "p_val")])

        if (Parametric) {
          res <- tryCatch(stats::pairwise.t.test(df_vg[[var]], df_vg[[CompVariable]],
                                                 p.adjust.method = PairwiseMethod, pool.sd = FALSE),
                          error = function(e) NULL)
          if (is.null(res) || is.null(res$p.value)) return(out_template[, c("variable", "contrast_label", "p_val")])

          r <- as.data.frame(as.table(res$p.value))
          r <- r[!is.na(r$Freq), , drop = FALSE]

          got <- purrr::map_dfr(combos, function(cp) {
            pv <- NA_real_
            hit1 <- r$Freq[r$Var1 == cp[2] & r$Var2 == cp[1]]
            hit2 <- r$Freq[r$Var1 == cp[1] & r$Var2 == cp[2]]
            if (length(hit1)) pv <- hit1[1]
            if (is.na(pv) && length(hit2)) pv <- hit2[1]
            tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = as.numeric(pv))
          })

          return(out_template %>%
                   dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                   dplyr::select("variable", "contrast_label", "p_val"))
        }

        out <- purrr::map_dfr(combos, function(cp) {
          sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
          p <- NA_real_
          if (nlevels(sub[[CompVariable]]) >= 2) {
            p <- tryCatch(stats::wilcox.test(.fmla(var, CompVariable), data = sub)$p.value, error = function(e) NA_real_)
          }
          tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = as.numeric(p))
        })

        if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)

        return(out_template %>%
                 dplyr::mutate(p_val = out$p_val[match(.data$key, out$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Unadjusted categorical pairwise
      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) return(out_template[, c("variable", "contrast_label", "p_val")])

      out <- purrr::map_dfr(combos, function(cp) {
        sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
        sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
        tab <- table(.as_factor(sub[[var]]), .as_factor(sub[[CompVariable]]))
        tst <- .cat_global_test(tab, method = CatMethod)
        tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = as.numeric(tst$p))
      })

      if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)

      out_template %>%
        dplyr::mutate(p_val = out$p_val[match(.data$key, out$key)]) %>%
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
                                           ~ dplyr::if_else(.data$.is_main_row, as.numeric(.), NA_real_)
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
              "Categorical method: ", CatMethod, ".",
              if (!is.null(Referent)) paste0(" Referent = ", Referent, ".") else "",
              if (!is.null(Covariates)) paste0(" Adjusted continuous uses emmeans; falls back to unadjusted if adjusted model is not viable.") else ""
            )
          )
      }
    }
  }

  cap <- sprintf(
    "Comparison Table (display: %s; continuous tests: %s). Pairwise: %s (%s). Categorical method: %s.",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    if (Parametric) "parametric" else "non-parametric (robust ANCOVA used when covariates present)",
    if (AddPairwise) "included" else "not included",
    PairwiseMethod,
    CatMethod
  )

  tbl %>% gtsummary::modify_caption(cap)
}
