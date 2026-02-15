#' Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts
#'
#' Create a label friendly comparison table (via **gtsummary**) that summarizes variables
#' across groups, adds global tests, optional covariate adjustment, effect sizes, and
#' optional pairwise comparisons. Designed for “publication table” workflows where you
#' want the **Test** column to reflect the method that generated the p-values.
#'
#' This function is built to survive real world biomedical data, including:
#' non-syntactic column names (spaces, slashes, parentheses), sparse categorical data,
#' and models that occasionally refuse to cooperate.
#'
#' @details
#' ## Defaults (the ones you asked about)
#' - `CatMethod` defaults to `"auto"`.
#' - `MultiCatAdjusted` defaults to `"multinomial_LR"` (so multi-category covariate-adjusted
#'   globals work even when you do not specify it).
#'
#' ## What gets tested (global p-value)
#' The global p-value shown for each variable depends on variable type and options.
#'
#' ### Continuous outcomes
#' **No covariates**
#' - `Parametric = TRUE`
#'   - 2 groups: Welch two-sample t test (`stats::t.test`, `var.equal = FALSE`).
#'   - 3+ groups: one-way ANOVA (`stats::aov`).
#' - `Parametric = FALSE`
#'   - 2 groups: Wilcoxon rank-sum (`stats::wilcox.test`).
#'   - 3+ groups: Kruskal-Wallis (`stats::kruskal.test`).
#'
#' **With covariates**
#' - `Parametric = TRUE`: ANCOVA via linear model (`stats::lm`) with Type II test for the
#'   grouping term (`car::Anova(type = 2)`).
#' - `Parametric = FALSE`: Robust ANCOVA using `stats::lm` plus HC3 robust covariance
#'   (`sandwich::vcovHC(type = "HC3")`) and a Wald F test for the grouping term using
#'   `car::linearHypothesis`.
#'
#' ### Categorical outcomes
#' **No covariates**
#' - Default global test is Pearson chi-squared (`stats::chisq.test(correct = FALSE)`).
#' - If `CatMethod = "auto"`, Fisher’s exact is used when expected counts are small.
#'   For RxC tables, Fisher is computed via simulation (`simulate.p.value = TRUE`).
#'
#' **With covariates**
#' - Binary outcome: logistic regression LR global test (`stats::glm(..., binomial)`,
#'   then `stats::drop1(test="Chisq")` for the grouping term).
#' - Multi-category outcome:
#'   - If `MultiCatAdjusted = "multinomial_LR"` (default): multinomial logistic regression
#'     (`nnet::multinom`) with LR Chi-squared test for the grouping term via
#'     `car::Anova(type = 2)` when available. If the multinomial machinery fails, the
#'     function falls back to the unadjusted chi-squared or Fisher global test.
#'   - If `MultiCatAdjusted = "none"`: no adjusted global is attempted; unadjusted global is used.
#'
#' ## What gets tested (pairwise p-values)
#' Pairwise p-values are intended to match the method implied by the **Test** column,
#' but there is one deliberate compromise for multi-category categorical outcomes with
#' covariates:
#'
#' ### Continuous outcomes
#' - No covariates, `Parametric = TRUE`: `stats::pairwise.t.test(pool.sd = FALSE)`.
#' - No covariates, `Parametric = FALSE`: `stats::pairwise.wilcox.test`.
#' - With covariates:
#'   - `Parametric = TRUE`: pairwise comparisons of adjusted means using `emmeans::emmeans`
#'     and `emmeans::contrast` (pairwise or treatment-vs-control if `Referent` is set).
#'   - `Parametric = FALSE`: same as above, but using HC3 robust covariance in emmeans
#'     (`emmeans(..., vcov. = V_hc3)`), so pairwise columns populate under robust mode.
#'
#' ### Categorical outcomes
#' - No covariates: pairwise Pearson chi-squared or Fisher (depending on `CatMethod`)
#'   on 2-group subtables.
#' - With covariates:
#'   - Binary: pairwise Wald z tests of linear predictor contrasts from logistic regression.
#'   - Multi-category: **pairwise are unadjusted 2-group chi-squared or Fisher subtables**
#'     (Option 1). The global p-value can be adjusted (multinomial LR), but the pairwise
#'     comparisons are intentionally kept as simple contingency-table tests because a fully
#'     adjusted multinomial pairwise framework is more complex and far less stable in sparse
#'     biomedical tables.
#'
#' ## Effect sizes
#' - Continuous, 2 groups, unadjusted: |Cohen’s d| (`effectsize::cohens_d`).
#' - Continuous, 3+ groups, unadjusted parametric: eta-squared (`effectsize::eta_squared`).
#' - Continuous, adjusted: partial eta-squared from Type II ANOVA table.
#' - Continuous, nonparametric: epsilon-squared approximation from Kruskal-Wallis.
#' - Categorical: Cramer’s V.
#'
#' ## Why the backticks matter
#' This function supports non-syntactic column names (for example `"sCD163(ng/mL)"`,
#' `"PSQI TOTAL"`, `"CD4/CD8"`). Any internally constructed formulas always backtick
#' variable names so model calls do not silently fail.
#'
#' ## Practical interpretation of the p-values
#' - If you supply `Covariates`, the printed p-value column is **adjusted when available**
#'   (ANCOVA Type II for continuous parametric, robust Wald for continuous nonparametric,
#'   logistic LR for binary, multinomial LR for multi-category if enabled).
#' - If the adjusted model cannot be fit reliably (separation, rank deficiency, etc.),
#'   the function falls back to the unadjusted global test and labels it accordingly.
#'
#' ## References
#' - Agresti A. *Categorical Data Analysis* (chi-squared and Fisher testing).
#' - Fox J, Weisberg S. *car* package documentation for Type II/III tests (`car::Anova`,
#'   `car::linearHypothesis`).
#' - Lenth RV. *emmeans* package documentation (estimated marginal means and contrasts).
#' - Long JS, Ervin LH. (2000). Using heteroscedasticity-consistent standard errors in the
#'   linear regression model (HC3 rationale). See also `sandwich` package docs.
#' - Venables WN, Ripley BD. *Modern Applied Statistics with S* (multinomial logit via `nnet`).
#'
#' @param DataFrame A data frame.
#' @param CompVariable Grouping variable name (character scalar). If `NULL` or not present and
#'   `IncludeOverallStats = TRUE`, an overall-only table is returned.
#' @param Variables Character vector of variables to include. You may also pass additional
#'   variable names as unnamed arguments via `...` (convenience).
#' @param ... Optional additional variable names (character scalars). Useful when you accidentally
#'   typed `Variables="A", "B"` instead of `Variables=c("A","B")`.
#' @param Covariates Optional covariate names (character vector) used for adjusted models.
#'   Covariates are automatically removed from `Variables` so they do not appear as rows.
#' @param ValueDigits Digits for summary statistics (default 2).
#' @param pDigits Digits to display in formatted p-values. Formatting uses `format.pval`.
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
#' @param BinaryPairwiseScale For adjusted binary outcomes: "logit" (default) returns coefficient-contrast
#'   p-values on the logit scale; "prob" uses a margins-style probability scale approximation.
#' @param CatMethod Categorical test method: "chisq" (Pearson), "fisher", or "auto" (default).
#' @param MultiCatAdjusted For multi-category categorical outcomes with covariates:
#'   "multinomial_LR" (default) uses multinomial logistic LR global test when possible;
#'   "none" skips adjusted multi-category globals.
#'
#' @return A `gtsummary::tbl_summary` object (with added columns in `table_body`).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   Cluster = factor(sample(1:4, 200, TRUE)),
#'   `sCD163(ng/mL)` = rlnorm(200, 6, 0.5),
#'   `PSQI TOTAL` = rnorm(200, 8, 3),
#'   sex = factor(sample(c("F","M"), 200, TRUE)),
#'   age = rnorm(200, 50, 10)
#' )
#'
#' MakeComparisonTable(df, "Cluster", c("sCD163(ng/mL)", "PSQI TOTAL", "sex"),
#'   AddEffectSize = TRUE, AddPairwise = TRUE,
#'   PairwiseMethod = "none", Parametric = TRUE,
#'   Covariates = c("age", "sex")
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
    CatMethod            = c("auto", "chisq", "fisher"),
    MultiCatAdjusted     = c("multinomial_LR", "none")
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
  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale)
  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)

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

  # If covariates are in Variables, remove them so they do not appear as rows
  if (!is.null(Covariates)) {
    drop_cov_from_vars <- intersect(Variables, Covariates)
    if (length(drop_cov_from_vars)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov_from_vars, collapse = ", "))
      Variables <- setdiff(Variables, drop_cov_from_vars)
    }
    if (!length(Variables)) stop("No Variables left after removing covariates.")
  }

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
        return(list(p = as.numeric(fish$p.value),
                    label = if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"))
      }
      return(list(p = NA_real_, label = "Fisher failed"))
    }

    # auto: decide by expected counts
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

    # Guard: aliased coefficients -> linearHypothesis can error
    ali <- tryCatch(stats::alias(fit_lm)$Complete, error = function(e) NULL)
    if (!is.null(ali) && length(ali)) return(NA_real_)

    L <- matrix(0, nrow = length(idx), ncol = length(cn))
    L[cbind(seq_along(idx), idx)] <- 1

    lh <- tryCatch(
      car::linearHypothesis(fit_lm, hypothesis.matrix = L, vcov. = V, test = "F"),
      error = function(e) NULL
    )
    if (is.null(lh) || !"Pr(>F)" %in% colnames(lh)) return(NA_real_)
    as.numeric(lh[2, "Pr(>F)"])
  }

  .strip_group_prefix <- function(x, compvar) {
    x <- .norm(x)
    x <- gsub(paste0("^", compvar), "", x)
    x <- gsub("^\\s+|\\s+$", "", x)
    x
  }

  # Pairwise adjusted logistic (logit scale): Wald z tests of coefficient contrasts
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
      i <- idx[[cp[1]]]; j <- idx[[cp[2]]]
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

    if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    out
  }

  # Pairwise adjusted logistic (probability scale): margins-style
  .pairwise_logit_prob <- function(gm, df_cc, compvar, PairwiseMethod, Referent = NULL) {

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

    lvl_stats <- purrr::map(comp_levels, function(L) {
      nd <- df_cc
      nd[[compvar]] <- factor(L, levels = comp_levels)

      X <- tryCatch(stats::model.matrix(tt, data = nd, xlev = gm$xlevels), error = function(e) NULL)
      if (is.null(X)) return(list(pbar = NA_real_, grad = rep(NA_real_, length(beta))))

      eta <- as.numeric(X %*% beta)
      p <- stats::plogis(eta)
      pbar <- mean(p, na.rm = TRUE)

      w <- p * (1 - p)
      WX <- X * w
      grad <- colMeans(WX, na.rm = TRUE)

      list(pbar = pbar, grad = grad)
    })
    names(lvl_stats) <- comp_levels

    out <- purrr::map_dfr(combos, function(cp) {
      a <- cp[1]; b <- cp[2]
      da <- lvl_stats[[a]]; db <- lvl_stats[[b]]
      d <- da$pbar - db$pbar
      g <- da$grad - db$grad
      se <- sqrt(as.numeric(t(g) %*% V %*% g))
      z <- d / se
      pval <- 2 * stats::pnorm(-abs(z))

      tibble::tibble(
        contrast_label = .pair_label(a, b),
        key = .pair_key(a, b),
        p_val = pval
      )
    })

    if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    out
  }

  # Multinomial adjusted global for multi-category outcomes
  .multicat_global_multinom_lr <- function(df_cc, outcome, compvar, covars) {
    if (!requireNamespace("nnet", quietly = TRUE)) return(list(p = NA_real_, label = "Multinomial LR (nnet missing)"))
    y <- df_cc[[outcome]]
    y <- droplevels(.as_factor(y))
    if (nlevels(y) < 3) return(list(p = NA_real_, label = "Multinomial LR (not multi-category)"))

    f <- .fmla(outcome, c(compvar, covars))
    fit <- tryCatch(nnet::multinom(f, data = df_cc, trace = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(list(p = NA_real_, label = "Multinomial LR failed"))

    a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
    if (!is.null(a2) && compvar %in% rownames(a2) && "Pr(>Chisq)" %in% colnames(a2)) {
      return(list(p = as.numeric(a2[compvar, "Pr(>Chisq)"]), label = "Multinomial logistic (LR)"))
    }

    # Fallback: compare reduced vs full by LR if possible
    fit0 <- tryCatch(nnet::multinom(.fmla(outcome, covars), data = df_cc, trace = FALSE), error = function(e) NULL)
    if (is.null(fit0)) return(list(p = NA_real_, label = "Multinomial LR failed"))

    ll1 <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
    ll0 <- tryCatch(as.numeric(stats::logLik(fit0)), error = function(e) NA_real_)
    df1 <- tryCatch(attr(stats::logLik(fit), "df"), error = function(e) NA_integer_)
    df0 <- tryCatch(attr(stats::logLik(fit0), "df"), error = function(e) NA_integer_)

    if (is.na(ll1) || is.na(ll0) || is.na(df1) || is.na(df0) || (df1 <= df0)) {
      return(list(p = NA_real_, label = "Multinomial logistic (LR)"))
    }

    LR <- 2 * (ll1 - ll0)
    ddf <- df1 - df0
    p <- stats::pchisq(LR, df = ddf, lower.tail = FALSE)
    list(p = as.numeric(p), label = "Multinomial logistic (LR)")
  }

  # Prepare data ------------------------------------------------------------

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

  # Statistic templates -----------------------------------------------------

  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  # OVERALL-ONLY MODE -------------------------------------------------------

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

  # GROUPED SUMMARY ---------------------------------------------------------

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
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = "Insufficient groups"))
    }

    if (is_cont) {
      k <- nlevels(df_vg[[CompVariable]])

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Insufficient data (adjusted)"))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Model failed (adjusted)"))
        }

        p_un <- tryCatch(
          summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) as.numeric(a2[CompVariable, "Pr(>F)"]) else NA_real_
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = "ANCOVA (Type II)"))
        }

        p_rb <- .robust_group_p(fit, CompVariable)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_rb, test_label = "Robust ANCOVA (HC3 Wald)"))
      }

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
    p_un <- tst$p
    test_label <- tst$label
    p_adj <- NA_real_

    # Adjusted binary categorical: logistic LR global test
    if (!is.null(Covariates) && nlevels(droplevels(x)) == 2) {
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- .as_factor(df_cc[[var]])

      if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5 && nlevels(df_cc[[var]]) == 2) {
        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        gm <- tryCatch(stats::glm(.fmla(var, c(CompVariable, Covariates)),
                                  data = df_cc, family = stats::binomial()),
                       error = function(e) NULL)
        if (!is.null(gm)) {
          d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
          if (!is.null(d1) && CompVariable %in% rownames(d1)) {
            p_adj <- as.numeric(d1[CompVariable, "Pr(>Chi)"])
            test_label <- "Logistic regression (LR)"
          }
        }
      }
    }

    # Adjusted multi-category categorical: multinomial LR global (default enabled)
    if (!is.null(Covariates) && nlevels(droplevels(x)) >= 3 && MultiCatAdjusted == "multinomial_LR") {
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))
      if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 10 && nlevels(df_cc[[var]]) >= 3) {
        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        mc <- .multicat_global_multinom_lr(df_cc, var, CompVariable, Covariates)
        if (is.finite(mc$p)) {
          p_adj <- mc$p
          test_label <- mc$label
        }
      }
    }

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = test_label)
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
          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "partial η²"))
          }
          fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
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
        H <- tryCatch(stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$statistic,
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

      # Continuous with covariates: emmeans on lm (parametric or robust)
      if (is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        # Robust covariance for Parametric=FALSE
        V <- NULL
        if (!Parametric) {
          V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
        }

        emm <- tryCatch(
          if (!is.null(V)) emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", .btick(CompVariable))), vcov. = V)
          else emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", .btick(CompVariable)))),
          error = function(e) NULL
        )
        if (is.null(emm)) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        ctr <- tryCatch(
          if (!is.null(Referent)) emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)
          else emmeans::contrast(emm, method = "pairwise"),
          error = function(e) NULL
        )
        if (is.null(ctr)) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        res <- tryCatch(as.data.frame(summary(ctr, adjust = PairwiseMethod)), error = function(e) NULL)
        if (is.null(res) || !("p.value" %in% names(res))) {
          res <- tryCatch(as.data.frame(summary(ctr)), error = function(e) NULL)
          if (is.null(res) || !("p.value" %in% names(res))) {
            return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
          }
          if (!identical(PairwiseMethod, "none")) res$p.value <- stats::p.adjust(res$p.value, method = PairwiseMethod)
        }

        res$contrast <- as.character(res$contrast)

        got <- tibble::tibble(
          key = vapply(res$contrast, function(s) {
            parts <- strsplit(.norm(s), "\\s*-\\s*")[[1]]
            if (length(parts) != 2) return(NA_character_)
            a <- .strip_group_prefix(parts[1], CompVariable)
            b <- .strip_group_prefix(parts[2], CompVariable)
            .pair_key(a, b)
          }, character(1)),
          p_val = as.numeric(res$p.value)
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Binary categorical with covariates: logistic pairwise
      if (!is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        df_cc[[var]] <- .as_factor(df_cc[[var]])

        # Multi-category outcomes with covariates:
        # Option 1: do NOT attempt adjusted multinomial pairwise; use unadjusted 2-group subtables
        if (nlevels(df_cc[[var]]) >= 3) {
          df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
          df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
          out <- purrr::map_dfr(combos, function(cp) {
            sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
            sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
            tab <- table(.as_factor(sub[[var]]), .as_factor(sub[[CompVariable]]))
            tst <- .cat_global_test(tab, method = CatMethod)
            tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = as.numeric(tst$p))
          })
          if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
          return(out)
        }

        if (nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5 && nlevels(df_cc[[var]]) == 2) {

          if (!is.null(Referent)) df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)

          gm <- tryCatch(
            stats::glm(.fmla(var, c(CompVariable, Covariates)),
                       data = df_cc, family = stats::binomial()),
            error = function(e) NULL
          )

          if (!is.null(gm)) {
            got <- if (BinaryPairwiseScale == "logit") {
              .pairwise_logit_xb(gm, df_cc, CompVariable, PairwiseMethod, Referent)
            } else {
              .pairwise_logit_prob(gm, df_cc, CompVariable, PairwiseMethod, Referent)
            }

            if (!is.null(got) && nrow(got)) {
              return(out_template %>%
                       dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                       dplyr::select("variable", "contrast_label", "p_val"))
            }
          }
        }

        return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
      }

      # Continuous no covariates
      if (is_cont && is.null(Covariates)) {

        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(tibble::tibble(variable = character(), contrast_label = character(), p_val = numeric()))
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
            return(tibble::tibble(variable = character(), contrast_label = character(), p_val = numeric()))
          }

          r <- as.data.frame(as.table(res$p.value))
          r <- r[!is.na(r$Freq), , drop = FALSE]

          return(purrr::map_dfr(combos, function(cp) {
            lab <- .pair_label(cp[1], cp[2])
            pv <- NA_real_
            hit1 <- r$Freq[r$Var1 == cp[2] & r$Var2 == cp[1]]
            hit2 <- r$Freq[r$Var1 == cp[1] & r$Var2 == cp[2]]
            if (length(hit1)) pv <- hit1[1]
            if (is.na(pv) && length(hit2)) pv <- hit2[1]
            tibble::tibble(variable = var, contrast_label = lab, p_val = as.numeric(pv))
          }))
        }

        # nonparametric continuous
        out <- purrr::map_dfr(combos, function(cp) {
          sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
          if (nlevels(sub[[CompVariable]]) < 2) {
            return(tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = NA_real_))
          }
          p <- tryCatch(stats::wilcox.test(.fmla(var, CompVariable), data = sub)$p.value, error = function(e) NA_real_)
          tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = as.numeric(p))
        })

        if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
        return(out)
      }

      # Categorical no covariates: pairwise chi-squared/Fisher
      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(tibble::tibble(variable = character(), contrast_label = character(), p_val = numeric()))
      }

      out <- purrr::map_dfr(combos, function(cp) {
        sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
        sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
        tab <- table(.as_factor(sub[[var]]), .as_factor(sub[[CompVariable]]))
        tst <- .cat_global_test(tab, method = CatMethod)
        tibble::tibble(variable = var, contrast_label = .pair_label(cp[1], cp[2]), p_val = as.numeric(tst$p))
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
              if (!is.null(Covariates)) paste0(" Adjusted binary uses logistic contrasts on ", BinaryPairwiseScale, " scale.") else "",
              if (!is.null(Covariates)) " Multi-category categorical pairwise uses unadjusted 2-group subtables (Option 1)." else ""
            )
          )
      }
    }
  }

  cap <- sprintf(
    "Comparison Table (display: %s; continuous tests: %s). Pairwise: %s (%s). Categorical method: %s. Multi-cat adjusted: %s.",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    if (Parametric) "parametric" else "non-parametric (robust ANCOVA used when covariates present)",
    if (AddPairwise) "included" else "not included",
    PairwiseMethod,
    CatMethod,
    MultiCatAdjusted
  )

  tbl %>% gtsummary::modify_caption(cap)
}
