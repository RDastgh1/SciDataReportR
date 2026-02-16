#' Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts
#'
#' Create a label-friendly comparison table (via **gtsummary**) summarizing variables across
#' groups, with optional global hypothesis tests, covariate-adjusted tests, effect sizes,
#' and optional pairwise comparisons.
#'
#' This function is designed for publication tables where you want the **Test** column to
#' reflect the method used for the reported p-value(s), while still leveraging gtsummary's
#' formatting, labeling, and table structure.
#'
#' @details
#' ## Overview
#' `MakeComparisonTable()` is a wrapper around [gtsummary::tbl_summary()] that:
#' - creates group-wise descriptive summaries,
#' - computes global p-values (unadjusted or covariate-adjusted),
#' - optionally computes effect sizes,
#' - optionally computes pairwise p-values (with multiplicity control and optional referent),
#' - adds a **Notes** column to explain common edge cases (e.g., complete-case filtering dropping
#'   a group level, causing missing adjusted pairwise contrasts).
#'
#' The function aims to use gtsummary capabilities wherever possible (summary, formatting,
#' styling, captioning) and computes additional statistics not provided directly by gtsummary
#' (adjusted global tests, adjusted pairwise, effect sizes).
#'
#' ## Variable typing rules
#' - Continuous: numeric variables with > 2 unique non-missing values.
#' - Dichotomous numeric: numeric variables with exactly 2 unique non-missing values (treated as categorical).
#' - Categorical: factors, characters, logicals, and non-numeric variables (or dichotomous numeric).
#'
#' ## Global tests (p-values)
#' Global p-values are computed per variable, using complete-case data for the required fields.
#'
#' ### Continuous outcomes
#' - No covariates:
#'   - `Parametric = TRUE`:
#'     - 2 groups: Welch t-test ([stats::t.test()], unequal variances)
#'     - 3+ groups: one-way ANOVA ([stats::aov()])
#'   - `Parametric = FALSE`:
#'     - 2 groups: Wilcoxon rank-sum ([stats::wilcox.test()])
#'     - 3+ groups: Kruskal-Wallis ([stats::kruskal.test()])
#' - With covariates:
#'   - `Parametric = TRUE`: ANCOVA via linear model ([stats::lm()]) with Type II test for the
#'     grouping term ([car::Anova(type = 2)]).
#'   - `Parametric = FALSE`: robust ANCOVA via [stats::lm()] plus HC3 robust covariance
#'     ([sandwich::vcovHC(type = "HC3")]) and a Wald F-test for the grouping term
#'     ([car::linearHypothesis()]).
#'
#' ### Categorical outcomes (unadjusted)
#' - `CatMethod = "chisq"`: Pearson chi-squared with `correct = FALSE` ([stats::chisq.test()]).
#' - `CatMethod = "fisher"`: Fisher's exact test ([stats::fisher.test()]); for RxC tables,
#'   uses simulated p-value (B = 1e4).
#' - `CatMethod = "auto"` (default): uses chi-squared unless expected counts are small; then Fisher.
#'
#' ### Categorical outcomes (adjusted; covariates present)
#' - Binary outcome: logistic regression likelihood ratio test (LR) using
#'   [stats::glm(family = binomial)] and [stats::drop1(test = "Chisq")].
#' - Multi-category outcome (3+ levels): multinomial LR (default) via [nnet::multinom()] and
#'   an LR comparison of models with and without the grouping term.
#'   Controlled by `MultiCatAdjusted`, which defaults to `"multinomial_LR"`.
#'
#' ## Pairwise comparisons
#' Pairwise columns are added when `AddPairwise = TRUE`.
#'
#' ### Which pairs are compared
#' - If `Referent` is `NULL`: all pairwise group comparisons.
#' - If `Referent` is set: all groups are compared against the referent (treatment-vs-control).
#'
#' ### Multiplicity control
#' - `PairwiseMethod` defaults to `"bonferroni"` and may be any value in [stats::p.adjust.methods].
#' - Use `"none"` for no adjustment.
#'
#' ### Continuous outcomes
#' - No covariates:
#'   - Parametric: [stats::pairwise.t.test(pool.sd = FALSE)]
#'   - Nonparametric: [stats::pairwise.wilcox.test()]
#' - With covariates:
#'   - Uses adjusted means via [emmeans::emmeans()] on the ANCOVA model.
#'   - If `Parametric = FALSE`, the same ANCOVA model is used, but pairwise inference uses the
#'     HC3 robust covariance matrix supplied to emmeans (so adjusted pairwise still works).
#'
#' ### Categorical outcomes
#' - No covariates: pairwise chi-squared or Fisher (per `CatMethod`) on the 2-group subset.
#' - With covariates:
#'   - Binary outcome: for each pair, fit logistic models with and without group and use LR p-value.
#'   - Multi-category outcome: for each pair, fit multinomial models with and without group and use LR p-value.
#'   - Edge case: in a given pairwise subset, a multi-category outcome may collapse to 2 levels;
#'     in that case, the function automatically switches to logistic LR for that pair.
#'
#' ## Effect sizes
#' Effect sizes are provided when `AddEffectSize = TRUE`.
#'
#' - Continuous, 2 groups, unadjusted: absolute Cohen's d (|d|) via [effectsize::cohens_d()].
#' - Continuous, 3+ groups, unadjusted parametric: eta-squared (η²) via [effectsize::eta_squared()].
#' - Continuous, adjusted: partial eta-squared (partial η²) from Type II ANOVA table.
#' - Continuous, nonparametric: epsilon-squared approximation from Kruskal-Wallis.
#' - Categorical: Cramer's V (uses DescTools if available; otherwise a chi-squared-based approximation).
#'
#' ## Captions
#' The caption describes:
#' - display statistic for continuous variables (mean SD vs median IQR),
#' - whether covariates were used (and how continuous outcomes were tested),
#' - categorical test selection (`CatMethod`),
#' - adjusted multi-category method (`MultiCatAdjusted`),
#' - pairwise inclusion and multiplicity control (`PairwiseMethod`).
#'
#' @param DataFrame A data frame.
#' @param CompVariable Grouping variable name (character scalar). If NULL or missing and
#'   `IncludeOverallStats = TRUE`, an overall-only table is returned.
#' @param Variables Character vector of variables to include. You may also pass additional
#'   variable names as unnamed arguments via `...` (convenience).
#' @param ... Optional additional variable names (character scalars). Useful when you accidentally
#'   typed `Variables="A", "B"` instead of `Variables=c("A","B")`.
#' @param Covariates Optional covariate names (character vector) used for adjusted models.
#'   Any covariates appearing in `Variables` are removed from the displayed table with a warning.
#' @param ValueDigits Digits for summary statistics (default 2).
#' @param pDigits Digits to display formatted p-values using [format.pval()].
#' @param AddEffectSize Add effect sizes column (default FALSE).
#' @param EffectSizeDigits Digits for effect sizes (default 2).
#' @param AddPairwise Add pairwise p-value columns (default FALSE).
#' @param PairwiseMethod P-adjustment method (default "bonferroni"). Use "none" for no adjustment.
#'   Must be one of "none" or [stats::p.adjust.methods].
#' @param Parametric If TRUE, use parametric tests where relevant. If FALSE and covariates are present,
#'   robust ANCOVA is used for continuous outcomes and robust covariance is used for continuous pairwise.
#' @param ParametricDisplay If TRUE show mean (SD); if FALSE show median [IQR]. Defaults to `Parametric`.
#' @param IncludeOverallN If TRUE add overall N column.
#' @param IncludeMissing If TRUE include missing/unknown category rows in summaries.
#' @param suppress_warnings If TRUE suppress gtsummary warnings.
#' @param Referent Optional reference level for pairwise contrasts (character scalar). If set,
#'   pairwise comparisons are against this referent.
#' @param IncludeOverallStats If TRUE return overall-only summary (ignores grouping).
#' @param ShowPositiveBinaryOnLabel If TRUE, display only the “positive” level for binary categorical variables.
#' @param CatMethod Categorical test method: "chisq", "fisher", or "auto" (default).
#' @param MultiCatAdjusted Adjusted multi-category method when covariates are present:
#'   currently supports "multinomial_LR" (default) or "none".
#'
#' @return A `gtsummary::tbl_summary` object with added columns in `table_body`.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Cluster = factor(sample(0:3, 120, TRUE)),
#'   Age = rnorm(120, 50, 10),
#'   Sex = factor(sample(c("F","M"), 120, TRUE)),
#'   wrat3 = rnorm(120, 95, 12),
#'   Education_years = rnorm(120, 14, 2),
#'   Race = factor(sample(c("A","B","D"), 120, TRUE)),
#'   Hypertension = factor(sample(c("No","Yes"), 120, TRUE))
#' )
#'
#' MakeComparisonTable(
#'   df, "Cluster",
#'   Variables = c("Education_years", "Race", "Hypertension"),
#'   Covariates = c("Age","Sex","wrat3"),
#'   AddPairwise = TRUE,
#'   PairwiseMethod = "none",
#'   Parametric = FALSE
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

  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)

  if (is.null(PairwiseMethod)) PairwiseMethod <- "none"
  valid_methods <- c("none", stats::p.adjust.methods)
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1 || !(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "),
         ". Got: ", paste(PairwiseMethod, collapse = ", "))
  }

  req_pkgs <- c("gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
                "car", "emmeans", "effectsize", "sandwich")
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

  # Drop covariates from Variables (display) --------------------------------

  if (!is.null(Covariates)) {
    drop_cov_from_vars <- intersect(Variables, Covariates)
    if (length(drop_cov_from_vars)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov_from_vars, collapse = ", "))
      Variables <- setdiff(Variables, Covariates)
    }
  }

  if (!length(Variables)) stop("No variables left to summarise after removing covariates from Variables.")

  # Helpers ----------------------------------------------------------------

  as_factor_drop <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  btick <- function(x) paste0("`", gsub("`", "\\\\`", as.character(x)), "`")

  fmla <- function(lhs, rhs_terms) {
    rhs_terms <- as.character(rhs_terms)
    rhs <- paste(btick(rhs_terms), collapse = " + ")
    stats::as.formula(paste(btick(lhs), "~", rhs))
  }

  norm <- function(x) gsub("\\s+", " ", trimws(as.character(x)))
  pair_label <- function(a, b) paste(a, b, sep = " - ")
  pair_key <- function(a, b) paste(sort(c(norm(a), norm(b))), collapse = "||")

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
    if (is.null(chi_exp) || any(chi_exp < 5)) return(cat_global_test(tab, method = "fisher"))
    cat_global_test(tab, method = "chisq")
  }

  cramers_v <- function(tab) {
    tab <- clean_tab(tab)
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

  robust_group_p <- function(fit_lm, group_term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)
    cn <- names(stats::coef(fit_lm))
    idx <- grep(paste0("^", group_term), cn)
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

  multinom_lr_p <- function(df_cc, y, group, covs) {
    if (!requireNamespace("nnet", quietly = TRUE)) return(list(p = NA_real_, label = "Multinomial LR (nnet missing)"))

    df_cc[[y]] <- droplevels(as_factor_drop(df_cc[[y]]))
    if (nlevels(df_cc[[y]]) < 3) return(list(p = NA_real_, label = "Multinomial LR (collapsed)"))

    f_full <- fmla(y, c(group, covs))
    f_red  <- if (length(covs)) fmla(y, covs) else stats::as.formula(paste(btick(y), "~ 1"))

    m_full <- tryCatch(nnet::multinom(f_full, data = df_cc, trace = FALSE), error = function(e) NULL)
    m_red  <- tryCatch(nnet::multinom(f_red,  data = df_cc, trace = FALSE), error = function(e) NULL)
    if (is.null(m_full) || is.null(m_red)) return(list(p = NA_real_, label = "Multinomial LR (failed)"))

    lr <- 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))
    df_lr <- attr(stats::logLik(m_full), "df") - attr(stats::logLik(m_red), "df")
    if (!is.finite(lr) || !is.finite(df_lr) || df_lr <= 0) return(list(p = NA_real_, label = "Multinomial LR (invalid)"))

    list(p = stats::pchisq(lr, df = df_lr, lower.tail = FALSE), label = "Multinomial LR")
  }

  pairwise_cat_lr <- function(df_cc, y, group, covs, combos, PairwiseMethod) {

    out <- purrr::map_dfr(combos, function(cp) {
      a <- cp[1]; b <- cp[2]
      sub <- df_cc[df_cc[[group]] %in% c(a, b), , drop = FALSE]
      sub[[group]] <- droplevels(as_factor_drop(sub[[group]]))
      sub[[y]] <- droplevels(as_factor_drop(sub[[y]]))

      if (nlevels(sub[[group]]) < 2 || nrow(sub) < 5) {
        return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = NA_real_))
      }

      k <- nlevels(sub[[y]])
      if (k < 2) {
        return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = NA_real_))
      }

      if (k == 2) {
        f_full <- fmla(y, c(group, covs))
        f_red  <- if (length(covs)) fmla(y, covs) else stats::as.formula(paste(btick(y), "~ 1"))

        m_full <- tryCatch(stats::glm(f_full, data = sub, family = stats::binomial()), error = function(e) NULL)
        m_red  <- tryCatch(stats::glm(f_red,  data = sub, family = stats::binomial()), error = function(e) NULL)
        if (is.null(m_full) || is.null(m_red)) {
          return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = NA_real_))
        }

        an <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
        p <- if (!is.null(an) && "Pr(>Chi)" %in% names(an)) as.numeric(an$`Pr(>Chi)`[2]) else NA_real_
        return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = p))
      }

      if (!requireNamespace("nnet", quietly = TRUE)) {
        return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = NA_real_))
      }

      f_full <- fmla(y, c(group, covs))
      f_red  <- if (length(covs)) fmla(y, covs) else stats::as.formula(paste(btick(y), "~ 1"))

      m_full <- tryCatch(nnet::multinom(f_full, data = sub, trace = FALSE), error = function(e) NULL)
      m_red  <- tryCatch(nnet::multinom(f_red,  data = sub, trace = FALSE), error = function(e) NULL)
      if (is.null(m_full) || is.null(m_red)) {
        return(tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = NA_real_))
      }

      lr <- 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))
      df_lr <- attr(stats::logLik(m_full), "df") - attr(stats::logLik(m_red), "df")
      p <- if (is.finite(lr) && is.finite(df_lr) && df_lr > 0) stats::pchisq(lr, df = df_lr, lower.tail = FALSE) else NA_real_
      tibble::tibble(key = pair_key(a, b), contrast_label = pair_label(a, b), p_val = p)
    })

    if (!identical(PairwiseMethod, "none")) {
      out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    }
    out
  }

  # Prepare data ------------------------------------------------------------

  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- as_factor_drop(df[[CompVariable]])

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

    cap <- sprintf(
      "Overall summary (display: %s; categorical global test: %s).",
      if (ParametricDisplay) "mean (SD)" else "median [IQR]",
      CatMethod
    )
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
    cap <- sprintf("Overall summary (only one level present in '%s').", CompVariable)
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  # Pairwise combo definitions (global) -------------------------------------

  lvls_all <- levels(droplevels(as_factor_drop(df[[CompVariable]])))
  combos_all <- if (!is.null(Referent)) {
    if (!Referent %in% lvls_all) stop("Referent level not found: ", Referent)
    lapply(setdiff(lvls_all, Referent), function(x) c(Referent, x))
  } else {
    utils::combn(lvls_all, 2, simplify = FALSE)
  }

  # P-values, test labels, notes -------------------------------------------

  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])
    notes <- NA_character_

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))
    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = "Insufficient groups", Notes = "Insufficient group levels after filtering"))
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
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Insufficient data (adjusted)", Notes = notes))
        }

        fit <- tryCatch(stats::lm(fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Model failed (adjusted)", Notes = notes))
        }

        p_un <- tryCatch(
          summary(stats::aov(fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) as.numeric(a2[CompVariable, "Pr(>F)"]) else NA_real_
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                                test_label = "ANCOVA (Type II)", Notes = notes))
        }

        p_rb <- robust_group_p(fit, CompVariable)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_rb,
                              test_label = "Robust ANCOVA (HC3 Wald)", Notes = notes))
      }

      if (Parametric) {
        if (k == 2) {
          p_un <- tryCatch(stats::t.test(fmla(var, CompVariable), data = df_vg, var.equal = FALSE)$p.value,
                           error = function(e) NA_real_)
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                                test_label = "Welch t-test", Notes = notes))
        }
        p_un <- tryCatch(summary(stats::aov(fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                              test_label = "ANOVA", Notes = notes))
      }

      p_un <- tryCatch(
        if (k == 2) stats::wilcox.test(fmla(var, CompVariable), data = df_vg)$p.value
        else stats::kruskal.test(fmla(var, CompVariable), data = df_vg)$p.value,
        error = function(e) NA_real_
      )
      return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                            test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis",
                            Notes = notes))
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

      if (nlevels(df_cc[[var]]) == 2 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5) {
        gm <- tryCatch(stats::glm(fmla(var, c(CompVariable, Covariates)),
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

      if (is.na(p_adj) && nlevels(df_cc[[var]]) >= 3 && MultiCatAdjusted == "multinomial_LR") {
        mlr <- multinom_lr_p(df_cc, y = var, group = CompVariable, covs = Covariates)
        p_adj <- mlr$p
        test_label <- mlr$label
      }
    }

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                   test_label = test_label, Notes = notes)
  })

  # Apply p-values/test labels/notes into gtsummary table_body --------------

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
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
    gtsummary::modify_header(Notes ~ "**Notes**")

  tbl <- tbl %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows    = .data$.is_main_row & !is.na(.data$p.value) & .data$p.value <= 0.05,
      text_format = "bold"
    )

  # Effect sizes ------------------------------------------------------------

  if (AddEffectSize) {

    es_df <- purrr::map_dfr(Variables, function(var) {

      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "Insufficient groups"))
      }

      if (isTRUE(treat_as_continuous[var])) {

        if (!is.null(Covariates)) {
          cols_cc <- c(var, CompVariable, Covariates)
          df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
          df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))
          if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
            return(tibble::tibble(variable = var, effect_size = NA_real_, es_method = "partial η²"))
          }
          fit <- tryCatch(stats::lm(fmla(var, c(CompVariable, Covariates)), data = df_cc),
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
            val <- tryCatch(abs(effectsize::cohens_d(fmla(var, CompVariable), data = df_vg)$Cohens_d),
                            error = function(e) NA_real_)
            return(tibble::tibble(variable = var, effect_size = val, es_method = "|d|"))
          }
          val <- tryCatch(effectsize::eta_squared(stats::aov(fmla(var, CompVariable), data = df_vg),
                                                  partial = FALSE)$Eta2[1],
                          error = function(e) NA_real_)
          return(tibble::tibble(variable = var, effect_size = val, es_method = "η²"))
        }

        n <- nrow(df_vg)
        H <- tryCatch(stats::kruskal.test(fmla(var, CompVariable), data = df_vg)$statistic,
                      error = function(e) NA_real_)
        eps2 <- suppressWarnings(as.numeric((H - k + 1) / (n - k)))
        return(tibble::tibble(variable = var, effect_size = eps2, es_method = "ε²"))
      }

      x <- as_factor_drop(df_vg[[var]])
      g <- as_factor_drop(df_vg[[CompVariable]])
      tab <- table(x, g)
      v <- cramers_v(tab)
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

      # Continuous with covariates: emmeans on lm (robust vcov when Parametric=FALSE)
      if (is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))
        if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
          df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        }

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        fit <- tryCatch(stats::lm(fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        V <- NULL
        if (!Parametric) {
          V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
        }

        # IMPORTANT: emmeans specs must be a string or formula, not an rlang symbol.
        emm <- tryCatch(
          if (is.null(V)) emmeans::emmeans(fit, specs = CompVariable)
          else emmeans::emmeans(fit, specs = CompVariable, vcov. = V),
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
            parts <- strsplit(norm(s), "\\s*-\\s*")[[1]]
            if (length(parts) != 2) return(NA_character_)
            pair_key(parts[1], parts[2])
          }, character(1)),
          p_val = as.numeric(res$p.value)
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Binary/multi categorical with covariates: LR per pair subset
      if (!is_cont && !is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(as_factor_drop(df_cc[[CompVariable]]))
        df_cc[[var]] <- droplevels(as_factor_drop(df_cc[[var]]))

        if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
          df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        }

        got <- pairwise_cat_lr(df_cc, y = var, group = CompVariable, covs = Covariates,
                               combos = combos_all, PairwiseMethod = PairwiseMethod)

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Continuous no covariates
      if (is_cont && is.null(Covariates)) {

        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
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
            return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
          }

          r <- as.data.frame(as.table(res$p.value))
          r <- r[!is.na(r$Freq), , drop = FALSE]

          out <- purrr::map_dfr(combos_all, function(cp) {
            lab <- pair_label(cp[1], cp[2])
            pv <- NA_real_
            hit1 <- r$Freq[r$Var1 == cp[2] & r$Var2 == cp[1]]
            hit2 <- r$Freq[r$Var1 == cp[1] & r$Var2 == cp[2]]
            if (length(hit1)) pv <- hit1[1]
            if (is.na(pv) && length(hit2)) pv <- hit2[1]
            tibble::tibble(variable = var, contrast_label = lab, p_val = as.numeric(pv))
          })
          return(out)
        }

        out <- purrr::map_dfr(combos_all, function(cp) {
          sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(as_factor_drop(sub[[CompVariable]]))
          if (nlevels(sub[[CompVariable]]) < 2) {
            return(tibble::tibble(variable = var, contrast_label = pair_label(cp[1], cp[2]), p_val = NA_real_))
          }
          p <- tryCatch(stats::wilcox.test(fmla(var, CompVariable), data = sub)$p.value, error = function(e) NA_real_)
          tibble::tibble(variable = var, contrast_label = pair_label(cp[1], cp[2]), p_val = as.numeric(p))
        })

        if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
        return(out)
      }

      # Categorical no covariates: pairwise chi-squared/Fisher per CatMethod
      df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
      df_vg[[CompVariable]] <- droplevels(as_factor_drop(df_vg[[CompVariable]]))
      if (nlevels(df_vg[[CompVariable]]) < 2) {
        return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
      }

      out <- purrr::map_dfr(combos_all, function(cp) {
        sub <- df_vg[df_vg[[CompVariable]] %in% cp, , drop = FALSE]
        sub[[CompVariable]] <- droplevels(as_factor_drop(sub[[CompVariable]]))
        tab <- table(as_factor_drop(sub[[var]]), as_factor_drop(sub[[CompVariable]]))
        tst <- cat_global_test(tab, method = CatMethod)
        tibble::tibble(variable = var, contrast_label = pair_label(cp[1], cp[2]), p_val = as.numeric(tst$p))
      })

      if (!identical(PairwiseMethod, "none")) {
        out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
      }
      out
    })

    if (nrow(pw_long) > 0) {
      contrast_levels <- unique(pw_long$contrast_label)
      safe_names <- make.unique(paste0("pw_", make.names(contrast_levels)))
      map_cols <- tibble::tibble(contrast_label = contrast_levels, col_safe = safe_names)

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
        lab <- map_cols$contrast_label[match(col, map_cols$col_safe)]

        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(!!rlang::sym(col) ~ function(x) vapply(x, fmt_p, character(1), digits = pDigits)) %>%
          gtsummary::modify_table_styling(
            columns = col,
            rows    = .data$.is_main_row & !is.na(.data[[col]]) & .data[[col]] <= 0.05,
            text_format = "bold"
          ) %>%
          gtsummary::modify_header(!!rlang::sym(col) := paste0("**", lab, "**")) %>%
          gtsummary::modify_footnote(
            !!rlang::sym(col) ~ paste0(
              "Pairwise p-value adjustment: ", PairwiseMethod, ". ",
              "Categorical method: ", CatMethod, ".",
              if (!is.null(Referent)) paste0(" Referent = ", Referent, ".") else "",
              if (!is.null(Covariates)) paste0(
                " Adjusted continuous pairwise uses emmeans on ANCOVA (robust vcov if Parametric=FALSE). ",
                "Adjusted categorical uses LR per pair subset; MultiCatAdjusted = ", MultiCatAdjusted, "."
              ) else ""
            )
          )
      }
    }
  }

  # Caption ----------------------------------------------------------------

  cap <- sprintf(
    "Comparison table (display: %s). Global p-values: %s. Categorical global test: %s; adjusted multi-category: %s. Pairwise: %s (p-adjust: %s).",
    if (ParametricDisplay) "mean (SD)" else "median [IQR]",
    if (is.null(Covariates)) "unadjusted (no covariates)" else if (Parametric) "adjusted (ANCOVA Type II / LR)" else "adjusted (robust ANCOVA HC3 / LR)",
    CatMethod,
    MultiCatAdjusted,
    if (AddPairwise) "included" else "not included",
    PairwiseMethod
  )

  tbl %>% gtsummary::modify_caption(cap)
}
