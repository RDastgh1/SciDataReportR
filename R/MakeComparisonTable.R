#' Create a publication-ready comparison table using gtsummary
#'
#' Wraps \pkg{gtsummary} (\code{tbl_summary()} plus targeted post-processing) to produce a
#' predictable comparison table with (optionally) covariate-adjusted global tests,
#' covariate-adjusted and unadjusted pairwise contrasts, effect sizes, and a deterministic caption.
#' The function is label-friendly because the backbone is \pkg{gtsummary}, which respects variable
#' labels when present (e.g., via \pkg{labelled} attributes).
#'
#' @param DataFrame A data.frame or tibble containing all variables.
#' @param CompVariable Character scalar. Grouping variable name (e.g., "Cluster").
#'   If \code{NULL} or not present, the function returns an overall-only summary (no tests).
#' @param Variables Character vector of variable names to summarize and test.
#' @param ... Additional variable names (character) to include (appended to \code{Variables}).
#'   Unnamed arguments must be character variable names.
#' @param Covariates Optional character vector of covariate names used only for adjusted models.
#'   If any covariate is also present in \code{Variables}, it is dropped from display (with a warning).
#' @param ValueDigits Digits for continuous summary statistics.
#' @param pDigits Digits for formatted p-values (via \code{format.pval}).
#' @param AddEffectSize Logical. Add an "Effect size" column using \code{gtsummary::add_stat()}
#'   and the \pkg{effectsize} package.
#' @param EffectSizeDigits Digits for effect size formatting.
#' @param AddPairwise Logical. Add pairwise contrast p-value columns.
#' @param PairwiseMethod Character scalar. Multiple-comparisons adjustment method.
#'   Must be \code{"none"} or one of \code{stats::p.adjust.methods}.
#' @param Parametric Logical. Controls the *testing* mode for continuous outcomes.
#'   \itemize{
#'     \item \code{TRUE}: parametric tests and ANCOVA
#'     \item \code{FALSE}: nonparametric unadjusted tests, and robust ANCOVA (HC3) when adjusted
#'   }
#' @param ParametricDisplay Logical or \code{NULL}. Controls how continuous variables are displayed:
#'   mean (SD) when \code{TRUE}, median [Q1, Q3] when \code{FALSE}. If \code{NULL}, defaults to \code{Parametric}.
#' @param IncludeOverallN Logical. Add \code{gtsummary::add_n()} to show non-missing N per variable.
#' @param IncludeMissing Logical. Passed to \code{tbl_summary(missing=...)}. When \code{TRUE}, uses "ifany".
#' @param suppress_warnings Logical. If \code{TRUE}, suppresses warnings emitted by \pkg{gtsummary} calls.
#'   (Does not suppress warnings from your own tests unless explicitly noted below.)
#' @param Referent Optional character scalar level name for the grouping variable. If provided,
#'   pairwise contrasts become treatment-vs-control comparisons against \code{Referent}. If \code{NULL},
#'   all pairwise contrasts are computed.
#' @param IncludeOverallStats Logical. If \code{TRUE}, returns overall-only summary even if \code{CompVariable} is provided.
#' @param ShowPositiveBinaryOnLabel Logical. If \code{TRUE}, attempts to display the "positive" level for binary variables.
#' @param BinaryPairwiseScale Character. Reserved for future expansion (currently unused).
#' @param CatMethod Character. Global unadjusted categorical test method:
#'   \itemize{
#'     \item \code{"auto"} (default): chi-squared when assumptions are met; Fisher otherwise
#'     \item \code{"chisq"}: Pearson chi-squared
#'     \item \code{"fisher"}: Fisher exact for 2x2; Fisher with simulated p-value for larger tables
#'   }
#' @param MultiCatAdjusted Character. Adjusted global test for multi-category outcomes (3+ levels).
#'   Currently supports only \code{"multinomial_LR"} (default).
#'
#' @details
#' ## Summary statistics (display)
#' Continuous variables are summarized as mean (SD) when \code{ParametricDisplay=TRUE}, otherwise
#' median [Q1, Q3]. Categorical variables are summarized as n (p%).
#'
#' ## Global tests: continuous outcomes
#' Let \code{G} be the number of non-empty group levels after filtering.
#'
#' ### Unadjusted (no covariates)
#' \itemize{
#'   \item Parametric (\code{Parametric=TRUE}):
#'     \itemize{
#'       \item \code{G=2}: Welch t-test (\code{stats::t.test})
#'       \item \code{G>=3}: One-way ANOVA (\code{stats::aov})
#'     }
#'   \item Nonparametric (\code{Parametric=FALSE}):
#'     \itemize{
#'       \item \code{G=2}: Wilcoxon rank-sum test (\code{stats::wilcox.test})
#'       \item \code{G>=3}: Kruskal-Wallis test (\code{stats::kruskal.test})
#'     }
#' }
#'
#' ### Adjusted (with covariates)
#' "Adjusted" means: for each variable, a complete-case dataset is formed on
#' (outcome, group, covariates), and a regression model is fit with group and covariates.
#' The global p-value tests the group term while controlling for the covariates.
#'
#' \itemize{
#'   \item Parametric (\code{Parametric=TRUE}): ANCOVA via \code{lm}, Type II test of group term
#'     using \code{car::Anova(type=2)} (F test).
#'   \item Robust (\code{Parametric=FALSE}): Robust ANCOVA via \code{lm} with HC3 covariance
#'     (\code{sandwich::vcovHC(type="HC3")}), and a Type II Wald-style F test of the group term
#'     via \code{car::Anova(vcov.=..., test.statistic="F")}.
#' }
#'
#' ## Global tests: categorical outcomes
#' Unadjusted categorical tests use a contingency table of outcome by group after complete-case
#' filtering on (outcome, group).
#'
#' ### Unadjusted (no covariates)
#' \itemize{
#'   \item \code{CatMethod="chisq"}: Pearson chi-squared (\code{stats::chisq.test(correct=FALSE)}).
#'   \item \code{CatMethod="fisher"}:
#'     \itemize{
#'       \item 2x2: \code{stats::fisher.test}
#'       \item Larger tables: \code{stats::fisher.test(simulate.p.value=TRUE, B=1e4)}
#'     }
#'   \item \code{CatMethod="auto"} (default):
#'     \itemize{
#'       \item Compute chi-squared expected counts (warnings suppressed only in auto mode).
#'       \item Use Fisher if any expected count < 1 OR if >20% of expected counts are < 5.
#'       \item Otherwise use Pearson chi-squared.
#'     }
#' }
#'
#' ### Adjusted (with covariates)
#' Adjusted categorical tests are likelihood ratio (LR) tests comparing a full model
#' (group + covariates) to a reduced model (covariates only), using complete-case filtering on
#' (outcome, group, covariates).
#' \itemize{
#'   \item Binary outcome (2 levels): logistic regression via \code{glm(family=binomial)} and LR test
#'     for the group term via \code{drop1(test="Chisq")}.
#'   \item Multi-category outcome (3+ levels): multinomial logistic regression via \code{nnet::multinom}
#'     and LR test comparing reduced vs full via \code{anova(test="Chisq")}, with a log-likelihood
#'     difference fallback if needed. Default controlled by \code{MultiCatAdjusted="multinomial_LR"}.
#' }
#'
#' ## Pairwise contrasts
#' Pairwise contrasts are returned as additional \code{pw_...} columns on the main (label) row for each variable.
#'
#' ### Continuous outcomes, unadjusted (no covariates)
#' \itemize{
#'   \item If \code{Referent=NULL}:
#'     \itemize{
#'       \item Parametric: \code{stats::pairwise.t.test(pool.sd=FALSE)}
#'       \item Nonparametric: \code{stats::pairwise.wilcox.test}
#'     }
#'   \item If \code{Referent} is provided: run per-pair two-sample tests against the referent
#'     (Welch t-test or Wilcoxon), then apply \code{stats::p.adjust} unless \code{PairwiseMethod="none"}.
#' }
#'
#' ### Continuous outcomes, adjusted (with covariates)
#' Fit \code{lm(outcome ~ group + covariates)} on complete cases. Pairwise contrasts are computed using
#' \pkg{emmeans} on estimated marginal means for the group:
#' \itemize{
#'   \item If \code{Referent=NULL}: \code{emmeans::contrast(method="pairwise")}
#'   \item If \code{Referent} provided: \code{emmeans::contrast(method="trt.vs.ctrl")}
#' }
#' In robust mode (\code{Parametric=FALSE}), the HC3 covariance matrix is passed to \code{emmeans}
#' via \code{vcov.=...} so adjusted pairwise remains available.
#'
#' ### Categorical outcomes, adjusted (with covariates)
#' For each pairwise subset (e.g., group A vs group B), the outcome can collapse to fewer levels after
#' filtering. The engine detects this per pair:
#' \itemize{
#'   \item If outcome is binary within the subset: logistic regression LR test for group
#'   \item If outcome remains 3+ levels: multinomial LR test for group
#' }
#' This prevents missing pairwise results for variables that collapse within certain contrasts.
#'
#' ### Categorical outcomes, unadjusted (no covariates)
#' Not currently implemented (by design). If you need this, add it explicitly because the appropriate
#' choice depends heavily on sparse cells and whether you want chi-squared vs Fisher per pair.
#'
#' ## Complete-case filtering and dropped group levels
#' For adjusted tests and adjusted pairwise, complete-case filtering is performed per variable using
#' (outcome + group + covariates). This can drop one or more group levels for a given variable.
#' When this happens:
#' \itemize{
#'   \item The table includes a \code{Notes} column stating which group levels were dropped.
#'   \item Pairwise contrasts involving dropped levels are returned as \code{NA} (correct behavior).
#' }
#'
#' ## P-values
#' Numeric p-values are clamped to at least \code{.Machine$double.xmin} to avoid underflow yielding
#' literal 0. Formatted p-values use \code{format.pval} and may display values like \code{"<0.001"}.
#'
#' ## Effect sizes (optional)
#' If \code{AddEffectSize=TRUE}, an "Effect size" column is added using \pkg{effectsize}:
#' \itemize{
#'   \item Continuous, unadjusted:
#'     \itemize{
#'       \item 2 groups parametric: Hedges g (via \code{effectsize::cohens_d(..., hedges.correction=TRUE)})
#'       \item 2 groups nonparametric: rank-biserial r (\code{effectsize::rank_biserial})
#'       \item 3+ groups parametric: eta-squared (\code{effectsize::eta_squared})
#'       \item 3+ groups nonparametric: epsilon-squared (\code{effectsize::epsilon_squared})
#'     }
#'   \item Continuous, adjusted: partial eta-squared for the group term from the adjusted \code{lm}
#'     (\code{effectsize::eta_squared(partial=TRUE)}).
#'   \item Categorical, unadjusted: Cramer's V (\code{effectsize::cramers_v}).
#'   \item Binary categorical, adjusted: Nagelkerke R2 (\code{effectsize::r2_nagelkerke}).
#' }
#' Effect sizes are computed on complete-case data consistent with the adjusted model when covariates
#' are provided.
#'
#' ## References
#' \itemize{
#'   \item \pkg{gtsummary}: table engine and styling tools. \url{https://www.danieldsjoberg.com/gtsummary/}
#'   \item \pkg{emmeans}: estimated marginal means and contrasts (Lenth). \url{https://cran.r-project.org/package=emmeans}
#'   \item \pkg{car}: Type II tests (Fox and Weisberg). \url{https://cran.r-project.org/package=car}
#'   \item \pkg{sandwich}: HC3 robust covariance estimators (MacKinnon and White style). \url{https://cran.r-project.org/package=sandwich}
#'   \item \pkg{effectsize}: effect size computations (Ben-Shachar et al.). \url{https://cran.r-project.org/package=effectsize}
#'   \item \pkg{nnet}: multinomial regression via \code{multinom}. \url{https://cran.r-project.org/package=nnet}
#' }
#'
#' @return A \pkg{gtsummary} object (a \code{tbl_summary}) with additional columns:
#' \itemize{
#'   \item \code{p.value_fmt}: formatted p-value for the variable (main row only)
#'   \item \code{Test}: test label describing the global test used (main row only)
#'   \item \code{Notes}: notes about dropped group levels under complete-case filtering (main row only)
#'   \item Optional \code{effect_size}: effect size string (main row only) when \code{AddEffectSize=TRUE}
#'   \item Optional \code{pw_*}: pairwise p-value columns when \code{AddPairwise=TRUE}
#' }
#'
#' @examples
#' # Overall-only (no grouping)
#' # MakeComparisonTable(df, Variables = c("Age", "Sex"))
#'
#' # Unadjusted group comparison
#' # MakeComparisonTable(df, CompVariable="Group", Variables=c("Age","Sex"), AddPairwise=TRUE)
#'
#' # Adjusted robust ANCOVA + robust emmeans pairwise
#' # MakeComparisonTable(df, "Group", c("Outcome"), Covariates=c("Age","Sex"),
#' #                    Parametric=FALSE, AddPairwise=TRUE, AddEffectSize=TRUE)
#'
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
    MultiCatAdjusted     = c("multinomial_LR")
) {

  # ---- dependencies ----
  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "sandwich", "effectsize"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Please install: nnet (required for multinomial adjusted tests).")

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  extra_vars <- list(...)
  if (length(extra_vars) > 0) {
    extra_vars <- unlist(extra_vars, use.names = FALSE)
    if (!is.character(extra_vars)) stop("Unnamed arguments in ... must be character variable names.")
    Variables <- c(Variables, extra_vars)
  }
  if (!is.character(Variables)) stop("Variables must be a character vector.")
  Variables <- unique(Variables)

  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale)
  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)

  valid_methods <- c("none", stats::p.adjust.methods)
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1 || !(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "),
         ". Got: ", paste(PairwiseMethod, collapse = ", "))
  }

  comp_present <- !is.null(CompVariable) && is.character(CompVariable) &&
    length(CompVariable) == 1 && CompVariable %in% names(DataFrame)
  overall_mode <- isTRUE(IncludeOverallStats) || !comp_present

  if (!overall_mode && !CompVariable %in% names(DataFrame))
    stop("Grouping variable not found: ", CompVariable)

  if (!all(Variables %in% names(DataFrame)))
    stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse = ", "))

  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame)))
    stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1))
    stop("Referent must be a single character level name or NULL.")

  # ---- helpers ----
  # These helpers are local (not exported) because they isolate repeated edge-case handling:
  # p-value underflow, formula backticks, contingency-table validity, LR tests, and robust covariance wiring.

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

  .clamp_p <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    p[!is.finite(p)] <- NA_real_
    pmax(p, .Machine$double.xmin)
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

  .cat_global_test <- function(tab, method = c("auto","chisq","fisher")) {
    method <- match.arg(method)
    tab <- .clean_tab(tab)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(list(p = NA_real_, label = "Insufficient data"))

    if (method == "chisq") {
      chi <- tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL)
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

    # auto: suppress chisq warnings because we decide by expected counts anyway
    chi_obj <- tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)), error = function(e) NULL)
    if (is.null(chi_obj) || is.null(chi_obj$expected)) return(.cat_global_test(tab, method = "fisher"))

    exp <- chi_obj$expected
    prop_lt5 <- mean(exp < 5)

    if (any(exp < 1) || prop_lt5 > 0.20) return(.cat_global_test(tab, method = "fisher"))

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

      tibble::tibble(
        contrast_label = .pair_label(cp[1], cp[2]),
        key = .pair_key(cp[1], cp[2]),
        p_val = .clamp_p(p)
      )
    })

    if (!identical(PairwiseMethod, "none")) {
      out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    }
    out
  }

  .fmt_es_num <- function(x) {
    if (!is.finite(x)) return(NA_character_)
    formatC(x, format = "f", digits = EffectSizeDigits)
  }

  .pull_es <- function(df_es, patterns) {
    if (is.null(df_es) || !is.data.frame(df_es) || nrow(df_es) < 1) return(NA_real_)
    cn <- colnames(df_es)
    hit <- character(0)
    for (pat in patterns) {
      hit <- grep(pat, cn, value = TRUE)
      if (length(hit)) break
    }
    if (!length(hit)) return(NA_real_)
    suppressWarnings(as.numeric(df_es[[hit[1]]][1]))
  }

  .effect_size_string <- function(data, variable, by, Covariates, Parametric) {
    if (is.null(by) || !nzchar(by)) return(tibble::tibble(effect_size = NA_character_))

    cols <- c(variable, by, Covariates)
    dat <- data[stats::complete.cases(data[, cols, drop = FALSE]), , drop = FALSE]
    if (nrow(dat) < 3) return(tibble::tibble(effect_size = NA_character_))

    y <- dat[[variable]]
    g <- droplevels(.as_factor(dat[[by]]))
    if (nlevels(g) < 2) return(tibble::tibble(effect_size = NA_character_))

    ux <- unique(y[!is.na(y)])
    is_cont <- is.numeric(y) && length(ux) > 2

    if (is_cont) {
      if (!is.null(Covariates)) {
        fit <- tryCatch(stats::lm(.fmla(variable, c(by, Covariates)), data = dat), error = function(e) NULL)
        if (is.null(fit)) return(tibble::tibble(effect_size = NA_character_))
        es <- tryCatch(effectsize::eta_squared(fit, partial = TRUE), error = function(e) NULL)
        if (is.null(es) || !"Parameter" %in% names(es)) return(tibble::tibble(effect_size = NA_character_))
        row <- es[es$Parameter == by, , drop = FALSE]
        if (nrow(row) != 1) return(tibble::tibble(effect_size = NA_character_))
        val <- .pull_es(row, patterns = c("Eta2.*partial", "^Eta2"))
        if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
        return(tibble::tibble(effect_size = paste0("η²p=", .fmt_es_num(val))))
      }

      if (nlevels(g) == 2) {
        if (Parametric) {
          es <- tryCatch(effectsize::cohens_d(y ~ g, hedges.correction = TRUE), error = function(e) NULL)
          val <- .pull_es(es, patterns = c("^Cohens_d$", "Cohens_d"))
          if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
          return(tibble::tibble(effect_size = paste0("g=", .fmt_es_num(val))))
        } else {
          es <- tryCatch(effectsize::rank_biserial(y ~ g), error = function(e) NULL)
          val <- .pull_es(es, patterns = c("Rank_biserial"))
          if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
          return(tibble::tibble(effect_size = paste0("r_rb=", .fmt_es_num(val))))
        }
      }

      if (Parametric) {
        aov_obj <- tryCatch(stats::aov(.fmla(variable, by), data = dat), error = function(e) NULL)
        es <- tryCatch(effectsize::eta_squared(aov_obj), error = function(e) NULL)
        val <- .pull_es(es, patterns = c("^Eta2$", "Eta2"))
        if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
        return(tibble::tibble(effect_size = paste0("η²=", .fmt_es_num(val))))
      } else {
        kw <- tryCatch(stats::kruskal.test(.fmla(variable, by), data = dat), error = function(e) NULL)
        es <- tryCatch(effectsize::epsilon_squared(kw), error = function(e) NULL)
        val <- .pull_es(es, patterns = c("Epsilon2"))
        if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
        return(tibble::tibble(effect_size = paste0("ε²=", .fmt_es_num(val))))
      }
    }

    x <- droplevels(.as_factor(y))

    if (!is.null(Covariates) && nlevels(x) == 2) {
      fit <- tryCatch(stats::glm(.fmla(variable, c(by, Covariates)), data = dat, family = stats::binomial()),
                      error = function(e) NULL)
      if (is.null(fit)) return(tibble::tibble(effect_size = NA_character_))
      r2 <- tryCatch(effectsize::r2_nagelkerke(fit), error = function(e) NULL)
      val <- .pull_es(r2, patterns = c("^R2$", "R2"))
      if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
      return(tibble::tibble(effect_size = paste0("R²_N=", .fmt_es_num(val))))
    }

    tab <- table(x, g)
    es <- tryCatch(effectsize::cramers_v(tab), error = function(e) NULL)
    val <- .pull_es(es, patterns = c("Cramers_v"))
    if (!is.finite(val)) return(tibble::tibble(effect_size = NA_character_))
    tibble::tibble(effect_size = paste0("V=", .fmt_es_num(val)))
  }

  .pairwise_cont_unadjusted <- function(df_var, var, compvar, combos, Referent, Parametric, PairwiseMethod) {

    out_template <- purrr::map_dfr(combos, function(cp) {
      tibble::tibble(
        variable = var,
        contrast_label = .pair_label(cp[1], cp[2]),
        key = .pair_key(cp[1], cp[2]),
        p_val = NA_real_
      )
    })

    y <- df_var[[var]]
    g <- droplevels(.as_factor(df_var[[compvar]]))
    if (nlevels(g) < 2) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

    # referent: per-pair tests against referent
    if (!is.null(Referent)) {
      p_raw <- purrr::map_dbl(combos, function(cp) {
        sub <- df_var[df_var[[compvar]] %in% cp, , drop = FALSE]
        sub[[compvar]] <- droplevels(.as_factor(sub[[compvar]]))
        if (nlevels(sub[[compvar]]) < 2) return(NA_real_)

        p <- tryCatch(
          if (Parametric) stats::t.test(.fmla(var, compvar), data = sub)$p.value
          else stats::wilcox.test(.fmla(var, compvar), data = sub)$p.value,
          error = function(e) NA_real_
        )
        .clamp_p(p)
      })

      if (!identical(PairwiseMethod, "none")) {
        p_raw <- stats::p.adjust(p_raw, method = PairwiseMethod)
      }

      out_template$p_val <- p_raw
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    # all pairwise: use pairwise.*.test matrices
    if (Parametric) {
      pw <- tryCatch(
        stats::pairwise.t.test(x = y, g = g, pool.sd = FALSE, p.adjust.method = PairwiseMethod),
        error = function(e) NULL
      )
    } else {
      pw <- tryCatch(
        stats::pairwise.wilcox.test(x = y, g = g, p.adjust.method = PairwiseMethod),
        error = function(e) NULL
      )
    }

    if (is.null(pw) || is.null(pw$p.value)) {
      return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
    }

    mat <- pw$p.value
    got <- purrr::map_dfr(combos, function(cp) {
      a <- as.character(cp[1]); b <- as.character(cp[2])
      p <- NA_real_
      if (!is.null(rownames(mat)) && !is.null(colnames(mat))) {
        if (a %in% rownames(mat) && b %in% colnames(mat)) p <- mat[a, b]
        if (b %in% rownames(mat) && a %in% colnames(mat)) p <- mat[b, a]
      }
      tibble::tibble(key = .pair_key(a, b), p_val = .clamp_p(p))
    })

    out_template %>%
      dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
      dplyr::select("variable", "contrast_label", "p_val")
  }

  # ---- exclude covariates from Variables ----
  if (!is.null(Covariates)) {
    drop_cov <- intersect(Variables, Covariates)
    if (length(drop_cov)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov, collapse = ", "))
      Variables <- setdiff(Variables, drop_cov)
    }
  }

  # ---- data subset ----
  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- .as_factor(df[[CompVariable]])

  # ---- drop constants ----
  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]
  drop <- setdiff(Variables, keep)
  if (length(drop)) warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  Variables <- keep
  if (!length(Variables)) stop("No variables left to summarise after dropping constants.")

  # ---- type inference ----
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

  # ---- gtsummary display control ----
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

  stat_cont <- if (ParametricDisplay) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  # ---- overall only ----
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

    if (AddEffectSize) {
      warning("AddEffectSize ignored in overall-only mode (no CompVariable).")
    }

    if (suppress_warnings) tbl <- suppressWarnings(tbl)
    return(tbl)
  }

  # ---- grouped summary ----
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

  # ---- add effect sizes via add_stat() ----
  if (AddEffectSize) {
    tbl <- tbl %>%
      gtsummary::add_stat(
        fns = gtsummary::everything() ~ (function(data, variable, by, ...) {
          .effect_size_string(data, variable, by, Covariates = Covariates, Parametric = Parametric)
        }),
        location = gtsummary::everything() ~ "label"
      ) %>%
      gtsummary::modify_header(effect_size ~ "**Effect size**")
  }

  if (suppress_warnings) tbl <- suppressWarnings(tbl)
  if (nlevels(df[[CompVariable]]) < 2) return(tbl)

  full_group_levels <- levels(droplevels(.as_factor(df[[CompVariable]])))

  # ---- global p-values ----
  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = "Insufficient groups", Notes = NA_character_))
    }

    # continuous
    if (is_cont) {

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Insufficient data (adjusted)", Notes = note))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Model failed (adjusted)", Notes = note))
        }

        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        p_un <- .clamp_p(p_un)

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) .clamp_p(a2[CompVariable, "Pr(>F)"]) else NA_real_
          return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                                test_label = "ANCOVA (Type II)", Notes = note))
        }

        p_rb <- .robust_type2_p(fit, CompVariable)
        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_rb,
                              test_label = "Robust ANCOVA (HC3 Type II)", Notes = note))
      }

      # unadjusted
      k <- nlevels(df_vg[[CompVariable]])
      if (Parametric) {
        if (k == 2) {
          p_un <- tryCatch(stats::t.test(.fmla(var, CompVariable), data = df_vg)$p.value,
                           error = function(e) NA_real_)
          return(tibble::tibble(variable = var, p_unadj = .clamp_p(p_un), p_adj = NA_real_,
                                test_label = "Welch t-test", Notes = NA_character_))
        }
        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_vg))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        return(tibble::tibble(variable = var, p_unadj = .clamp_p(p_un), p_adj = NA_real_,
                              test_label = "ANOVA", Notes = NA_character_))
      }

      p_un <- tryCatch(
        if (k == 2) stats::wilcox.test(.fmla(var, CompVariable), data = df_vg)$p.value
        else stats::kruskal.test(.fmla(var, CompVariable), data = df_vg)$p.value,
        error = function(e) NA_real_
      )
      return(tibble::tibble(variable = var, p_unadj = .clamp_p(p_un), p_adj = NA_real_,
                            test_label = if (k == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis",
                            Notes = NA_character_))
    }

    # categorical
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- .cat_global_test(tab, method = CatMethod)
    p_un <- tst$p
    test_label <- tst$label
    p_adj <- NA_real_
    note <- NA_character_

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

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = test_label, Notes = note)
  })

  # ---- merge p-values into gtsummary table ----
  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::group_by(.data$variable) %>%
                                   dplyr::mutate(.is_main_row = dplyr::row_number() == 1) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(
                                     p.value     = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value     = dplyr::if_else(.data$.is_main_row, .clamp_p(.data$p.value), NA_real_),
                                     p.value_fmt = dplyr::if_else(.data$.is_main_row, .fmt_p(.data$p.value, digits = pDigits), NA_character_),
                                     Test        = dplyr::if_else(.data$.is_main_row, .data$test_label, NA_character_),
                                     Notes       = dplyr::if_else(.data$.is_main_row, .data$Notes, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**") %>%
    gtsummary::modify_header(Notes ~ "**Notes**")

  # ---- caption (deterministic metadata) ----
  cap_bits <- c()
  cap_bits <- c(cap_bits,
                if (ParametricDisplay) "Continuous summaries are mean (SD)" else "Continuous summaries are median [Q1, Q3]"
  )
  cap_bits <- c(cap_bits,
                if (Parametric) "Continuous global tests use Welch t-test (2 groups) or ANOVA (3+ groups)"
                else "Continuous global tests use Wilcoxon (2 groups) or Kruskal-Wallis (3+ groups)"
  )
  cap_bits <- c(cap_bits, paste0("Categorical global test: ", CatMethod))
  if (!is.null(Covariates)) {
    cap_bits <- c(cap_bits, paste0("Adjusted models include covariates: ", paste(Covariates, collapse = ", ")))
    cap_bits <- c(cap_bits,
                  if (Parametric) "Adjusted continuous test: ANCOVA (lm, Type II)"
                  else "Adjusted continuous test: robust ANCOVA (lm + HC3, Type II Wald F)"
    )
    cap_bits <- c(cap_bits, paste0("Adjusted multi-category categorical test: ", MultiCatAdjusted))
  }
  if (AddPairwise) {
    cap_bits <- c(cap_bits, paste0("Pairwise contrasts included; p-adjust method: ", PairwiseMethod))
    if (!is.null(Referent)) cap_bits <- c(cap_bits, paste0("Pairwise referent level: ", Referent))
  }
  if (AddEffectSize) cap_bits <- c(cap_bits, "Effect sizes computed via effectsize")
  tbl <- gtsummary::modify_caption(tbl, paste(cap_bits, collapse = ". "))

  # ---- pairwise ----
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

      # continuous pairwise: adjusted uses emmeans; unadjusted uses pairwise.* or per-pair tests
      if (is_cont && is.null(Covariates)) {
        df_var <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_var[[CompVariable]] <- droplevels(.as_factor(df_var[[CompVariable]]))
        got <- .pairwise_cont_unadjusted(
          df_var = df_var,
          var = var,
          compvar = CompVariable,
          combos = combos,
          Referent = Referent,
          Parametric = Parametric,
          PairwiseMethod = PairwiseMethod
        )
        return(got)
      }

      if (is_cont && !is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }
        if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
          df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

        V <- NULL
        if (!Parametric) V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)

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

        res <- tryCatch(
          as.data.frame(summary(ctr, adjust = if (PairwiseMethod == "none") "none" else PairwiseMethod)),
          error = function(e) NULL
        )
        if (is.null(res) || !("p.value" %in% names(res))) {
          res <- tryCatch(as.data.frame(summary(ctr)), error = function(e) NULL)
          if (is.null(res) || !("p.value" %in% names(res))) {
            return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
          }
          if (!identical(PairwiseMethod, "none")) res$p.value <- stats::p.adjust(res$p.value, method = PairwiseMethod)
        }

        got <- tibble::tibble(
          key = vapply(as.character(res$contrast), function(s) {
            parts <- strsplit(.norm(s), "\\s*-\\s*")[[1]]
            if (length(parts) != 2) return(NA_character_)
            a <- .strip_group_prefix(parts[1], CompVariable)
            b <- .strip_group_prefix(parts[2], CompVariable)
            .pair_key(a, b)
          }, character(1)),
          p_val = suppressWarnings(as.numeric(res$p.value))
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # categorical adjusted pairwise
      if (!is_cont && !is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))

        if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
          df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        }

        got <- .cat_pairwise_lr(df_cc, var, CompVariable, Covariates, combos, PairwiseMethod)
        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      dplyr::select(out_template, "variable", "contrast_label", "p_val")
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

  tbl
}
