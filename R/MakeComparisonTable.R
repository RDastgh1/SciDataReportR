#' Make a publication-ready comparison table (gtsummary backbone) with optional adjustment,
#' effect sizes, and pairwise contrasts
#'
#' Creates a label-friendly, predictable comparison table using \pkg{gtsummary} as the backbone
#' (\code{\link[gtsummary:tbl_summary]{gtsummary::tbl_summary()}}), and then adds deterministic
#' test metadata columns (global p-values, test names), optional effect sizes, optional pairwise
#' p-values, and an optional Notes column that documents complete-case group-level dropout.
#'
#' @param DataFrame A data frame.
#' @param CompVariable A single character string naming the grouping variable (2+ levels).
#'   If \code{NULL} or if \code{IncludeOverallStats = TRUE}, the function returns an overall-only
#'   table without tests/p-values.
#' @param Variables Character vector of variables to summarize and test.
#' @param ... Unnamed additional variable names to append to \code{Variables}.
#' @param Covariates Optional character vector of covariates used only for adjusted models
#'   (global adjusted p-values and adjusted pairwise p-values).
#'   If any covariate names appear in \code{Variables}, they are removed with a warning.
#' @param ValueDigits Digits for continuous summary values in the displayed table.
#' @param pDigits Digits used when formatting p-values via \code{format.pval()}.
#' @param AddEffectSize Logical; if \code{TRUE}, adds an \code{effect_size} column.
#'   Continuous: partial eta squared (η²p). Categorical: Cramér's V (V).
#' @param EffectSizeDigits Digits for effect sizes.
#' @param AddPairwise Logical; if \code{TRUE}, adds pairwise p-value columns \code{pw_*}.
#' @param PairwiseMethod Multiple-comparisons correction for pairwise p-values.
#'   Must be \code{"none"} or one of \code{stats::p.adjust.methods}.
#' @param Parametric Logical; if \code{TRUE}, uses parametric tests/models when applicable.
#'   If \code{FALSE}, uses nonparametric tests unadjusted and robust covariance for adjusted
#'   continuous models.
#' @param ParametricDisplay Logical; controls only how continuous summaries are displayed.
#'   If \code{NULL}, defaults to \code{Parametric}. \code{TRUE} shows mean (SD); \code{FALSE}
#'   shows median [Q1, Q3].
#' @param IncludeOverallN Logical; if \code{TRUE}, adds \code{gtsummary::add_n()}.
#' @param IncludeMissing Logical; passed to \code{gtsummary::tbl_summary(missing=)} as
#'   \code{"ifany"} when \code{TRUE}, otherwise \code{"no"}.
#' @param suppress_warnings Logical; if \code{TRUE}, suppresses warnings/messages from model fitting
#'   where reasonable (does not suppress your explicit warnings).
#' @param Referent Optional single character level of \code{CompVariable} to use as referent for
#'   treatment-versus-control pairwise contrasts. If \code{NULL}, all pairwise contrasts are computed.
#' @param IncludeOverallStats Logical; if \code{TRUE}, returns an overall-only table (no tests).
#' @param ShowPositiveBinaryOnLabel Logical; if \code{TRUE}, attempts to choose a "positive" level for
#'   dichotomous variables via \code{gtsummary::tbl_summary(value=)} so the row label reads as the
#'   positive category.
#' @param BinaryPairwiseScale Character; reserved for future scaling of binary pairwise results.
#'   Currently unused (kept for API stability).
#' @param CatMethod Categorical global test method: \code{"auto"}, \code{"chisq"}, or \code{"fisher"}.
#'   - \code{"chisq"}: Pearson chi-squared (no continuity correction).
#'   - \code{"fisher"}: Fisher exact for 2x2; simulated Fisher for larger tables.
#'   - \code{"auto"}: uses chi-squared when assumptions are acceptable (expected counts >= 5),
#'     otherwise Fisher.
#' @param MultiCatAdjusted Adjusted method for multi-category (3+ levels) outcomes when covariates
#'   are provided. Default is \code{"multinomial_LR"} (multinomial logistic likelihood ratio test).
#' @param IncludeNotes Logical; if \code{TRUE}, includes a Notes column and forces it to be the last
#'   column in the output. If \code{FALSE}, Notes is omitted.
#'
#' @return A \code{gtsummary} object (typically \code{tbl_summary}) with additional columns:
#'   \itemize{
#'     \item \code{effect_size} (optional)
#'     \item \code{p.value} (numeric, clamped to avoid underflow to 0)
#'     \item \code{p.value_fmt} (formatted character p-values)
#'     \item \code{Test} (global test label)
#'     \item \code{pw_*} (optional pairwise p-values per contrast)
#'     \item \code{Notes} (optional; documents dropped group levels under complete-case filtering)
#'   }
#'
#' @section Statistical tests (global):
#' Continuous outcomes (no covariates):
#' \itemize{
#'   \item Parametric: Welch t-test (2 groups) or ANOVA (3+ groups)
#'   \item Nonparametric: Wilcoxon rank-sum (2 groups) or Kruskal-Wallis (3+ groups)
#' }
#'
#' Continuous outcomes (with covariates):
#' \itemize{
#'   \item Parametric: ANCOVA via \code{lm(outcome ~ group + covariates)} and Type II F test for group term
#'     (\pkg{car} \code{Anova(type=2)}).
#'   \item Robust: same \code{lm} but with HC3 robust covariance (\pkg{sandwich} \code{vcovHC(type="HC3")})
#'     and Type II Wald F test for group term (\pkg{car} \code{Anova(vcov.=V, test.statistic="F")}).
#' }
#'
#' Categorical outcomes (no covariates):
#' \itemize{
#'   \item \code{CatMethod="chisq"}: Pearson chi-squared (no correction)
#'   \item \code{CatMethod="fisher"}: Fisher exact for 2x2, otherwise simulated Fisher (\code{B=10000})
#'   \item \code{CatMethod="auto"}: chi-squared if all expected counts >= 5, else Fisher
#' }
#'
#' Categorical outcomes (with covariates):
#' \itemize{
#'   \item Binary: logistic regression LR test for group term via \code{drop1(test="Chisq")}
#'   \item 3+ levels: multinomial LR test using \pkg{nnet} \code{multinom()} comparing reduced vs full model
#' }
#'
#' @section Pairwise tests:
#' Continuous outcomes (no covariates):
#' \itemize{
#'   \item Parametric: pairwise t-tests per contrast (Welch-style), then \code{p.adjust} unless \code{"none"}
#'   \item Nonparametric: pairwise Wilcoxon per contrast, then \code{p.adjust} unless \code{"none"}
#' }
#'
#' Continuous outcomes (with covariates):
#' \itemize{
#'   \item Uses \pkg{emmeans} on the adjusted \code{lm}. In robust mode, passes HC3 covariance to \code{emmeans}
#'     so adjusted pairwise results are still produced.
#' }
#'
#' Categorical outcomes (with covariates):
#' \itemize{
#'   \item For each pairwise subset, detects whether the outcome collapses to binary.
#'   \item If binary in that subset: logistic LR.
#'   \item If still 3+ levels: multinomial LR.
#'   \item P-values are adjusted across pairwise comparisons unless \code{PairwiseMethod="none"}.
#' }
#'
#' @section Complete-case filtering and dropped group levels:
#' For adjusted models (global adjusted tests and adjusted pairwise), complete-case filtering is performed
#' per variable on \code{outcome + group + covariates}. This can drop entire group levels for some variables.
#' Contrasts involving dropped levels return \code{NA}. When \code{IncludeNotes=TRUE}, the Notes column
#' documents which levels were dropped.
#'
#' @section P-values:
#' Numeric p-values are clamped to at least \code{.Machine$double.xmin} to avoid underflow to 0.
#' Formatted p-values use \code{format.pval()} and may display \code{"<0.001"}.
#'
#' @references
#' - Fox J, Weisberg S. An R Companion to Applied Regression (car::Anova Type II tests).
#' - Zeileis A. sandwich: Robust covariance matrix estimators (HC3).
#' - Lenth RV. emmeans: Estimated marginal means, aka least-squares means.
#' - Venables WN, Ripley BD. Modern Applied Statistics with S (nnet::multinom).
#'
#' @examples
#' \dontrun{
#' tbl <- MakeComparisonTable(
#'   DataFrame = df,
#'   CompVariable = "Cluster",
#'   Variables = c("Education_years", "Race"),
#'   Covariates = c("Age", "Sex"),
#'   AddPairwise = TRUE,
#'   Parametric = FALSE,
#'   PairwiseMethod = "holm",
#'   AddEffectSize = TRUE,
#'   IncludeNotes = TRUE
#' )
#' tbl %>% gtsummary::as_gt()
#' }
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
    MultiCatAdjusted     = c("multinomial_LR"),
    IncludeNotes         = TRUE
) {

  # Validate inputs

  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "sandwich", "effectsize"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Please install: nnet (required for multinomial adjusted tests).")
  }

  if (!is.data.frame(DataFrame)) stop("DataFrame must be a data.frame.")
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
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "), ".")
  }

  comp_present <- !is.null(CompVariable) &&
    is.character(CompVariable) &&
    length(CompVariable) == 1 &&
    CompVariable %in% names(DataFrame)

  overall_mode <- isTRUE(IncludeOverallStats) || !comp_present

  if (!overall_mode && !CompVariable %in% names(DataFrame)) {
    stop("Grouping variable not found: ", CompVariable)
  }
  if (!all(Variables %in% names(DataFrame))) {
    stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  }
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) {
    stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))
  }
  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  # Prepare data

  .as_factor <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  .norm <- function(x) gsub("\\s+", " ", trimws(as.character(x)))

  .clamp_p <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    p[!is.finite(p)] <- NA_real_
    pmax(p, .Machine$double.xmin, na.rm = FALSE)
  }

  .fmt_p <- function(p, digits = 3) {
    p <- .clamp_p(p)
    out <- rep(NA_character_, length(p))
    ok2 <- !is.na(p)
    out[ok2] <- format.pval(p[ok2], digits = digits, eps = 10^-digits)
    out
  }

  .btick <- function(x) paste0("`", gsub("`", "\\\\`", as.character(x)), "`")

  .fmla <- function(lhs, rhs_terms) {
    rhs_terms <- as.character(rhs_terms)
    rhs <- paste(.btick(rhs_terms), collapse = " + ")
    stats::as.formula(paste(.btick(lhs), "~", rhs))
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

    # auto: probe chi-squared silently, fall back to Fisher when expected counts < 5
    chi_obj <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL))
    if (is.null(chi_obj) || is.null(chi_obj$expected)) return(.cat_global_test(tab, method = "fisher"))
    if (any(chi_obj$expected < 5)) return(.cat_global_test(tab, method = "fisher"))
    list(p = .clamp_p(chi_obj$p.value), label = "Pearson chi-squared")
  }

  .robust_type2_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)
    a2 <- tryCatch(car::Anova(fit_lm, type = 2, vcov. = V, test.statistic = "F"), error = function(e) NULL)
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

  .pair_label <- function(a, b) paste(a, b, sep = " - ")
  .pair_key <- function(a, b) paste(sort(c(.norm(a), .norm(b))), collapse = "||")

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
    suppressWarnings(as.numeric(a[2, pcol[1]]))
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

  # Exclude covariates from Variables

  if (!is.null(Covariates)) {
    drop_cov <- intersect(Variables, Covariates)
    if (length(drop_cov)) {
      warning("Dropping covariate(s) from Variables: ", paste(drop_cov, collapse = ", "))
      Variables <- setdiff(Variables, drop_cov)
    }
  }

  # Subset columns early

  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- .as_factor(df[[CompVariable]])

  # Drop constants

  keep <- Variables[vapply(Variables, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]
  drop <- setdiff(Variables, keep)
  if (length(drop)) warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  Variables <- keep
  if (!length(Variables)) stop("No variables left to summarise after dropping constants.")

  # Type inference (numeric with >2 unique treated as continuous)

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

  # Label-friendly binary "positive" selection (gtsummary value=)

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

  # Build outputs

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

    if (AddEffectSize) warning("AddEffectSize ignored in overall-only mode (no CompVariable).")
    if (AddPairwise) warning("AddPairwise ignored in overall-only mode (no CompVariable).")

    return(tbl)
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

  if (nlevels(df[[CompVariable]]) < 2) return(tbl)

  full_group_levels <- levels(droplevels(.as_factor(df[[CompVariable]])))

  # Global p-values + test labels + Notes + effect sizes (per-variable)

  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))

    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(
        variable = var, p_unadj = NA_real_, p_adj = NA_real_,
        test_label = "Insufficient groups", Notes = NA_character_,
        effect_size = NA_character_
      ))
    }

    # defaults
    p_un <- NA_real_
    p_adj <- NA_real_
    test_label <- NA_character_
    note <- NA_character_
    es <- NA_character_

    # Continuous

    if (is_cont) {

      if (!is.null(Covariates)) {

        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))

        note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Insufficient data (adjusted)", Notes = note,
            effect_size = NA_character_
          ))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc), error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(
            variable = var, p_unadj = NA_real_, p_adj = NA_real_,
            test_label = "Model failed (adjusted)", Notes = note,
            effect_size = NA_character_
          ))
        }

        # unadjusted p on the same CC subset (context)
        p_un <- tryCatch(
          summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
          error = function(e) NA_real_
        )
        p_un <- .clamp_p(p_un)

        if (isTRUE(Parametric)) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          if (!is.null(a2) && CompVariable %in% rownames(a2)) {
            p_adj <- .clamp_p(a2[CompVariable, "Pr(>F)"])
          }
          test_label <- "ANCOVA (Type II)"
        } else {
          p_adj <- .robust_type2_p(fit, CompVariable)
          test_label <- "Robust ANCOVA (HC3 Type II)"
        }

        if (isTRUE(AddEffectSize)) {
          es_obj <- tryCatch(effectsize::eta_squared(fit, partial = TRUE), error = function(e) NULL)
          if (!is.null(es_obj) && "Parameter" %in% names(es_obj)) {
            hit <- es_obj[es_obj$Parameter == CompVariable, , drop = FALSE]
            if (nrow(hit) == 1 && "Eta2_partial" %in% names(hit)) {
              es <- paste0("η²p=", formatC(hit$Eta2_partial, digits = EffectSizeDigits, format = "f"))
            }
          }
        }

        return(tibble::tibble(
          variable = var, p_unadj = p_un, p_adj = p_adj,
          test_label = test_label, Notes = note,
          effect_size = es
        ))
      }

      # unadjusted continuous (no covariates)
      k <- nlevels(df_vg[[CompVariable]])

      if (isTRUE(Parametric)) {
        if (k == 2) {
          p_un <- tryCatch(stats::t.test(.fmla(var, CompVariable), data = df_vg)$p.value, error = function(e) NA_real_)
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
      p_un <- .clamp_p(p_un)

      if (isTRUE(AddEffectSize)) {
        if (k >= 3) {
          aov_fit <- tryCatch(stats::aov(.fmla(var, CompVariable), data = df_vg), error = function(e) NULL)
          es_obj <- tryCatch(effectsize::eta_squared(aov_fit, partial = TRUE), error = function(e) NULL)
          if (!is.null(es_obj) && "Parameter" %in% names(es_obj)) {
            hit <- es_obj[es_obj$Parameter == CompVariable, , drop = FALSE]
            if (nrow(hit) == 1 && "Eta2_partial" %in% names(hit)) {
              es <- paste0("η²p=", formatC(hit$Eta2_partial, digits = EffectSizeDigits, format = "f"))
            }
          }
        } else {
          # 2-group: Cohen's d (Welch-style interpretation, unpooled SD)
          es_obj <- tryCatch(effectsize::cohens_d(.fmla(var, CompVariable), data = df_vg, pooled_sd = FALSE),
                             error = function(e) NULL)
          if (!is.null(es_obj) && "Cohens_d" %in% names(es_obj)) {
            es <- paste0("d=", formatC(es_obj$Cohens_d[1], digits = EffectSizeDigits, format = "f"))
          }
        }
      }

      return(tibble::tibble(
        variable = var, p_unadj = p_un, p_adj = NA_real_,
        test_label = test_label, Notes = NA_character_,
        effect_size = es
      ))
    }

    # Categorical

    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- .cat_global_test(tab, method = CatMethod)
    p_un <- tst$p
    test_label <- tst$label

    if (isTRUE(AddEffectSize)) {
      es_obj <- tryCatch(effectsize::cramers_v(tab), error = function(e) NULL)
      if (!is.null(es_obj) && "Cramers_v" %in% names(es_obj)) {
        es <- paste0("V=", formatC(es_obj$Cramers_v[1], digits = EffectSizeDigits, format = "f"))
      }
    }

    if (!is.null(Covariates)) {

      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))

      note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

      if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
        df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
      }

      # binary logistic LR
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

      # multinomial LR
      if (nlevels(df_cc[[var]]) >= 3 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 10) {
        if (MultiCatAdjusted == "multinomial_LR") {
          g_lr <- .multinom_global_lr(df_cc, var, CompVariable, Covariates)
          if (!is.na(g_lr$p)) {
            p_adj <- .clamp_p(g_lr$p)
            test_label <- g_lr$label
          }
        }
      }

      # effect size on CC subset (consistent with adjusted tests)
      if (isTRUE(AddEffectSize)) {
        tab_cc <- table(.as_factor(df_cc[[var]]), .as_factor(df_cc[[CompVariable]]))
        es_obj <- tryCatch(effectsize::cramers_v(tab_cc), error = function(e) NULL)
        if (!is.null(es_obj) && "Cramers_v" %in% names(es_obj)) {
          es <- paste0("V=", formatC(es_obj$Cramers_v[1], digits = EffectSizeDigits, format = "f"))
        }
      }
    }

    tibble::tibble(
      variable = var, p_unadj = p_un, p_adj = p_adj,
      test_label = test_label, Notes = note,
      effect_size = es
    )
  })

  # Merge p-values/tests/notes/effect sizes into gtsummary

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
                                     effect_size = dplyr::if_else(.data$.is_main_row, .data$effect_size, NA_character_),
                                     Notes       = dplyr::if_else(.data$.is_main_row, .data$Notes, NA_character_)
                                   )
    )

  # Styling / headers

  if (isTRUE(AddEffectSize)) {
    tbl <- tbl %>% gtsummary::modify_header(effect_size ~ "**Effect size**")
  } else {
    # remove placeholder effect_size values if AddEffectSize=FALSE
    tbl <- tbl %>% gtsummary::modify_table_body(~ dplyr::select(.x, -any_of("effect_size")))
    if (!is.null(tbl$table_styling$header)) {
      tbl$table_styling$header <- dplyr::filter(tbl$table_styling$header, .data$column != "effect_size")
    }
  }

  tbl <- tbl %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  if (isTRUE(IncludeNotes)) {
    tbl <- tbl %>% gtsummary::modify_header(Notes ~ "**Notes**")
  }

  # Pairwise contrasts (optional)

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

      # Continuous pairwise WITHOUT covariates: per-combo tests + p.adjust
      if (is_cont && is.null(Covariates)) {
        got <- purrr::map_dfr(combos, function(cp) {
          sub <- df[df[[CompVariable]] %in% cp, , drop = FALSE]
          sub[[CompVariable]] <- droplevels(.as_factor(sub[[CompVariable]]))
          if (nlevels(sub[[CompVariable]]) < 2) {
            return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
          }
          p <- tryCatch(
            if (isTRUE(Parametric)) stats::t.test(.fmla(var, CompVariable), data = sub)$p.value
            else stats::wilcox.test(.fmla(var, CompVariable), data = sub)$p.value,
            error = function(e) NA_real_
          )
          tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = .clamp_p(p))
        })

        if (!identical(PairwiseMethod, "none")) {
          got$p_val <- stats::p.adjust(got$p_val, method = PairwiseMethod)
        }

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # Continuous adjusted: emmeans (robust passes HC3 vcov)
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

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc), error = function(e) NULL)
        if (is.null(fit)) return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))

        V <- NULL
        if (!isTRUE(Parametric)) V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)

        emm <- tryCatch(
          suppressMessages(
            if (!is.null(V)) emmeans::emmeans(fit, specs = CompVariable, vcov. = V)
            else emmeans::emmeans(fit, specs = CompVariable)
          ),
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
          suppressMessages(as.data.frame(summary(ctr, adjust = PairwiseMethod))),
          error = function(e) NULL
        )
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

      # Categorical adjusted: logistic-vs-multinomial per subset with collapse detection
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

  # Caption (deterministic metadata)

  cont_disp <- if (isTRUE(ParametricDisplay)) "mean (SD)" else "median [Q1, Q3]"
  cont_global <- if (isTRUE(Parametric)) {
    "Welch t-test (2 groups) or ANOVA (3+ groups)"
  } else {
    "Wilcoxon (2 groups) or Kruskal-Wallis (3+ groups)"
  }
  cat_global <- paste0("Categorical global test: ", CatMethod)

  adj_txt <- ""
  if (!is.null(Covariates)) {
    adj_txt <- paste0(
      " Adjusted models include covariates: ",
      paste(Covariates, collapse = ", "),
      ". Adjusted continuous test: ",
      if (isTRUE(Parametric)) "ANCOVA (lm, Type II)." else "robust ANCOVA (lm + HC3, Type II).",
      " Multi-category adjusted categorical test: ",
      MultiCatAdjusted, "."
    )
  }

  pw_txt <- ""
  if (isTRUE(AddPairwise)) {
    pw_txt <- paste0(
      " Pairwise contrasts included; p-adjust method: ", PairwiseMethod, ".",
      if (!is.null(Referent)) paste0(" Referent: ", Referent, ".") else ""
    )
  }

  es_txt <- ""
  if (isTRUE(AddEffectSize)) {
    es_txt <- " Effect sizes computed (continuous: η²p; categorical: V)."
  }

  tbl <- tbl %>%
    gtsummary::modify_caption(
      paste0(
        "Continuous summaries are ", cont_disp, ". ",
        "Continuous global tests use ", cont_global, ". ",
        cat_global, ".",
        adj_txt, pw_txt, es_txt
      )
    )

  # Notes column: optional and forced-last

  if (!("Notes" %in% names(tbl$table_body))) {
    return(tbl)
  }

  if (!isTRUE(IncludeNotes)) {
    tbl <- tbl %>%
      gtsummary::modify_table_body(~ dplyr::select(.x, -any_of("Notes")))
    if (!is.null(tbl$table_styling$header)) {
      tbl$table_styling$header <- dplyr::filter(tbl$table_styling$header, .data$column != "Notes")
    }
    return(tbl)
  }

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ dplyr::relocate(.x, .data$Notes, .after = dplyr::last_col()))

  if (!is.null(tbl$table_styling$header) && "Notes" %in% tbl$table_styling$header$column) {
    hdr <- tbl$table_styling$header
    tbl$table_styling$header <- dplyr::bind_rows(
      dplyr::filter(hdr, .data$column != "Notes"),
      dplyr::filter(hdr, .data$column == "Notes")
    )
  }

  tbl
}
