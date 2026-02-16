#' Make a publication-ready comparison table (gtsummary backbone)
#'
#' Wraps `gtsummary::tbl_summary()` and appends deterministic, label-friendly columns for
#' global p-values (unadjusted and optionally covariate-adjusted), optional effect sizes,
#' and optional pairwise contrasts with multiple-comparisons adjustment and optional referent.
#'
#' The function is designed to be predictable: defaults are set for all arguments except
#' `DataFrame` and `CompVariable`, and outputs are stable across runs.
#'
#' ## Defaults and “just run it” behavior
#' - If `Variables` is `NULL` (default), all columns in `DataFrame` are summarized **except**
#'   `CompVariable` and any `Covariates`.
#' - If covariates appear in `Variables`, they are removed with a warning.
#' - If `PairwiseMethod` is `NULL`, `character(0)`, `""`, or `NA`, it defaults to `"bonferroni"`.
#'   This avoids failures when calling the function via `do.call()` with partially specified args.
#'
#' ## Continuous outcomes: global tests
#' **No covariates (`Covariates = NULL`)**
#' - Parametric (`Parametric = TRUE`):
#'   - 2 groups: Welch t-test (`stats::t.test`)
#'   - 3+ groups: ANOVA (`stats::aov`)
#' - Nonparametric (`Parametric = FALSE`):
#'   - 2 groups: Wilcoxon rank-sum (`stats::wilcox.test`)
#'   - 3+ groups: Kruskal-Wallis (`stats::kruskal.test`)
#'
#' **With covariates (`Covariates` provided)**
#' - Parametric (`Parametric = TRUE`): ANCOVA via `stats::lm` with Type II test for the group term
#'   using `car::Anova(type = 2)`.
#' - Robust (`Parametric = FALSE`): robust ANCOVA via `stats::lm` with HC3 covariance
#'   (`sandwich::vcovHC(type = "HC3")`) and a Type II Wald F test via `car::Anova(vcov. = V, test.statistic="F")`.
#'
#' Complete-case filtering is applied per variable:
#' - Unadjusted tests filter on `outcome + CompVariable`
#' - Adjusted tests filter on `outcome + CompVariable + Covariates`
#' This can drop whole group levels for a variable. Those dropped levels are reported in `Notes`
#' (if `IncludeNotes = TRUE`) and any pairwise contrasts involving the dropped level are `NA`.
#'
#' ## Continuous outcomes: pairwise tests
#' **No covariates**
#' - Parametric: `stats::pairwise.t.test(pool.sd = FALSE)` (Welch-style) with `p.adjust.method`
#' - Nonparametric: `stats::pairwise.wilcox.test` with `p.adjust.method`
#'
#' **With covariates**
#' - Always uses `emmeans` on the adjusted `lm` model.
#' - Robust mode passes the robust covariance matrix into `emmeans::emmeans(..., vcov. = V)`
#'   so adjusted pairwise still works in robust mode.
#'
#' Referent handling:
#' - If `Referent` is provided, uses treatment-vs-control (`emmeans::contrast(..., method="trt.vs.ctrl")`)
#' - Otherwise uses all pairwise (`method="pairwise"`)
#'
#' ## Categorical outcomes: global tests (unadjusted)
#' Controlled by `CatMethod`:
#' - `"chisq"`: Pearson chi-squared (`stats::chisq.test`)
#' - `"fisher"`: Fisher exact for 2x2; simulation-based Fisher for larger tables
#' - `"auto"` (default): chi-squared when assumptions are acceptable; otherwise Fisher
#'   (fallback triggered by small expected counts or failures).
#'
#' ## Categorical outcomes: adjusted tests (with covariates)
#' - Binary outcome: logistic regression likelihood ratio test for the group term
#'   using `stats::glm(..., family=binomial())` and `stats::drop1(test="Chisq")`.
#' - Multi-category outcome (3+ levels): multinomial LR test via `nnet::multinom`
#'   comparing reduced vs full model using `stats::anova(test="Chisq")` (with a logLik-based fallback).
#' - `MultiCatAdjusted` currently supports `"multinomial_LR"` only and defaults to it
#'   so users don’t accidentally fall into a missing-method pit.
#'
#' Adjusted pairwise for categorical outcomes:
#' - For each pairwise subset, the outcome can collapse (e.g., 3 levels overall becomes 2 levels in a pair).
#' - The function detects this per pair:
#'   - If collapsed to 2 levels: logistic LR
#'   - If still 3+ levels: multinomial LR
#'
#' ## Effect sizes (optional)
#' When `AddEffectSize = TRUE`, an `effect_size` column is added on the label row.
#' - Continuous: partial eta-squared (η²p) for the group term using `effectsize::eta_squared`
#'   (adjusted model if covariates exist; otherwise ANOVA model).
#' - Categorical: Cramer’s V from the unadjusted contingency table using `effectsize::cramers_v`.
#'
#' Effect sizes are intended as descriptive magnitudes; they do not change the hypothesis tests.
#'
#' ## P-values
#' To avoid numeric underflow returning exactly 0, numeric p-values are clamped to
#' `.Machine$double.xmin` before formatting.
#'
#' ## Caption
#' A deterministic caption is built with `gtsummary::modify_caption()` summarizing:
#' - continuous display (mean (SD) vs median [IQR])
#' - parametric vs robust mode
#' - covariate adjustment (and which covariates)
#' - categorical method (`CatMethod`)
#' - multi-category adjusted method (`MultiCatAdjusted`)
#' - whether pairwise was included and the p-adjust method.
#'
#' @param DataFrame A data frame.
#' @param CompVariable A single character name of the grouping variable.
#' @param Variables Character vector of variables to summarize. Default `NULL` means all columns
#'   except `CompVariable` and `Covariates`.
#' @param ... Optional unnamed character variable names to append to `Variables`.
#' @param Covariates Character vector of covariate names used only for adjusted tests. Default `NULL`.
#' @param ValueDigits Digits for continuous summary statistics. Default `2`.
#' @param pDigits Digits for formatted p-values. Default `3`.
#' @param AddEffectSize Logical. Add an effect size column. Default `FALSE`.
#' @param EffectSizeDigits Digits for effect size formatting. Default `2`.
#' @param AddPairwise Logical. Add pairwise p-value columns. Default `FALSE`.
#' @param PairwiseMethod Multiple comparisons adjustment. One of `"none"` or `stats::p.adjust.methods`.
#'   Default `"bonferroni"`. If passed as `NULL`/empty it falls back to `"bonferroni"`.
#' @param Parametric Logical. Use parametric tests (or robust ANCOVA when covariates exist) vs nonparametric.
#'   Default `TRUE`.
#' @param ParametricDisplay Logical. Controls how continuous variables are displayed:
#'   mean (SD) when `TRUE`, median [IQR] when `FALSE`. Default `NULL` which uses `Parametric`.
#' @param IncludeOverallN Logical. Add an N column via `gtsummary::add_n()`. Default `FALSE`.
#' @param IncludeMissing Logical. Include missing row in summaries. Default `FALSE`.
#' @param suppress_warnings Logical. Suppress warnings from `gtsummary::tbl_summary`. Default `FALSE`.
#' @param Referent Optional referent level for pairwise contrasts. Default `NULL` (all pairwise).
#' @param IncludeOverallStats Logical. If `TRUE` or `CompVariable` is `NULL`, return overall-only table.
#'   Default `FALSE`.
#' @param ShowPositiveBinaryOnLabel Logical. If `TRUE`, attempt to show “positive” level for dichotomous.
#'   Default `TRUE`.
#' @param BinaryPairwiseScale Not used in current implementation (reserved). Default `c("logit","prob")`.
#' @param CatMethod Global categorical method: `"auto"`, `"chisq"`, or `"fisher"`. Default `"auto"`.
#' @param MultiCatAdjusted Adjusted multi-category method. Default `"multinomial_LR"`.
#' @param IncludeNotes Logical. Include Notes column (dropped-group info). Default `TRUE`.
#'
#' @return A `gtsummary` object (`tbl_summary`) with added columns for p-values, tests,
#' optional effect sizes, optional pairwise p-values, and optional Notes.
#'
#' @examples
#' \dontrun{
#' tbl <- MakeComparisonTable(df, "Cluster")
#' tbl2 <- MakeComparisonTable(df, "Cluster", Covariates = c("Age","Sex"), AddPairwise = TRUE)
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

  # ---- dependencies ----
  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "sandwich", "effectsize"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Please install: nnet (required for multinomial adjusted tests).")

  # ---- validate inputs ----
  if (!is.data.frame(DataFrame)) stop("DataFrame must be a data.frame.")
  if (missing(CompVariable) || is.null(CompVariable) || !is.character(CompVariable) || length(CompVariable) != 1) {
    stop("CompVariable must be a single character column name.")
  }
  if (!CompVariable %in% names(DataFrame)) stop("Grouping variable not found: ", CompVariable)

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  # Variables default: everything except group and covariates
  if (is.null(Variables)) {
    Variables <- setdiff(names(DataFrame), c(CompVariable, Covariates))
  }

  # allow ... to append variable names
  extra_vars <- list(...)
  if (length(extra_vars) > 0) {
    extra_vars <- unlist(extra_vars, use.names = FALSE)
    if (!is.character(extra_vars)) stop("Unnamed arguments in ... must be character variable names.")
    Variables <- c(Variables, extra_vars)
  }
  if (!is.character(Variables)) stop("Variables must be a character vector or NULL.")
  Variables <- unique(Variables)

  if (!all(Variables %in% names(DataFrame))) {
    stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  }
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) {
    stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))
  }

  BinaryPairwiseScale <- match.arg(BinaryPairwiseScale)
  CatMethod <- match.arg(CatMethod)
  MultiCatAdjusted <- match.arg(MultiCatAdjusted)

  # PairwiseMethod: tolerate NULL/empty coming from do.call grids
  if (is.null(PairwiseMethod) || length(PairwiseMethod) == 0 ||
      (is.character(PairwiseMethod) && length(PairwiseMethod) == 1 && (is.na(PairwiseMethod) || PairwiseMethod == ""))) {
    PairwiseMethod <- "bonferroni"
  }
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1) {
    stop("PairwiseMethod must be a single character string.")
  }
  valid_methods <- unique(c("none", stats::p.adjust.methods))
  if (!(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "), ". Got: ", PairwiseMethod)
  }

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  # ---- helpers (used in multiple places to reduce repeated fragile code) ----
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

  # clamp p to avoid numeric 0 due to underflow (vectorized)
  .clamp_p <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    p[!is.finite(p)] <- NA_real_
    p <- pmax(p, .Machine$double.xmin, na.rm = FALSE)
    p
  }

  .fmt_p <- function(p, digits = 3) {
    ifelse(is.na(p), NA_character_, format.pval(p, digits = digits, eps = 10^-digits))
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
      chi <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL))
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

    # auto: use chi-squared only when expected counts look OK
    chi_obj <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL))
    if (is.null(chi_obj) || is.null(chi_obj$expected)) return(.cat_global_test(tab, method = "fisher"))
    if (any(chi_obj$expected < 5) || any(chi_obj$expected == 0)) return(.cat_global_test(tab, method = "fisher"))
    list(p = .clamp_p(chi_obj$p.value), label = "Pearson chi-squared")
  }

  .robust_type2_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)
    a2 <- tryCatch(car::Anova(fit_lm, type = 2, vcov. = V, test.statistic = "F"),
                   error = function(e) NULL)
    if (is.null(a2) || !(term %in% rownames(a2))) return(NA_real_)
    .clamp_p(a2[term, "Pr(>F)"])[1]
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
    .clamp_p(stats::pchisq(LR, df = ddf, lower.tail = FALSE))[1]
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

    list(p = .clamp_p(p)[1], label = "Multinomial LR")
  }

  .pw_matrix_to_long <- function(mat, combos) {
    # pairwise.*.test returns a lower-tri matrix with row/colnames
    # convert to named long keyed by sorted pair
    out <- purrr::map_dfr(combos, function(cp) {
      a <- cp[1]; b <- cp[2]
      pa <- NA_real_
      if (!is.null(mat)) {
        if (a %in% rownames(mat) && b %in% colnames(mat)) pa <- mat[a, b]
        if (b %in% rownames(mat) && a %in% colnames(mat)) pa <- mat[b, a]
      }
      tibble::tibble(key = .pair_key(a, b), p_val = as.numeric(pa))
    })
    out
  }

  .cat_pairwise_unadj <- function(df_vg, var, compvar, combos, method, PairwiseMethod) {
    out <- purrr::map_dfr(combos, function(cp) {
      sub <- df_vg[df_vg[[compvar]] %in% cp, , drop = FALSE]
      sub[[compvar]] <- droplevels(.as_factor(sub[[compvar]]))
      sub[[var]] <- droplevels(.as_factor(sub[[var]]))
      if (nlevels(sub[[compvar]]) < 2 || nlevels(sub[[var]]) < 2) {
        return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
      }
      tab <- table(sub[[var]], sub[[compvar]])
      tst <- .cat_global_test(tab, method = method)
      tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = tst$p)
    })
    if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    out
  }

  .cat_pairwise_lr <- function(df_cc, var, compvar, covs, combos, PairwiseMethod) {

    out <- purrr::map_dfr(combos, function(cp) {

      sub <- df_cc[df_cc[[compvar]] %in% cp, , drop = FALSE]
      sub[[compvar]] <- droplevels(.as_factor(sub[[compvar]]))
      sub[[var]] <- droplevels(.as_factor(sub[[var]]))

      if (nlevels(sub[[compvar]]) < 2 || nrow(sub) < 5 || nlevels(sub[[var]]) < 2) {
        return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
      }

      # binary logistic LR
      if (nlevels(sub[[var]]) == 2) {
        f_full <- .fmla(var, c(compvar, covs))
        gm <- tryCatch(stats::glm(f_full, data = sub, family = stats::binomial()), error = function(e) NULL)
        if (is.null(gm)) return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
        d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
        p <- NA_real_
        if (!is.null(d1) && compvar %in% rownames(d1) && "Pr(>Chi)" %in% colnames(d1)) {
          p <- as.numeric(d1[compvar, "Pr(>Chi)"])
        }
        return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = .clamp_p(p)[1]))
      }

      # multinomial LR
      f_full <- .fmla(var, c(compvar, covs))
      f_red  <- .fmla(var, covs)

      m_full <- tryCatch(nnet::multinom(f_full, data = sub, trace = FALSE), error = function(e) NULL)
      m_red  <- tryCatch(nnet::multinom(f_red,  data = sub, trace = FALSE), error = function(e) NULL)
      if (is.null(m_full) || is.null(m_red)) return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))

      a <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
      p <- .extract_anova_lr_p(a)
      if (is.na(p)) p <- .lrt_loglik_p(m_red, m_full)
      tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = .clamp_p(p)[1])
    })

    if (!identical(PairwiseMethod, "none")) out$p_val <- stats::p.adjust(out$p_val, method = PairwiseMethod)
    out
  }

  .fmt_es <- function(es_label, value, digits) {
    if (is.na(value) || !is.finite(value)) return(NA_character_)
    paste0(es_label, "=", formatC(value, digits = digits, format = "f"))
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
  cols <- unique(c(CompVariable, Variables, Covariates))
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  df[[CompVariable]] <- .as_factor(df[[CompVariable]])

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
  overall_mode <- isTRUE(IncludeOverallStats) || is.null(CompVariable)
  if (overall_mode) {
    tbl0 <- gtsummary::tbl_summary(
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
    if (IncludeOverallN) tbl0 <- tbl0 %>% gtsummary::add_n()
    if (suppress_warnings) tbl0 <- suppressWarnings(tbl0)
    return(tbl0)
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
  if (suppress_warnings) tbl <- suppressWarnings(tbl)
  if (nlevels(df[[CompVariable]]) < 2) return(tbl)

  full_group_levels <- levels(droplevels(.as_factor(df[[CompVariable]])))

  # ---- global tests + optional effect sizes ----
  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    # unadjusted CC for global
    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))

    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(
        variable = var, p_unadj = NA_real_, p_adj = NA_real_,
        test_label = "Insufficient groups", Notes = NA_character_, effect_size = NA_character_
      ))
    }

    p_un <- NA_real_
    p_adj <- NA_real_
    test_label <- NA_character_
    note <- NA_character_
    es <- NA_character_

    # continuous
    if (is_cont) {

      # effect size helper (continuous): eta^2_p for group term
      if (isTRUE(AddEffectSize)) {
        # unadjusted model (always available if 2+ groups)
        aov_un <- tryCatch(stats::aov(.fmla(var, CompVariable), data = df_vg), error = function(e) NULL)
        if (!is.null(aov_un)) {
          es_tab <- tryCatch(effectsize::eta_squared(aov_un, partial = TRUE), error = function(e) NULL)
          if (!is.null(es_tab) && "Eta2_partial" %in% names(es_tab) && any(grepl(CompVariable, es_tab$Parameter))) {
            v <- es_tab$Eta2_partial[which(grepl(CompVariable, es_tab$Parameter))[1]]
            es <- .fmt_es("\u03b7\u00b2p", as.numeric(v), EffectSizeDigits)
          }
        }
      }

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        note <- .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]]))

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Insufficient data (adjusted)", Notes = note, effect_size = es))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Model failed (adjusted)", Notes = note, effect_size = es))
        }

        # unadjusted p on same CC set (for transparency)
        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        p_un <- .clamp_p(p_un)[1]

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) .clamp_p(a2[CompVariable, "Pr(>F)"])[1] else NA_real_
          test_label <- "ANCOVA (Type II)"
        } else {
          p_adj <- .robust_type2_p(fit, CompVariable)
          test_label <- "Robust ANCOVA (HC3 Type II)"
        }

        # if effect sizes requested, prefer adjusted eta^2_p (group term)
        if (isTRUE(AddEffectSize)) {
          es_adj <- tryCatch(effectsize::eta_squared(fit, partial = TRUE, type = 2), error = function(e) NULL)
          if (!is.null(es_adj) && any(grepl(CompVariable, es_adj$Parameter))) {
            col_eta <- intersect(c("Eta2_partial", "Eta2_partial_CI_low"), names(es_adj))
            if ("Eta2_partial" %in% names(es_adj)) {
              v <- es_adj$Eta2_partial[which(grepl(CompVariable, es_adj$Parameter))[1]]
              es <- .fmt_es("\u03b7\u00b2p", as.numeric(v), EffectSizeDigits)
            }
          }
        }

        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                              test_label = test_label, Notes = note, effect_size = es))
      }

      # unadjusted continuous global
      k <- nlevels(df_vg[[CompVariable]])
      if (Parametric) {
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

      p_un <- .clamp_p(p_un)[1]
      return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = NA_real_,
                            test_label = test_label, Notes = NA_character_, effect_size = es))
    }

    # categorical
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    # effect size categorical: Cramer's V from unadjusted table
    if (isTRUE(AddEffectSize)) {
      v <- tryCatch(effectsize::cramers_v(tab)$Cramers_v, error = function(e) NA_real_)
      if (!is.na(v)) es <- .fmt_es("V", as.numeric(v), EffectSizeDigits)
    }

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

      # binary logistic LR
      if (nlevels(df_cc[[var]]) == 2 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5) {
        gm <- tryCatch(stats::glm(.fmla(var, c(CompVariable, Covariates)),
                                  data = df_cc, family = stats::binomial()),
                       error = function(e) NULL)
        if (!is.null(gm)) {
          d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
          if (!is.null(d1) && CompVariable %in% rownames(d1) && "Pr(>Chi)" %in% colnames(d1)) {
            p_adj <- .clamp_p(d1[CompVariable, "Pr(>Chi)"])[1]
            test_label <- "Logistic regression (LR)"
          }
        }
      }

      # multi-category multinomial LR
      if (nlevels(df_cc[[var]]) >= 3 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 10) {
        if (MultiCatAdjusted == "multinomial_LR") {
          g_lr <- .multinom_global_lr(df_cc, var, CompVariable, Covariates)
          if (!is.na(g_lr$p)) {
            p_adj <- .clamp_p(g_lr$p)[1]
            test_label <- g_lr$label
          }
        }
      }
    }

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                   test_label = test_label, Notes = note, effect_size = es)
  })

  # ---- merge global results into gtsummary table ----
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
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  if (isTRUE(AddEffectSize)) {
    tbl <- tbl %>% gtsummary::modify_header(effect_size ~ "**Effect size**")
  } else {
    # remove column if present from earlier joins
    tbl <- tbl %>% gtsummary::modify_table_body(~ dplyr::select(.x, -dplyr::any_of("effect_size")))
  }

  if (isTRUE(IncludeNotes)) {
    tbl <- tbl %>% gtsummary::modify_header(Notes ~ "**Notes**")
  } else {
    tbl <- tbl %>% gtsummary::modify_table_body(~ dplyr::select(.x, -dplyr::any_of("Notes")))
  }

  # Ensure Notes is last column if included
  tbl <- tbl %>% gtsummary::modify_table_body(~{
    tb <- .x
    if ("Notes" %in% names(tb)) tb <- dplyr::relocate(tb, .data$Notes, .after = dplyr::last_col())
    tb
  })

  # ---- pairwise ----
  if (isTRUE(AddPairwise) && nlevels(df[[CompVariable]]) > 1) {

    lvls <- levels(droplevels(.as_factor(df[[CompVariable]])))
    if (!is.null(Referent) && !Referent %in% lvls) stop("Referent level not found: ", Referent)

    combos <- if (!is.null(Referent)) {
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

      # Continuous, adjusted: emmeans
      if (is_cont && !is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))
        }
        if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
          df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc),
                        error = function(e) NULL)
        if (is.null(fit)) return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))

        V <- NULL
        if (!Parametric) V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)

        emm <- tryCatch(
          if (!is.null(V)) emmeans::emmeans(fit, specs = CompVariable, vcov. = V)
          else emmeans::emmeans(fit, specs = CompVariable),
          error = function(e) NULL
        )
        if (is.null(emm)) return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))

        ctr <- tryCatch(
          if (!is.null(Referent)) emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)
          else emmeans::contrast(emm, method = "pairwise"),
          error = function(e) NULL
        )
        if (is.null(ctr)) return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))

        # emmeans supports adjust, but "none" sometimes behaves inconsistently across versions; handle safely
        res <- tryCatch(as.data.frame(summary(ctr, adjust = if (PairwiseMethod == "none") "none" else PairwiseMethod)),
                        error = function(e) NULL)
        if (is.null(res) || !("p.value" %in% names(res))) {
          res <- tryCatch(as.data.frame(summary(ctr)), error = function(e) NULL)
          if (is.null(res) || !("p.value" %in% names(res))) {
            return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))
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
          p_val = as.numeric(res$p.value)
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "key", "p_val"))
      }

      # Continuous, unadjusted: pairwise.t.test / pairwise.wilcox.test
      if (is_cont && is.null(Covariates)) {
        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, "variable", "contrast_label", "key", "p_val"))
        }
        x <- df_vg[[var]]
        g <- df_vg[[CompVariable]]

        if (Parametric) {
          pw <- tryCatch(stats::pairwise.t.test(x, g, pool.sd = FALSE, p.adjust.method = PairwiseMethod)$p.value,
                         error = function(e) NULL)
        } else {
          pw <- tryCatch(stats::pairwise.wilcox.test(x, g, p.adjust.method = PairwiseMethod)$p.value,
                         error = function(e) NULL)
        }
        got <- .pw_matrix_to_long(pw, combos)

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "key", "p_val"))
      }

      # Categorical adjusted: logistic LR if binary in subset, else multinomial LR
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
                 dplyr::select("variable", "contrast_label", "key", "p_val"))
      }

      # Categorical unadjusted: per-pair chi-square/fisher using CatMethod
      if (!is_cont && is.null(Covariates)) {
        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        df_vg[[var]] <- droplevels(.as_factor(df_vg[[var]]))
        got <- .cat_pairwise_unadj(df_vg, var, CompVariable, combos, CatMethod, PairwiseMethod)

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "key", "p_val"))
      }

      dplyr::select(out_template, "variable", "contrast_label", "key", "p_val")
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

  # ---- caption (deterministic metadata string) ----
  cont_disp <- if (ParametricDisplay) "mean (SD)" else "median [IQR]"
  mode_txt <- if (Parametric) "parametric tests" else "robust or nonparametric tests"
  cov_txt <- if (is.null(Covariates)) "unadjusted" else paste0("adjusted for ", paste(Covariates, collapse = ", "))
  cat_txt <- paste0("categorical: ", CatMethod)
  pw_txt <- if (isTRUE(AddPairwise)) paste0("pairwise: ", PairwiseMethod, if (!is.null(Referent)) paste0(" (ref=", Referent, ")") else "") else "pairwise: none"
  cap <- paste0(
    "Continuous summaries shown as ", cont_disp, "; ", mode_txt, "; ", cov_txt, "; ",
    cat_txt, "; multi-category adjusted: ", MultiCatAdjusted, "; ", pw_txt, "."
  )
  tbl <- tbl %>% gtsummary::modify_caption(cap)

  tbl
}
