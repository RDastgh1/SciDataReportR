#' Make a publication-ready comparison table (gtsummary backbone)
#'
#' Builds a `gtsummary::tbl_summary()` table comparing groups defined by `CompVariable`,
#' and appends deterministic, label-friendly columns for global tests, optional adjusted tests,
#' optional effect sizes, optional pairwise contrasts, and an optional Notes column.
#'
#' The function is designed to be usable with just:
#' \preformatted{MakeComparisonTable(df, "Cluster")}
#'
#' @section Defaults (what happens if you do not set an argument):
#' - `Variables = NULL` (default): compare all columns except `CompVariable` and `Covariates`.
#' - `Covariates = NULL` (default): no covariate adjustment.
#' - `Parametric = TRUE` (default): parametric continuous tests (Welch t-test, ANOVA, ANCOVA).
#' - `ParametricDisplay = NULL` (default): if NULL, matches `Parametric` for display
#'   (mean (SD) if TRUE; median [IQR] if FALSE).
#' - `CatMethod = "auto"` (default): chi-squared when assumptions are adequate, otherwise Fisher.
#' - `MultiCatAdjusted = "multinomial_LR"` (default): multinomial LR for adjusted multi-category outcomes.
#' - `AddPairwise = FALSE` (default): no pairwise columns unless enabled.
#' - `PairwiseMethod = "bonferroni"` (default): multiplicity correction for pairwise p-values.
#'   If `PairwiseMethod` is `NULL`, `""`, `NA`, or `character(0)`, it falls back to `"bonferroni"`.
#' - `Referent = NULL` (default): all pairwise; if set, do treatment-vs-control vs referent.
#' - `AddEffectSize = FALSE` (default): effect sizes off unless enabled.
#' - `IncludeNotes = TRUE` (default): include Notes (set to FALSE to remove). Notes is forced last if included.
#'
#' @section Continuous outcomes: global tests:
#' **No covariates**
#' - `Parametric = TRUE`: 2 groups Welch t-test; 3+ groups ANOVA.
#' - `Parametric = FALSE`: 2 groups Wilcoxon rank-sum; 3+ groups Kruskal-Wallis.
#'
#' **With covariates**
#' - `Parametric = TRUE`: ANCOVA via `lm()` and Type II F test for group term via `car::Anova(type = 2)`.
#' - `Parametric = FALSE`: robust ANCOVA via `lm()` with HC3 covariance (`sandwich::vcovHC(type="HC3")`)
#'   and a Type II Wald F test via `car::Anova(vcov.=V, test.statistic="F")`.
#'
#' Complete-case filtering is per-variable:
#' - Unadjusted tests filter on `outcome + CompVariable`
#' - Adjusted tests filter on `outcome + CompVariable + Covariates`
#' This may drop entire group levels.
#'
#' @section Continuous outcomes: pairwise tests:
#' **No covariates**
#' - Parametric: `stats::pairwise.t.test(pool.sd = FALSE)`.
#' - Nonparametric: `stats::pairwise.wilcox.test`.
#'
#' **With covariates**
#' - Uses `emmeans` on the adjusted `lm()` model.
#' - Robust mode passes HC3 covariance into `emmeans::emmeans(..., vcov. = V)` so adjusted pairwise works.
#' - `Referent` set: treatment-vs-control vs referent; otherwise all pairwise.
#'
#' @section Categorical outcomes: global tests:
#' **Unadjusted**
#' - `CatMethod="chisq"`: Pearson chi-squared.
#' - `CatMethod="fisher"`: Fisher exact for 2x2; simulated Fisher for larger tables.
#' - `CatMethod="auto"`: chi-squared unless expected counts are small or zero, then Fisher.
#'
#' **Adjusted (with covariates)**
#' - Binary outcome: logistic regression likelihood ratio test (LR) for group term.
#' - 3+ outcome: multinomial LR (`nnet::multinom`) comparing reduced vs full model.
#'
#' @section Categorical outcomes: adjusted pairwise:
#' For each pairwise subset, outcome levels may collapse. The engine detects per pair:
#' - If outcome becomes binary: logistic LR.
#' - If outcome remains 3+ levels: multinomial LR.
#'
#' @section Notes column:
#' With covariates, complete-case filtering can drop group levels for specific variables.
#' Pairwise contrasts involving dropped levels are returned as `NA`. If `IncludeNotes=TRUE`,
#' Notes will describe which group levels were dropped. Notes is always the last column if included.
#'
#' @section P-values:
#' Numeric p-values are clamped to at least `.Machine$double.xmin` to avoid underflow returning 0.
#'
#' @section Caption:
#' Caption is deterministic and includes: display choice, parametric vs robust mode, covariates (if any),
#' categorical method, multi-category adjusted method, pairwise inclusion, and p-adjust method.
#'
#' @param DataFrame A data.frame or tibble.
#' @param CompVariable A single character name of the grouping variable.
#' @param Variables Character vector of variables to summarize. Default NULL compares all other columns.
#' @param ... Optional unnamed character variable names to append to Variables.
#' @param Covariates Character vector of covariates used only for adjusted tests. Default NULL.
#' @param ValueDigits Digits for continuous summaries. Default 2.
#' @param pDigits Digits for p-value formatting. Default 3.
#' @param AddEffectSize Logical. Add an effect size column. Default FALSE.
#' @param EffectSizeDigits Digits for effect size formatting. Default 2.
#' @param AddPairwise Logical. Add pairwise columns. Default FALSE.
#' @param PairwiseMethod P-adjust method for pairwise tests. Default "bonferroni".
#'   Must be one of `stats::p.adjust.methods` or "none". Empty inputs fall back to default.
#' @param Parametric Logical. Parametric tests when TRUE; robust or nonparametric when FALSE. Default TRUE.
#' @param ParametricDisplay Logical or NULL. If NULL, matches Parametric. Default NULL.
#' @param IncludeOverallN Logical. Add overall N via `gtsummary::add_n()`. Default FALSE.
#' @param IncludeMissing Logical. Include missing rows in summaries. Default FALSE.
#' @param suppress_warnings Logical. Suppress warnings from `tbl_summary()`. Default FALSE.
#' @param Referent Optional single referent level for pairwise. Default NULL.
#' @param IncludeOverallStats Logical. If TRUE, returns overall-only summary (ignores grouping). Default FALSE.
#' @param ShowPositiveBinaryOnLabel Logical. Attempt to set a “positive” level for dichotomous display. Default TRUE.
#' @param BinaryPairwiseScale Character. Reserved. Default "logit". Choices: "logit", "prob".
#' @param CatMethod Character. Default "auto". Choices: "auto", "chisq", "fisher".
#' @param MultiCatAdjusted Character. Default "multinomial_LR". Choices: "multinomial_LR".
#' @param IncludeNotes Logical. Include Notes column and force it last. Default TRUE.
#'
#' @return A `gtsummary` `tbl_summary` object with added columns in `table_body`.
#'
#' @references
#' - gtsummary: https://www.danieldsjoberg.com/gtsummary/
#' - emmeans: https://cran.r-project.org/package=emmeans
#' - car::Anova (Type II): https://cran.r-project.org/package=car
#' - sandwich::vcovHC (HC3): https://cran.r-project.org/package=sandwich
#' - effectsize: https://cran.r-project.org/package=effectsize
#' - nnet::multinom: https://cran.r-project.org/package=nnet
#'
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
    BinaryPairwiseScale = "logit",
    CatMethod = "auto",
    MultiCatAdjusted = "multinomial_LR",
    IncludeNotes = TRUE
) {

  # ---- dependencies ----
  req_pkgs <- c(
    "gtsummary", "dplyr", "tidyr", "purrr", "tibble", "rlang", "tidyselect",
    "car", "emmeans", "sandwich"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Please install: nnet (required for multinomial adjusted tests).")
  if (isTRUE(AddEffectSize) && !requireNamespace("effectsize", quietly = TRUE)) {
    stop("Please install: effectsize (required when AddEffectSize = TRUE).")
  }

  # ---- normalize possibly-empty “grid” inputs into real defaults ----
  .norm_choice <- function(x, default) {
    if (is.null(x) || length(x) == 0) return(default)
    x <- as.character(x)[1]
    if (is.na(x) || x == "") return(default)
    x
  }

  if (is.null(ParametricDisplay)) ParametricDisplay <- Parametric

  BinaryPairwiseScale <- match.arg(.norm_choice(BinaryPairwiseScale, "logit"), c("logit", "prob"))
  CatMethod <- match.arg(.norm_choice(CatMethod, "auto"), c("auto", "chisq", "fisher"))
  MultiCatAdjusted <- match.arg(.norm_choice(MultiCatAdjusted, "multinomial_LR"), c("multinomial_LR"))

  PairwiseMethod <- .norm_choice(PairwiseMethod, "bonferroni")
  valid_methods <- unique(c("none", stats::p.adjust.methods))
  if (!is.character(PairwiseMethod) || length(PairwiseMethod) != 1 || !(PairwiseMethod %in% valid_methods)) {
    stop("PairwiseMethod must be one of: ", paste(valid_methods, collapse = ", "), ". Got: ", PairwiseMethod)
  }

  # ---- validate core args ----
  if (!is.data.frame(DataFrame)) stop("DataFrame must be a data.frame.")
  if (!is.character(CompVariable) || length(CompVariable) != 1) stop("CompVariable must be a single character column name.")
  if (!CompVariable %in% names(DataFrame)) stop("Grouping variable not found: ", CompVariable)

  if (!is.null(Covariates)) {
    if (!is.character(Covariates)) stop("Covariates must be NULL or a character vector.")
    if (!all(Covariates %in% names(DataFrame))) {
      stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))
    }
  }

  if (!is.null(Referent) && (!is.character(Referent) || length(Referent) != 1)) {
    stop("Referent must be a single character level name or NULL.")
  }

  # ---- Variables default: all except CompVariable and Covariates ----
  if (is.null(Variables)) {
    Variables <- setdiff(names(DataFrame), c(CompVariable, Covariates))
  }

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

  # ---- helpers ----
  .as_factor <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)))
    factor(x)
  }

  .norm <- function(x) gsub("\\s+", " ", trimws(as.character(x)))
  .pair_label <- function(a, b) paste(a, b, sep = " - ")
  .pair_key <- function(a, b) paste(sort(c(.norm(a), .norm(b))), collapse = "||")

  .btick <- function(x) {
    x <- gsub("`", "\\\\`", as.character(x))
    paste0("`", x, "`")
  }

  .fmla <- function(lhs, rhs_terms) {
    rhs_terms <- as.character(rhs_terms)
    rhs <- paste(.btick(rhs_terms), collapse = " + ")
    stats::as.formula(paste(.btick(lhs), "~", rhs))
  }

  # vectorized clamp p (fixes the earlier length > 1 crash)
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
      if (!is.null(chi) && !is.null(chi$p.value)) return(list(p = .clamp_p(chi$p.value)[1], label = "Pearson chi-squared"))
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
          p = .clamp_p(fish$p.value)[1],
          label = if (nrow(tab) == 2 && ncol(tab) == 2) "Fisher exact" else "Fisher (sim.)"
        ))
      }
      return(list(p = NA_real_, label = "Fisher failed"))
    }

    # auto
    chi_obj <- suppressWarnings(tryCatch(stats::chisq.test(tab, correct = FALSE), error = function(e) NULL))
    if (is.null(chi_obj) || is.null(chi_obj$expected)) return(.cat_global_test(tab, method = "fisher"))
    if (any(chi_obj$expected < 5) || any(chi_obj$expected == 0)) return(.cat_global_test(tab, method = "fisher"))
    list(p = .clamp_p(chi_obj$p.value)[1], label = "Pearson chi-squared")
  }

  .robust_type2_p <- function(fit_lm, term) {
    V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) return(NA_real_)
    a2 <- tryCatch(car::Anova(fit_lm, type = 2, vcov. = V, test.statistic = "F"), error = function(e) NULL)
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

  .extract_anova_lr_p <- function(a) {
    if (is.null(a) || nrow(a) < 2) return(NA_real_)
    pcol <- intersect(c("Pr(>Chi)", "Pr(Chi)", "Pr(>Chisq)", "Pr(Chisq)"), colnames(a))
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
      return(list(p = NA_real_, label = "Insufficient data (multinomial)", LR = NA_real_, n = nrow(df_cc)))
    }

    f_full <- .fmla(var, c(compvar, covs))
    f_red  <- .fmla(var, covs)

    m_full <- tryCatch(nnet::multinom(f_full, data = df_cc, trace = FALSE), error = function(e) NULL)
    m_red  <- tryCatch(nnet::multinom(f_red,  data = df_cc, trace = FALSE), error = function(e) NULL)
    if (is.null(m_full) || is.null(m_red)) return(list(p = NA_real_, label = "Multinomial failed", LR = NA_real_, n = nrow(df_cc)))

    a <- tryCatch(stats::anova(m_red, m_full, test = "Chisq"), error = function(e) NULL)
    p <- .extract_anova_lr_p(a)
    if (is.na(p)) p <- .lrt_loglik_p(m_red, m_full)

    # best-effort LR stat if available
    LR <- NA_real_
    if (!is.null(a) && "LR stat." %in% colnames(a)) LR <- as.numeric(a[2, "LR stat."])

    list(p = .clamp_p(p)[1], label = "Multinomial LR", LR = LR, n = nrow(df_cc))
  }

  # LR-based V (used in your prior output for adjusted multinomial/logistic)
  .V_from_LR <- function(LR, n, n_outcome_levels, n_group_levels) {
    if (is.na(LR) || is.na(n) || n <= 0) return(NA_real_)
    k <- min(n_outcome_levels - 1, n_group_levels - 1)
    if (is.na(k) || k <= 0) return(NA_real_)
    sqrt(LR / (n * k))
  }

  # adjusted categorical pairwise: logistic if collapses to 2 levels, else multinomial
  .cat_pairwise_lr <- function(df_cc, var, compvar, covs, combos, PairwiseMethod) {

    out <- purrr::map_dfr(combos, function(cp) {

      sub <- df_cc[df_cc[[compvar]] %in% cp, , drop = FALSE]
      sub[[compvar]] <- droplevels(.as_factor(sub[[compvar]]))
      sub[[var]] <- droplevels(.as_factor(sub[[var]]))

      if (nlevels(sub[[compvar]]) < 2 || nrow(sub) < 5 || nlevels(sub[[var]]) < 2) {
        return(tibble::tibble(key = .pair_key(cp[1], cp[2]), p_val = NA_real_))
      }

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

  .strip_group_prefix <- function(x, compvar) {
    x <- .norm(x)
    x <- gsub(paste0("^", compvar, "\\s*=?\\s*"), "", x)
    x <- gsub(paste0("^", compvar), "", x)
    .norm(x)
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
  cols <- unique(c(CompVariable, cols))
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

  # ---- gtsummary display value control for dichotomous ----
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
  if (isTRUE(IncludeOverallStats)) {
    tbl <- gtsummary::tbl_summary(
      df,
      include   = tidyselect::all_of(Variables),
      missing   = if (IncludeMissing) "ifany" else "no",
      statistic = list(
        gtsummary::all_continuous()  ~ stat_cont,
        gtsummary::all_categorical() ~ stat_cat
      ),
      digits    = list(gtsummary::all_continuous() ~ ValueDigits),
      value     = value_list
    )
    if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
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
    value     = value_list
  )

  if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if (suppress_warnings) tbl <- suppressWarnings(tbl)
  if (nlevels(df[[CompVariable]]) < 2) return(tbl)

  full_group_levels <- levels(droplevels(df[[CompVariable]]))

  # ---- global p-values (+ adjusted p when covariates) and effect sizes ----
  pdat <- purrr::map_dfr(Variables, function(var) {

    is_cont <- isTRUE(treat_as_continuous[var])

    df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
    df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
    if (nlevels(df_vg[[CompVariable]]) < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = "Insufficient groups", Notes = NA_character_, effect_size = NA_character_))
    }

    p_un <- NA_real_
    p_adj <- NA_real_
    test_label <- NA_character_
    note <- NA_character_
    es <- NA_character_

    # --- continuous ---
    if (is_cont) {

      if (!is.null(Covariates)) {
        cols_cc <- c(var, CompVariable, Covariates)
        df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
        df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
        note <- if (isTRUE(IncludeNotes)) .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]])) else NA_character_

        if (nlevels(df_cc[[CompVariable]]) < 2 || nrow(df_cc) < 3) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Insufficient data (adjusted)", Notes = note, effect_size = es))
        }

        fit <- tryCatch(stats::lm(.fmla(var, c(CompVariable, Covariates)), data = df_cc), error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_,
                                test_label = "Model failed (adjusted)", Notes = note, effect_size = es))
        }

        # unadjusted ANOVA on the CC subset (so users can compare)
        p_un <- tryCatch(summary(stats::aov(.fmla(var, CompVariable), data = df_cc))[[1]][CompVariable, "Pr(>F)"],
                         error = function(e) NA_real_)
        p_un <- .clamp_p(p_un)[1]

        if (Parametric) {
          a2 <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
          p_adj <- if (!is.null(a2) && CompVariable %in% rownames(a2)) .clamp_p(a2[CompVariable, "Pr(>F)"])[1] else NA_real_
          test_label <- "ANCOVA (Type II)"

          if (isTRUE(AddEffectSize)) {
            eobj <- tryCatch(effectsize::eta_squared(a2, partial = TRUE), error = function(e) NULL)
            if (!is.null(eobj) && "Eta2_partial" %in% names(eobj)) {
              idx <- which(eobj$Parameter == CompVariable)[1]
              if (length(idx) && !is.na(idx)) {
                es <- paste0("\u03b7\u00b2p=", formatC(eobj$Eta2_partial[idx], digits = EffectSizeDigits, format = "f"))
              }
            }
          }

        } else {
          p_adj <- .robust_type2_p(fit, CompVariable)
          test_label <- "Robust ANCOVA (HC3 Type II)"

          if (isTRUE(AddEffectSize)) {
            V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
            a2r <- tryCatch(car::Anova(fit, type = 2, vcov. = V, test.statistic = "F"), error = function(e) NULL)
            eobj <- tryCatch(effectsize::eta_squared(a2r, partial = TRUE), error = function(e) NULL)
            if (!is.null(eobj) && "Eta2_partial" %in% names(eobj)) {
              idx <- which(eobj$Parameter == CompVariable)[1]
              if (length(idx) && !is.na(idx)) {
                es <- paste0("\u03b7\u00b2p=", formatC(eobj$Eta2_partial[idx], digits = EffectSizeDigits, format = "f"))
              }
            }
          }
        }

        return(tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj,
                              test_label = test_label, Notes = note, effect_size = es))
      }

      # unadjusted
      k <- nlevels(df_vg[[CompVariable]])

      if (Parametric) {
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

      if (isTRUE(AddEffectSize)) {
        aov_un <- tryCatch(stats::aov(.fmla(var, CompVariable), data = df_vg), error = function(e) NULL)
        eobj <- tryCatch(effectsize::eta_squared(aov_un, partial = TRUE), error = function(e) NULL)
        if (!is.null(eobj) && "Eta2_partial" %in% names(eobj)) {
          idx <- which(eobj$Parameter == CompVariable)[1]
          if (length(idx) && !is.na(idx)) {
            es <- paste0("\u03b7\u00b2p=", formatC(eobj$Eta2_partial[idx], digits = EffectSizeDigits, format = "f"))
          }
        }
      }

      return(tibble::tibble(variable = var, p_unadj = .clamp_p(p_un)[1], p_adj = NA_real_,
                            test_label = test_label, Notes = NA_character_, effect_size = es))
    }

    # --- categorical ---
    x <- .as_factor(df_vg[[var]])
    g <- .as_factor(df_vg[[CompVariable]])
    tab <- table(x, g)

    tst <- .cat_global_test(tab, method = CatMethod)
    p_un <- tst$p
    test_label <- tst$label
    p_adj <- NA_real_
    note <- NA_character_

    if (isTRUE(AddEffectSize)) {
      vobj <- tryCatch(effectsize::cramers_v(tab), error = function(e) NULL)
      vv <- tryCatch(as.numeric(vobj), error = function(e) NA_real_)
      if (!is.na(vv)) es <- paste0("V=", formatC(vv, digits = EffectSizeDigits, format = "f"))
    }

    if (!is.null(Covariates)) {
      cols_cc <- c(var, CompVariable, Covariates)
      df_cc <- df[stats::complete.cases(df[, cols_cc, drop = FALSE]), , drop = FALSE]
      df_cc[[CompVariable]] <- droplevels(.as_factor(df_cc[[CompVariable]]))
      df_cc[[var]] <- droplevels(.as_factor(df_cc[[var]]))
      note <- if (isTRUE(IncludeNotes)) .cc_levels_note(full_group_levels, levels(df_cc[[CompVariable]])) else NA_character_

      if (!is.null(Referent) && Referent %in% levels(df_cc[[CompVariable]])) {
        df_cc[[CompVariable]] <- stats::relevel(df_cc[[CompVariable]], ref = Referent)
      }

      # binary: logistic LR
      if (nlevels(df_cc[[var]]) == 2 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 5) {
        gm <- tryCatch(stats::glm(.fmla(var, c(CompVariable, Covariates)), data = df_cc, family = stats::binomial()),
                       error = function(e) NULL)
        if (!is.null(gm)) {
          d1 <- tryCatch(stats::drop1(gm, test = "Chisq"), error = function(e) NULL)
          if (!is.null(d1) && CompVariable %in% rownames(d1) && "Pr(>Chi)" %in% colnames(d1)) {
            p_adj <- .clamp_p(d1[CompVariable, "Pr(>Chi)"])[1]
            test_label <- "Logistic regression (LR)"

            # best-effort LR stat for V
            if (isTRUE(AddEffectSize) && "LRT" %in% colnames(d1)) {
              LR <- as.numeric(d1[CompVariable, "LRT"])
              vv <- .V_from_LR(LR, nrow(df_cc), nlevels(df_cc[[var]]), nlevels(df_cc[[CompVariable]]))
              if (!is.na(vv)) es <- paste0("V=", formatC(vv, digits = EffectSizeDigits, format = "f"))
            }
          }
        }
      }

      # multi-category: multinomial LR
      if (nlevels(df_cc[[var]]) >= 3 && nlevels(df_cc[[CompVariable]]) >= 2 && nrow(df_cc) >= 10) {
        g_lr <- .multinom_global_lr(df_cc, var, CompVariable, Covariates)
        if (!is.na(g_lr$p)) {
          p_adj <- .clamp_p(g_lr$p)[1]
          test_label <- g_lr$label

          if (isTRUE(AddEffectSize) && !is.na(g_lr$LR)) {
            vv <- .V_from_LR(g_lr$LR, g_lr$n, nlevels(df_cc[[var]]), nlevels(df_cc[[CompVariable]]))
            if (!is.na(vv)) es <- paste0("V=", formatC(vv, digits = EffectSizeDigits, format = "f"))
          }
        }
      }
    }

    tibble::tibble(variable = var, p_unadj = p_un, p_adj = p_adj, test_label = test_label, Notes = note, effect_size = es)
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
                                     effect_size = dplyr::if_else(.data$.is_main_row, .data$effect_size, NA_character_),
                                     Notes       = dplyr::if_else(.data$.is_main_row, .data$Notes, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

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

  # ---- pairwise ----
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

      # continuous pairwise: covariate-adjusted via emmeans
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

        # "none" means no adjustment in emmeans
        adj_arg <- if (identical(PairwiseMethod, "none")) "none" else PairwiseMethod
        res <- tryCatch(as.data.frame(summary(ctr, adjust = adj_arg)), error = function(e) NULL)
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
          p_val = as.numeric(res$p.value)
        )

        return(out_template %>%
                 dplyr::mutate(p_val = got$p_val[match(.data$key, got$key)]) %>%
                 dplyr::select("variable", "contrast_label", "p_val"))
      }

      # continuous pairwise: unadjusted (no covariates)
      if (is_cont && is.null(Covariates)) {
        df_vg <- df[stats::complete.cases(df[, c(var, CompVariable), drop = FALSE]), , drop = FALSE]
        df_vg[[CompVariable]] <- droplevels(.as_factor(df_vg[[CompVariable]]))
        if (nlevels(df_vg[[CompVariable]]) < 2) {
          return(dplyr::select(out_template, "variable", "contrast_label", "p_val"))
        }

        x <- df_vg[[var]]
        g <- df_vg[[CompVariable]]

        method_arg <- if (identical(PairwiseMethod, "none")) "none" else PairwiseMethod

        pw_mat <- tryCatch(
          if (Parametric) stats::pairwise.t.test(x, g, pool.sd = FALSE, p.adjust.method = method_arg)$p.value
          else stats::pairwise.wilcox.test(x, g, p.adjust.method = method_arg)$p.value,
          error = function(e) NULL
        )

        got <- out_template %>%
          dplyr::mutate(
            p_val = purrr::map_dbl(strsplit(.data$contrast_label, " - "), function(cp) {
              a <- cp[1]; b <- cp[2]
              if (is.null(pw_mat)) return(NA_real_)
              if (!is.null(rownames(pw_mat)) && !is.null(colnames(pw_mat))) {
                if (a %in% rownames(pw_mat) && b %in% colnames(pw_mat)) return(as.numeric(pw_mat[a, b]))
                if (b %in% rownames(pw_mat) && a %in% colnames(pw_mat)) return(as.numeric(pw_mat[b, a]))
              }
              NA_real_
            })
          )

        return(got %>% dplyr::select("variable", "contrast_label", "p_val"))
      }

      # categorical pairwise: adjusted only (per requirements)
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
                                           ~ dplyr::if_else(.data$.is_main_row, .clamp_p(.), NA_real_)
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

  # ---- ensure Notes is last if included ----
  tbl <- tbl %>%
    gtsummary::modify_table_body(~{
      tb <- .x
      if ("Notes" %in% names(tb)) tb <- dplyr::relocate(tb, .data$Notes, .after = dplyr::last_col())
      tb
    })

  # ---- caption ----
  cont_disp <- if (ParametricDisplay) "mean (SD)" else "median [IQR]"
  mode_txt <- if (Parametric) "parametric tests" else "robust or nonparametric tests"
  cov_txt <- if (is.null(Covariates)) "unadjusted" else paste0("adjusted for ", paste(Covariates, collapse = ", "))
  pw_txt <- if (AddPairwise) paste0("pairwise: ", PairwiseMethod, if (!is.null(Referent)) paste0(" (ref=", Referent, ")") else "") else "pairwise: none"

  cap <- paste0(
    "Continuous summaries shown as ", cont_disp, "; ",
    mode_txt, "; ",
    cov_txt, "; ",
    "categorical: ", CatMethod, "; ",
    "multi-category adjusted: ", MultiCatAdjusted, "; ",
    pw_txt, "."
  )

  tbl <- tbl %>% gtsummary::modify_caption(cap)

  tbl
}
