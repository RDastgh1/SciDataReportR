#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, and Pairwise Contrasts
#'
#' @description
#' Creates a comparison table using {gtsummary} summarizing continuous and categorical variables
#' by a grouping factor, with optional covariate adjustment for continuous outcomes only (ANCOVA),
#' effect‐size calculations, and pairwise contrasts with p‐value adjustments (or logistic regression for binary outcomes),
#' with optional inclusion of overall N and missing (“Unknown”) rows.  Variables with fewer than
#' two unique non‐NA values are dropped automatically.  Categorical tests that would error (e.g. perfect separation)
#' are caught and yield NA p‐values.
#'
#' @param DataFrame A `data.frame` containing the raw data.
#' @param Variables Character vector of column names to include in the comparison.
#' @param CompVariable Character string specifying the grouping (comparison) variable.
#' @param Covariates Optional character vector of covariate column names for adjustment (default `NULL`).
#' @param ValueDigits Integer; number of digits for summary statistics (default `2`).
#' @param pDigits Integer; number of digits for formatted p‐values (default `3`).
#' @param AddEffectSize Logical; include an effect‐size column? (default `FALSE`).
#' @param EffectSizeDigits Integer; number of digits for effect sizes (default `2`).
#' @param AddPairwise Logical; include pairwise contrast columns? (default `FALSE`).
#' @param PairwiseMethod Character; p‐value adjustment method for contrasts (default `"bonferroni"`; `"none"` for unadjusted).
#' @param Parametric Logical; use parametric tests if `TRUE` (ANOVA/ANCOVA), otherwise nonparametric (Kruskal–Wallis) (default `TRUE`).
#' @param IncludeOverallN Logical; include a column with the overall N in the table? (default `FALSE`).
#' @param IncludeMissing Logical; include a row summarizing missing data (“Unknown”)? (default `FALSE`).
#' @param suppress_warnings Logical; suppress intermediate warnings from gtsummary (default `FALSE`).
#'
#' @return A `gtsummary::tbl_summary` object.  Variables or group‐levels with fewer than two observations yield NA p‐values.
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
    suppress_warnings = FALSE
) {

  req_pkgs <- c("gtsummary","dplyr","car","emmeans","broom","effectsize","purrr","tidyr")
  available <- vapply(req_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
  missing_pkgs <- req_pkgs[!available]
  if(length(missing_pkgs)) {
    stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))
  }

  # Validate inputs
  if(!CompVariable %in% names(DataFrame)) {
    stop("Grouping variable not found: ", CompVariable)
  }
  if(!all(Variables %in% names(DataFrame))) {
    stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  }
  if(!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) {
    stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))
  }

  # Prepare data
  df <- DataFrame %>%
    dplyr::select(all_of(c(CompVariable, Variables, Covariates)))
  df[[CompVariable]] <- factor(df[[CompVariable]])

  # Drop variables with <2 unique non-NA values
  keep_vars <- Variables[sapply(df[Variables], function(x) length(unique(x[!is.na(x)]))) >= 2]
  dropped   <- setdiff(Variables, keep_vars)
  if(length(dropped)) {
    warning("Dropping variable(s) with <2 unique values: ", paste(dropped, collapse = ", "))
  }
  Variables <- keep_vars
  if(length(Variables) == 0) {
    stop("No variables left after dropping those with <2 unique values.")
  }

  # Build base tbl_summary ----
  tbl <- gtsummary::tbl_summary(
    df,
    by      = CompVariable,
    missing = if(IncludeMissing) "ifany" else "no",
    statistic = list(
      all_continuous()  ~ if(Parametric) "{mean} ({sd})" else "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(all_continuous() ~ ValueDigits)
  )
  if(IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if(suppress_warnings) tbl <- suppressWarnings(tbl)

  # 4a) If grouping factor has <2 levels, return summary only
  if(nlevels(df[[CompVariable]]) < 2) {
    warning("Grouping variable '", CompVariable, "' has fewer than 2 levels; skipping comparisons.")
    return(tbl)
  }

  # Compute p-values ----
  pvals <- purrr::map_dfr(Variables, function(var) {
    # Subset non-missing for this var
    df_var <- df[!is.na(df[[var]]), , drop = FALSE]
    df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])

    # If <2 groups in subset, return NA
    if(nlevels(df_var[[CompVariable]]) < 2) {
      return(tibble::tibble(variable = var, p_unadj = NA_real_, p_adj = NA_real_))
    }

    # Unadjusted p
    if(is.numeric(df_var[[var]])) {
      pu <- if(Parametric) {
        summary(aov(as.formula(paste(var, "~", CompVariable)), data = df_var))[[1]][CompVariable, "Pr(>F)"]
      } else {
        kruskal.test(as.formula(paste(var, "~", CompVariable)), data = df_var)$p.value
      }
    } else {
      tbl0 <- table(df_var[[var]], df_var[[CompVariable]])
      # drop empty margins
      tbl0 <- tbl0[rowSums(tbl0) > 0, colSums(tbl0) > 0, drop = FALSE]
      pu <- if(nrow(tbl0) < 2 || ncol(tbl0) < 2) {
        NA_real_
      } else {
        tryCatch({
          if(any(suppressWarnings(chisq.test(tbl0)$expected) < 5)) {
            fisher.test(tbl0, simulate.p.value = TRUE, B = 1e5)$p.value
          } else {
            chisq.test(tbl0)$p.value
          }
        }, error = function(e) NA_real_)
      }
    }

    # Adjusted p
    pad <- NA_real_
    if(Parametric && !is.null(Covariates)) {
      if(is.numeric(df_var[[var]])) {
        m1  <- lm(as.formula(paste(var, "~", CompVariable, "+", paste(Covariates, collapse = "+"))),
                  data = df_var)
        pad <- car::Anova(m1, type = 2)[CompVariable, "Pr(>F)"]
      } else if(nlevels(df_var[[var]]) == 2) {
        gm  <- glm(as.formula(paste(var, "~", CompVariable, "+", paste(Covariates, collapse = "+"))),
                   data = df_var, family = binomial)
        pad <- drop1(gm, test = "Chisq")[CompVariable, "Pr(>Chi)"]
      }
    }

    tibble::tibble(variable = var, p_unadj = pu, p_adj = pad)
  })

  test_note <- if(Parametric) {
    "Continuous: ANOVA/ANCOVA; Categorical: chi-square/Fisher (or logistic adj)"
  } else {
    "Continuous: Kruskal–Wallis; Categorical: chi-square/Fisher"
  }

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   left_join(pvals, by = "variable") %>%
                                   mutate(
                                     p.value     = coalesce(p_adj, p_unadj),
                                     p.value_fmt = gtsummary::style_pvalue(p.value, digits = pDigits),
                                     p.value_fmt = ifelse(row_type == "label", p.value_fmt, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_table_styling(
      columns     = "p.value_fmt",
      rows        = p.value <= 0.05 & row_type == "label",
      text_format = "bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_footnote(p.value_fmt ~ test_note)


  tbl %>% gtsummary::modify_caption(
    sprintf(
      "Comparison Table (%s analysis%s)",
      if(Parametric) "parametric" else "non-parametric",
      if(!is.null(Covariates) && Parametric)
        paste0("; adjusted for ", paste(Covariates, collapse = ","))
      else ""
    )
  )
}
