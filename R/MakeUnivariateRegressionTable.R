#' Univariate Regression Table
#'
#' Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified).
#'
#' @details
#' Fits one model per outcome-predictor pair and collects every result into a
#' single table. "Univariate" means each predictor is tested on its own: for
#' three outcomes and eight predictors this is 24 separate models, not one
#' model with eight terms. That is the screening step - finding which
#' associations exist at all - and it is deliberately different from
#' [MultivariableRegressionTable()], which puts all the predictors in one model
#' and reports each one adjusted for the others.
#'
#' The model family is chosen per outcome. With `Method = "auto"` (the
#' default), numeric outcomes get linear regression and two-level outcomes get
#' logistic regression, so a mixed set of outcomes can be screened in one call
#' and each is still modeled correctly. `"lm"` or `"logistic"` force one family
#' for everything.
#'
#' @section Three views of the same results:
#' The return value carries the same estimates in three shapes, for three
#' different jobs:
#'
#' * `FormattedTable` - a wide `gt` table with one spanner per outcome,
#'   confidence intervals, stars, and significant results emphasized.
#' * `LargeTable` - the same wide layout with the individual model statistics,
#'   when the raw numbers matter more than the presentation.
#' * `Results` - one tidy row per estimated term, with `Estimate`, `StdError`,
#'   confidence bounds, `PValue`, and the labels. This is the one to work with
#'   programmatically: filter it, correct it with [ApplyFDRCorrection()], or
#'   pass it straight to [PlotForestFromTable()].
#'
#' Categorical predictors contribute one row per non-reference level, with the
#' level in `Level` and the level being compared against in `ReferenceValue`.
#'
#' @section Options that change the estimates:
#' `Standardize = TRUE` z-scores the numeric variables first, so coefficients
#' come back in standard deviations per standard deviation. Use it when
#' predictors are on scales that are not comparable to each other - it is what
#' makes "which predictor matters more" a meaningful question.
#'
#' `covariates` adds the named variables to every model, turning a screen of
#' raw associations into a screen of associations adjusted for them. This is
#' usually the difference between "these biomarkers track with the outcome" and
#' "these biomarkers track with the outcome beyond age and sex".
#'
#' `LogisticExponentiate = TRUE` (the default) reports logistic results as odds
#' ratios rather than log-odds, so the null value in the table is 1 rather
#' than 0. `ReturnModels = TRUE` keeps the fitted model objects for
#' residual checks; it is off by default because a large screen would otherwise
#' hold hundreds of models in memory.
#'
#' Running many models at once means many p-values. The default fdr_scope =
#' "none" leaves PValueFDR unavailable. Select a correction scope to compute
#' BH-adjusted p-values within each Group.
#'
#' @seealso [PlotForestFromTable()] to visualize `Results`,
#'   [MultivariableRegressionTable()] for mutually adjusted models, and
#'   [ApplyFDRCorrection()] for multiple-comparison correction.
#'
#' @param data Dataframe containing the variables
#' @param outcome_vars Character vector of outcome variable names
#' @param predictor_vars Character vector of predictor variable names
#' @param covariates Character vector of covariate variable names (default: NULL)
#' @param Standardize Logical indicating whether to standardize numeric variables (default: FALSE)
#' @param Method Character. Regression method to use. `"auto"` detects linear
#'   regression for numeric outcomes and logistic regression for two-level
#'   outcomes. `"lm"` and `"logistic"` force one model family for all outcomes.
#' @param LogisticExponentiate Logical. If `TRUE`, logistic regression estimates
#'   are exponentiated and reported as odds ratios.
#' @param ReturnModels Logical. If `TRUE`, return fitted model objects in
#'   `ModelSummaries`. Default is `FALSE` to keep large screening runs lighter.
#' @param Relabel Logical; if TRUE (default), display attached variable labels.
#' @param TreatOrdinalAs How ordinal outcomes and predictors are handled.
#' @importFrom sjlabelled get_label set_label
#' @return A list containing:
#'   - FormattedTable: A `gt` table with formatted regression results
#'   - LargeTable: A `gt` table with unformatted regression results
#'   - Results: A tidy dataframe with one row per estimated term. Columns:
#'     `Outcome`, `OutcomeLabel`, `OutcomeFamily`, `EffectType`, `Predictor`,
#'     `PredictorLabel`, `Term`, `Level`, `TermLabel`, `N`, `Estimate`,
#'     `StdError`, `ConfLow`, `ConfHigh`, `PValue`, `Significant`, and
#'     `ReferenceValue`. This dataframe can be filtered and passed directly to
#'     [PlotForestFromTable()].
#'   - ModelSummaries: A list of fitted model objects when `ReturnModels = TRUE`,
#'     otherwise `NULL`
#'   - Metadata: Outcome families and analysis settings
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param OutcomeVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @param PredictorVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param Covars \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' # Make the reference example deterministic: age is visibly associated with
#' # AXL, so the displayed table includes a clear significant result.
#' ExampleData <- Labelled
#' set.seed(20260810)
#' ExampleData$AXL <- as.numeric(scale(ExampleData$age)) * 2 +
#'   stats::rnorm(nrow(ExampleData), sd = 0.25)
#' attr(ExampleData$AXL, "label") <- sjlabelled::get_label(Labelled$AXL)
#'
#' vars_Outcomes <- c("AXL", "tau", "p_tau")
#' vars_Predictors <- c("age", "sex", "Adiponectin", "Cortisol")
#'
#' # Three outcomes by four predictors: twelve separate models, one table
#' urt <- MakeUnivariateRegressionTable(
#'   data = ExampleData,
#'   outcome_vars = vars_Outcomes,
#'   predictor_vars = vars_Predictors
#' )
#'
#' # Format 1: the report-ready table
#' urt$FormattedTable
#'
#' # Format 2: the same results unformatted
#' urt$LargeTable
#'
#' # Format 3: the tidy dataframe, one row per estimated term
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(
#'     dplyr::mutate(
#'       dplyr::select(
#'         urt$Results, Outcome, Predictor, TermLabel, Level, N,
#'         Estimate, StdError, ConfLow, ConfHigh, PValue, Significant
#'       ),
#'       dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 3))
#'     ),
#'     height = "320px", full_width = TRUE
#'   )
#' )))
#'
#' # Filter, correct, and plot
#' urt$Results$PValueFDR <- ApplyFDRCorrection(urt$Results$PValue)
#' PlotForestFromTable(urt$Results[urt$Results$PValueFDR < 0.1, ])
#'
#' # Standardized coefficients, in SD units
#' MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = vars_Outcomes,
#'   predictor_vars = vars_Predictors,
#'   Standardize = TRUE
#' )$FormattedTable
#'
#' # Every model adjusted for age and sex
#' MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = vars_Outcomes,
#'   predictor_vars = c("Adiponectin", "Cortisol", "Ferritin"),
#'   covariates = c("age", "sex")
#' )$FormattedTable
#'
#' # A binary outcome: logistic regression, reported as odds ratios
#' urt_Logistic <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = "Diagnosis",
#'   predictor_vars = c("age", "sex", "AXL", "tau", "p_tau")
#' )
#' urt_Logistic$Metadata
#' urt_Logistic$FormattedTable
#'
#' # Log-odds instead of odds ratios
#' MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = "Diagnosis",
#'   predictor_vars = c("age", "AXL", "tau"),
#'   LogisticExponentiate = FALSE
#' )$FormattedTable
#'
#' # Keep the fitted models for residual checks
#' urt_Models <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = "AXL",
#'   predictor_vars = c("age", "Cortisol"),
#'   ReturnModels = TRUE
#' )
#' summary(urt_Models$ModelSummaries[[1]])
#' }
#'
#' @param group_var Optional categorical stratification variable. Fits each
#'   observed level separately, excluding missing group values. Group is NA
#'   without stratification. The name "Pooled" is reserved.
#' @param include_pooled Include a pooled model adjusted for group_var (TRUE).
#' @param fdr_scope BH correction family: "none" (default, no adjustment),
#'   "per_predictor" (across outcomes), "per_outcome" (across predictors), or
#'   "matrix" (all tests). Each Group is corrected separately.
#' @section Stratified scaling:
#' With group_var and Standardize = TRUE, numeric model variables use the same
#' pooled complete-case means and SDs in every group. Logistic outcomes are
#' not scaled. Metadata records pooled scaling, model formulas, and group N.
#' @param IncludeTrendSymbol Logical; add a plus sign for `0.05 < p <= 0.10`.
#' @export
#'
#' @inheritParams inference-results
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
MakeUnivariateRegressionTable <- function(data,
    outcome_vars,
    predictor_vars,
    covariates = NULL,
    Standardize = FALSE,
    Method = c("auto", "lm", "logistic"),
    LogisticExponentiate = TRUE,
    ReturnModels = FALSE,
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    IncludeTrendSymbol = FALSE,
    Data = lifecycle::deprecated(),
    OutcomeVars = lifecycle::deprecated(),
    PredictorVars = lifecycle::deprecated(),
    Covars = lifecycle::deprecated(),
    Pthresh = 0.05,
    PVar = "PValue",
    group_var = NULL,
    include_pooled = TRUE,
    fdr_scope = c("none", "per_predictor", "per_outcome", "matrix")) {
  .scidr_inference <- ScidrStartInference("MakeUnivariateRegressionTable", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "MakeUnivariateRegressionTable(Data)", "MakeUnivariateRegressionTable(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(OutcomeVars)) {
    lifecycle::deprecate_warn("19.15.0", "MakeUnivariateRegressionTable(OutcomeVars)", "MakeUnivariateRegressionTable(outcome_vars)")
    outcome_vars <- OutcomeVars
  }
  if (!missing(outcome_vars)) OutcomeVars <- outcome_vars
  if (lifecycle::is_present(PredictorVars)) {
    lifecycle::deprecate_warn("19.15.0", "MakeUnivariateRegressionTable(PredictorVars)", "MakeUnivariateRegressionTable(predictor_vars)")
    predictor_vars <- PredictorVars
  }
  if (!missing(predictor_vars)) PredictorVars <- predictor_vars
  if (lifecycle::is_present(Covars)) {
    lifecycle::deprecate_warn("19.15.0", "MakeUnivariateRegressionTable(Covars)", "MakeUnivariateRegressionTable(covariates)")
    covariates <- Covars
  }
  Covars <- covariates

  Method <- match.arg(Method)
  fdr_scope <- match.arg(fdr_scope)
  if (!is.logical(include_pooled) || length(include_pooled) != 1L || is.na(include_pooled)) stop("include_pooled must be TRUE or FALSE.")
  group_levels <- NA_character_
  if (!is.null(group_var)) {
    ScidrValidateVariable(Data, group_var, "group_var")
    if (is.numeric(Data[[group_var]])) stop("group_var must be categorical.")
    if (group_var %in% c(OutcomeVars, PredictorVars)) stop("group_var cannot be a tested outcome or predictor.")
    Data[[group_var]] <- factor(Data[[group_var]])
    group_levels <- levels(droplevels(Data[[group_var]]))
    if (!length(group_levels) || "Pooled" %in% group_levels) stop("group_var needs observed levels other than the reserved name 'Pooled'.")
    Covars <- setdiff(Covars, group_var)
  }
  model_groups <- if (is.null(group_var)) NA_character_ else c(group_levels, if (include_pooled) "Pooled")
  scaling_records <- list()
  group_counts <- list()
  model_records <- list()

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }
  if (!is.character(OutcomeVars) || length(OutcomeVars) == 0) {
    stop("OutcomeVars must be a non-empty character vector.")
  }
  if (!is.character(PredictorVars) || length(PredictorVars) == 0) {
    stop("PredictorVars must be a non-empty character vector.")
  }
  if (!is.null(Covars) && !is.character(Covars)) {
    stop("Covars must be NULL or a character vector.")
  }
  if (!is.logical(Standardize) || length(Standardize) != 1) {
    stop("Standardize must be TRUE or FALSE.")
  }
  if (!is.logical(LogisticExponentiate) || length(LogisticExponentiate) != 1) {
    stop("LogisticExponentiate must be TRUE or FALSE.")
  }
  if (!is.logical(ReturnModels) || length(ReturnModels) != 1) {
    stop("ReturnModels must be TRUE or FALSE.")
  }
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop(
      "Package 'gt' is required by MakeUnivariateRegressionTable(). ",
      "Install it with install.packages('gt')."
    )
  }

  all_model_vars <- unique(c(OutcomeVars, PredictorVars, Covars))
  missing_vars <- setdiff(all_model_vars, names(Data))
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in Data: ", paste(missing_vars, collapse = ", "))
  }
  TreatOrdinalAs <- match.arg(TreatOrdinalAs, c("Categorical", "Continuous", "Both", "Exclude"))
  if (TreatOrdinalAs == "Both") {
    stop("TreatOrdinalAs = 'Both' is not meaningful for MakeUnivariateRegressionTable().", call. = FALSE)
  }

  ordinal_reference <- ConvertOrdinalToNumeric(
    Data, all_model_vars, TreatOrdinalAs = "Categorical", ReturnMetadata = TRUE
  )
  if (identical(TreatOrdinalAs, "Exclude") && length(ordinal_reference$ordinal_variables)) {
    stop("TreatOrdinalAs = 'Exclude' cannot be used when ordinal model variables are explicitly supplied.", call. = FALSE)
  }
  prep <- ConvertOrdinalToNumeric(
    Data, all_model_vars, TreatOrdinalAs = TreatOrdinalAs, ReturnMetadata = TRUE
  )
  Data <- prep$data
  display_labels <- ScidrDisplayLabels(Data, all_model_vars, Relabel)

  OutcomeVars <- ScidrDropUnavailableOutcomes(Data, OutcomeVars, PredictorVars, Covars)

  Wide_mod_list <- list()
  results_list <- list()
  outcome_metadata <- list()

  for (yVarIndex in seq_along(OutcomeVars)) {
    YVar <- OutcomeVars[yVarIndex]
    outcome_family <- ScidrUnivariateOutcomeFamily(Data[[YVar]], YVar, Method)
    outcome_info <- ScidrUnivariatePrepareOutcome(Data[[YVar]], YVar, outcome_family)
    model_data_full <- Data
    model_data_full[[YVar]] <- outcome_info$Value

    mod_list <- list()
    for (xVarIndex in seq_along(PredictorVars)) {
      xVar <- PredictorVars[xVarIndex]
      for (current_group in model_groups) {
      base_covars <- Covars
      Covars <- if (!is.null(group_var) && current_group == "Pooled") unique(c(base_covars, group_var)) else base_covars
      tryCatch(expr = {
        f <- ScidrRegressionFormula(c(xVar, Covars), outcome = YVar)
        ModelData <- model_data_full
        if (!is.null(group_var)) {
          vars_Model <- unique(c(YVar, xVar, base_covars, group_var))
          ModelData <- ModelData[stats::complete.cases(ModelData[, vars_Model, drop = FALSE]), , drop = FALSE]
        }
        if (Standardize) {
          if (is.null(group_var)) ModelData <- ModelData %>% dplyr::select(dplyr::all_of(unique(c(YVar, xVar, Covars))))
          numeric_cols <- vapply(ModelData, is.numeric, logical(1)) & names(ModelData) %in% c(YVar, xVar, Covars)
          if (outcome_family == "logistic") numeric_cols[YVar] <- FALSE
          if (any(numeric_cols)) {
            if (!is.null(group_var)) {
              scaling_records[[paste(YVar, xVar, sep = "\r")]] <- tibble::tibble(
                Outcome = YVar, Predictor = xVar, Variable = names(ModelData)[numeric_cols],
                SD = unname(vapply(ModelData[, numeric_cols, drop = FALSE], stats::sd, numeric(1))),
                Center = unname(vapply(ModelData[, numeric_cols, drop = FALSE], mean, numeric(1))), N = nrow(ModelData))
            }
            ModelData[, numeric_cols] <- ScidrObserveInference(scale(ModelData[, numeric_cols, drop = FALSE]))
          }
        }
        if (!is.null(group_var) && current_group != "Pooled") {
          ModelData <- ModelData[ModelData[[group_var]] == current_group, , drop = FALSE]
        }
        group_counts[[length(group_counts) + 1L]] <- tibble::tibble(Outcome = YVar, Predictor = xVar,
          Group = current_group, N = nrow(ModelData[stats::complete.cases(ModelData[, unique(c(YVar, xVar, Covars)), drop = FALSE]), ]))
        mod <- if (outcome_family == "logistic") {
          ScidrObserveInference(stats::glm(formula = f, data = ModelData, family = stats::binomial()))
        } else {
          ScidrObserveInference(stats::lm(formula = f, data = ModelData))
        }
        model_records[[length(model_records) + 1L]] <- tibble::tibble(
          Outcome = YVar, Predictor = xVar, Group = current_group,
          Formula = paste(deparse(stats::formula(mod)), collapse = " "), N = stats::nobs(mod))
        if (ReturnModels) {
          mod_list[[if (is.null(group_var)) xVar else paste(xVar, current_group, sep = ":")]] <- mod
        }

        model_summary <- summary(mod)
        coefficient_table <- as.data.frame(model_summary$coefficients)
        coefficient_table$Term <- rownames(coefficient_table)
        rownames(coefficient_table) <- NULL

        if (outcome_family == "logistic") {
          estimate_col <- "Estimate"
          std_error_col <- "Std. Error"
          p_col <- "Pr(>|z|)"
          critical_value <- stats::qnorm(0.975)
        } else {
          estimate_col <- "Estimate"
          std_error_col <- "Std. Error"
          p_col <- "Pr(>|t|)"
          critical_value <- stats::qt(0.975, df = stats::df.residual(mod))
        }

        labels <- display_labels[c(PredictorVars, OutcomeVars)]
        predictor_label <- labels[[xVar]]
        effect_type <- ifelse(
          outcome_family == "logistic" && LogisticExponentiate,
          "Odds ratio",
          "Estimate"
        )

        model_terms <- attr(stats::terms(mod), "term.labels")
        covariate_terms <- ScidrQuoteFormulaNames(Covars)
        predictor_term <- ScidrQuoteFormulaNames(xVar)
        keep_terms <- model_terms[model_terms %!in% covariate_terms]
        keep_rows <- coefficient_table$Term %in% keep_terms
        if (length(keep_terms) > 0) {
          keep_rows <- keep_rows | startsWith(coefficient_table$Term, paste0(keep_terms, ""))
        }
        coefficient_table <- coefficient_table[keep_rows, , drop = FALSE]

        if (nrow(coefficient_table) == 0) stop("No estimable predictor coefficients.")

        raw_estimate <- coefficient_table[[estimate_col]]
        std_error <- coefficient_table[[std_error_col]]
        conf_low <- raw_estimate - critical_value * std_error
        conf_high <- raw_estimate + critical_value * std_error
        estimate <- raw_estimate

        if (outcome_family == "logistic" && LogisticExponentiate) {
          estimate <- exp(raw_estimate)
          conf_low <- exp(conf_low)
          conf_high <- exp(conf_high)
        }

        term_label <- coefficient_table$Term
        level <- rep(NA_character_, nrow(coefficient_table))
        if (is.factor(ModelData[[xVar]]) || is.character(ModelData[[xVar]])) {
          level <- substring(
            coefficient_table$Term,
            first = nchar(predictor_term) + 1L
          )
          level[level == coefficient_table$Term | level == ""] <- NA_character_
          term_label <- ifelse(
            !is.na(level),
            paste0(predictor_label, " : ", level),
            predictor_label
          )
        } else {
          term_label <- rep(predictor_label, nrow(coefficient_table))
        }

        results_list[[length(results_list) + 1]] <- data.frame(
          Group = current_group,
          Outcome = YVar,
          OutcomeLabel = display_labels[[YVar]],
          OutcomeFamily = outcome_family,
          EffectType = effect_type,
          Predictor = xVar,
          PredictorLabel = predictor_label,
          Term = coefficient_table$Term,
          Level = level,
          TermLabel = term_label,
          N = stats::nobs(mod),
          Estimate = estimate,
          StdError = std_error,
          ConfLow = conf_low,
          ConfHigh = conf_high,
          PValue = coefficient_table[[p_col]],
          Significant = coefficient_table[[p_col]] < 0.05,
          ReferenceValue = ifelse(effect_type == "Odds ratio", 1, 0),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }, error = function(e) {
        warning(paste("Error processing", YVar, "and", xVar,
                   ": ", conditionMessage(e)), call. = FALSE)
      })
      Covars <- base_covars
      }
    }
    if (ReturnModels) {
      Wide_mod_list[[YVar]] <- mod_list
    }
    outcome_label <- display_labels[[YVar]]
    outcome_metadata[[YVar]] <- data.frame(
      Outcome = YVar,
      OutcomeLabel = outcome_label,
      OutcomeFamily = outcome_family,
      ReferenceLevel = outcome_info$ReferenceLevel,
      EventLevel = outcome_info$EventLevel,
      stringsAsFactors = FALSE
    )
  }

  results <- dplyr::bind_rows(results_list)
  if (!nrow(results)) results <- ScidrEmptyInferenceResults()
  results$PValueFDR <- NA_real_
  if (fdr_scope != "none") {
    families <- if (fdr_scope == "per_predictor") results$Predictor else if (fdr_scope == "per_outcome") results$Outcome else rep("matrix", nrow(results))
    families <- paste(results$Group, families, sep = "\r")
    for (indices in split(seq_len(nrow(results)), families)) results$PValueFDR[indices] <- stats::p.adjust(results$PValue[indices], "BH")
  }
  if (is.null(group_var)) results <- results[, c(setdiff(names(results), "Group"), "Group")]
  results$Significant <- (if (PVar == "PValue") results$PValue else results$PValueFDR) < Pthresh
  FinalTable <- ScidrUnivariateGtTable(results, formatted = FALSE)
  FinalFormattedTable <- ScidrUnivariateGtTable(
    results, formatted = TRUE, IncludeTrendSymbol = IncludeTrendSymbol
  )

  metadata <- list(
    GroupLevels = if (is.null(group_var)) NULL else group_levels,
    GroupN = dplyr::bind_rows(group_counts),
    GroupModels = dplyr::bind_rows(model_records),
    PooledScaling = dplyr::bind_rows(scaling_records),
    ScalingScope = if (is.null(group_var)) "Existing ungrouped scaling" else "Pooled complete model sample, shared across groups",
    Outcomes = dplyr::bind_rows(outcome_metadata),
    AnalysisSettings = list(
      Method = Method,
      OutcomeVars = OutcomeVars,
      PredictorVars = PredictorVars,
      Covars = Covars,
      Standardize = Standardize,
      LogisticExponentiate = LogisticExponentiate,
      ReturnModels = ReturnModels
    )
  )

  return(list(FormattedTable = FinalFormattedTable, LargeTable = FinalTable,
              Results = results,
              ModelSummaries = if (ReturnModels) Wide_mod_list else NULL,
              Metadata = metadata))
}

#' @description `UnivariateRegressionTable()` was renamed to
#'   `MakeUnivariateRegressionTable()` in SciDataReportR 20.5.0 to match the
#'   package's `Make*` naming convention. It remains available as a
#'   backwards-compatible synonym.
#' @rdname MakeUnivariateRegressionTable
#' @export
UnivariateRegressionTable <- function(data,
    outcome_vars,
    predictor_vars,
    covariates = NULL,
    Standardize = FALSE,
    Method = c("auto", "lm", "logistic"),
    LogisticExponentiate = TRUE,
    ReturnModels = FALSE,
    Data = lifecycle::deprecated(),
    OutcomeVars = lifecycle::deprecated(),
    PredictorVars = lifecycle::deprecated(),
    Covars = lifecycle::deprecated()) {
  lifecycle::deprecate_soft("20.5.0", "UnivariateRegressionTable()", "MakeUnivariateRegressionTable()")
  call <- match.call()
  call[[1L]] <- MakeUnivariateRegressionTable
  eval.parent(call)
}

ScidrUnivariateGtTable <- function(results, formatted = TRUE, IncludeTrendSymbol = FALSE) {
  # Results remains the source of truth. These two tables deliberately only
  # reshape it for reporting, so users can always filter Results for plotting.
  results <- results %>%
    dplyr::mutate(
      RowKey = paste(.data$Predictor, .data$Term, if ("Group" %in% names(results)) results$Group else "", sep = "\r")
    )

  if ("Group" %in% names(results) && any(!is.na(results$Group))) {
    results$TermLabel <- paste0(results$TermLabel, " (", results$Group, ")")
  }
  row_data <- results %>%
    dplyr::distinct(.data$RowKey, .data$TermLabel) %>%
    dplyr::rename(Variable = "TermLabel")
  outcome_order <- unique(results$Outcome)

  table_data <- row_data %>% dplyr::select(-tidyselect::all_of("RowKey"))
  column_labels <- list(Variable = "Variable")
  spanners <- list()
  numeric_columns <- character(0)
  significant_rows <- list()

  for (outcome_index in seq_along(outcome_order)) {
    outcome <- outcome_order[[outcome_index]]
    outcome_results <- results %>%
      dplyr::filter(.data$Outcome == outcome) %>%
      dplyr::select(-tidyselect::all_of("Outcome"))
    outcome_results <- outcome_results[match(row_data$RowKey, outcome_results$RowKey), , drop = FALSE]

    column_prefix <- paste0("Outcome", outcome_index)
    outcome_label <- outcome_results$OutcomeLabel[[which(!is.na(outcome_results$OutcomeLabel))[1]]]
    effect_type <- outcome_results$EffectType[[which(!is.na(outcome_results$EffectType))[1]]]
    effect_label <- if (identical(effect_type, "Odds ratio")) {
      "Odds ratio (95% CI)"
    } else {
      "Estimate (95% CI)"
    }

    if (formatted) {
      effect_column <- paste0(column_prefix, "_EffectCI")
      p_column <- paste0(column_prefix, "_P")
      stars <- ifelse(
        outcome_results$Significant,
        ScidrPValueStars(
          ScidrInferenceP(outcome_results$PValue, outcome_results$PValueFDR),
          ns_label = "",
          IncludeTrendSymbol = IncludeTrendSymbol
        ),
        ""
      )
      stars[is.na(stars)] <- ""
      table_data[[effect_column]] <- paste0(
        formatC(outcome_results$Estimate, digits = 3, format = "fg"),
        " (",
        formatC(outcome_results$ConfLow, digits = 3, format = "fg"),
        ", ",
        formatC(outcome_results$ConfHigh, digits = 3, format = "fg"),
        ")",
        stars
      )
      table_data[[p_column]] <- dplyr::case_when(
        is.na(outcome_results$PValue) ~ NA_character_,
        outcome_results$PValue < 0.001 ~ "<0.001",
        TRUE ~ formatC(outcome_results$PValue, digits = 2, format = "fg")
      )
      column_labels[[effect_column]] <- effect_label
      column_labels[[p_column]] <- "p-value"
      spanners[[column_prefix]] <- c(effect_column, p_column)
      significant_rows[[effect_column]] <- which(outcome_results$Significant)
      significant_rows[[p_column]] <- which(outcome_results$Significant)
    } else {
      large_columns <- c("N", "Estimate", "StdError", "ConfLow", "ConfHigh", "PValue")
      large_labels <- c(
        N = "N",
        Estimate = if (identical(effect_type, "Odds ratio")) "Odds ratio" else "Estimate",
        StdError = "SE",
        ConfLow = "95% CI Low",
        ConfHigh = "95% CI High",
        PValue = "p-value"
      )
      outcome_columns <- paste0(column_prefix, "_", large_columns)
      for (column_index in seq_along(large_columns)) {
        table_data[[outcome_columns[[column_index]]]] <- outcome_results[[large_columns[[column_index]]]]
      }
      column_labels <- c(column_labels, as.list(stats::setNames(large_labels, outcome_columns)))
      spanners[[column_prefix]] <- outcome_columns
      numeric_columns <- c(numeric_columns, outcome_columns)
    }

    attr(spanners[[column_prefix]], "label") <- outcome_label
  }

  out <- gt::gt(table_data, rowname_col = "Variable") %>%
    gt::cols_label(.list = column_labels)

  for (spanner_id in names(spanners)) {
    out <- gt::tab_spanner(
      out,
      label = attr(spanners[[spanner_id]], "label"),
      columns = spanners[[spanner_id]],
      id = spanner_id
    )
  }

  if (formatted) {
    for (column_name in names(significant_rows)) {
      if (length(significant_rows[[column_name]]) > 0) {
        out <- gt::tab_style(
          out,
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_body(
            columns = column_name,
            rows = significant_rows[[column_name]]
          )
        )
      }
    }
  } else {
    out <- gt::fmt_number(out, columns = numeric_columns, decimals = 3)
  }

  out
}

ScidrUnivariateOutcomeFamily <- function(x, outcome, method) {
  if (method == "lm") {
    return("linear")
  }
  if (method == "logistic") {
    ScidrUnivariateValidateBinaryOutcome(x, outcome)
    return("logistic")
  }

  non_missing <- x[!is.na(x)]
  if (is.numeric(non_missing)) {
    return("linear")
  }
  if (is.logical(non_missing)) {
    return("logistic")
  }
  if (is.factor(non_missing)) {
    n_levels <- length(levels(droplevels(non_missing)))
    if (n_levels == 2) {
      return("logistic")
    }
    if (n_levels > 2) {
      stop(
        "Outcome ", outcome, " has ", n_levels, " levels. ",
        "Multinomial regression is not yet supported; logistic regression requires a two-level outcome."
      )
    }
  }
  if (is.character(non_missing)) {
    n_levels <- length(unique(non_missing))
    if (n_levels == 2) {
      return("logistic")
    }
    if (n_levels > 2) {
      stop(
        "Outcome ", outcome, " has ", n_levels, " levels. ",
        "Multinomial regression is not yet supported; logistic regression requires a two-level outcome."
      )
    }
  }

  stop("Outcome ", outcome, " must be numeric, logical, or a two-level factor or character vector.")
}

ScidrUnivariateValidateBinaryOutcome <- function(x, outcome) {
  non_missing <- x[!is.na(x)]
  n_levels <- if (is.factor(non_missing)) {
    length(levels(droplevels(non_missing)))
  } else {
    length(unique(non_missing))
  }

  if (n_levels != 2) {
    stop("Outcome ", outcome, " must have exactly two levels for logistic regression.")
  }
}

ScidrUnivariatePrepareOutcome <- function(x, outcome, outcome_family) {
  if (outcome_family == "linear") {
    return(list(
      Value = x,
      ReferenceLevel = NA_character_,
      EventLevel = NA_character_
    ))
  }

  ScidrUnivariateValidateBinaryOutcome(x, outcome)
  if (is.logical(x)) {
    value <- as.integer(x)
    return(list(
      Value = value,
      ReferenceLevel = "FALSE",
      EventLevel = "TRUE"
    ))
  }

  value <- if (is.factor(x)) {
    droplevels(x)
  } else {
    factor(x)
  }
  outcome_levels <- levels(value)

  list(
    Value = value,
    ReferenceLevel = outcome_levels[[1]],
    EventLevel = outcome_levels[[2]]
  )
}
