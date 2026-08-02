#' Multivariable regression table
#'
#' Fit one multivariable regression model per outcome and return a stable,
#' label-aware regression object for downstream tables, diagnostics, and plots.
#'
#' @param data Data frame containing outcomes, predictors, and covariates.
#' @param outcome_vars Character vector of outcome variable names.
#' @param predictor_vars Character vector of predictor variable names.
#' @param covariates Optional character vector of covariate variable names.
#'   Covariates are treated as mandatory adjustments: for penalized methods
#'   (`"ridge"`, `"lasso"`, `"elasticnet"`) they are exempted from the penalty
#'   (`penalty.factor = 0`), so they are never shrunk or selected out of the
#'   model.
#' @param Standardize Logical. If `TRUE`, ordinary models are fit on
#'   standardized continuous variables for the primary estimate. Standardized
#'   coefficients are always calculated separately regardless of this setting.
#' @param Relabel Logical. If `TRUE`, use variable labels from `sjlabelled`
#'   when available.
#' @param TreatOrdinalAs How ordinal outcomes and predictors are handled.
#' @param FDR Logical. If `TRUE`, calculate FDR-adjusted p-values for ordinary
#'   regression terms.
#' @param FDRAlpha Numeric FDR threshold retained in metadata.
#' @param Method Regression method. One of `"lm"`, `"ridge"`, `"lasso"`, or
#'   `"elasticnet"`.
#' @param CVFolds Number of cross-validation folds for penalized models.
#' @param Lambda Lambda selection rule for penalized models. One of
#'   `"lambda.min"` or `"lambda.1se"`.
#' @param Seed Random seed used for deterministic cross-validation folds.
#' @param MissingDataStrategy Missing-data handling strategy. The default,
#'   `"drop_sparse_impute"`, drops sparse predictors and covariates, then
#'   imputes remaining predictor missingness.
#' @param MaxMissingPredictor Maximum allowed missingness proportion for
#'   predictors and covariates before they are dropped by sparse-drop
#'   strategies. Default is `0.30`.
#' @param ImputeMethod Imputation method for predictor/covariate missingness.
#'   Currently `"median_mode"`: median for numeric variables and mode for
#'   factor, character, and logical variables.
#' @param MinCompleteCases Optional minimum number of modeling rows required
#'   after missing-data handling.
#' @param outcome_modes Multi-category outcome strategy. Supply a single
#'   `"auto"` or a named character vector whose values are `"auto"`,
#'   `"multinomial"`, `"ordinal"`, `"one_vs_rest"`, `"binary_subset"`, or
#'   `"skip"`. In automatic mode, ordered factors use proportional-odds
#'   regression and unordered factors use multinomial regression.
#' @param reference_levels Optional named character vector giving reference
#'   levels for categorical outcomes. Unspecified outcomes use their first
#'   retained factor level.
#' @param binary_subsets Optional named list. Each outcome assigned
#'   `"binary_subset"` must have exactly two level names here, ordered as
#'   reference then event.
#'
#' @return A named list with stable components: `Models`, `FormattedTable`,
#'   `LargeTable`, `RegressionMatrix`, `VariableImportanceMatrix`,
#'   `Predictions`, `Diagnostics`, `ModelSummary`, `Multicollinearity`,
#'   `Plots`, and `Metadata`. `FormattedTable` is a report-facing `gt` table
#'   grouped by outcome, matching the style of
#'   [MakeUnivariateRegressionTable()]: predictor rows only, a combined
#'   `Estimate (95% CI)` cell, and bold significant p-values. `LargeTable` is a
#'   data frame holding the full per-term detail (including covariate rows) for
#'   programmatic use, plus an `Aliased` flag marking perfectly collinear terms
#'   the model dropped. `ModelSummary` reports per-outcome `Converged`,
#'   `SeparationDetected`, and `AliasedTermCount`. For ordinary (`"lm"`)
#'   logistic fits, quasi-complete separation is detected (fitted probabilities
#'   pinned at 0/1, exploded standardized coefficients, or non-convergence); the
#'   affected model's estimates are blanked (`NA`) and `Converged` is set to
#'   `FALSE` so unreliable coefficients do not propagate into tables or plots.
#'   `ModelSummary` also carries an omnibus model test per outcome
#'   (`ModelStat`, `ModelStatType`, `ModelPValue`): an F-test for linear models
#'   and a likelihood-ratio test for logistic models (`NA` for penalized fits,
#'   which have no valid classical omnibus test). `Plots` contains ggplot
#'   objects built from the stored result tables and predictions without
#'   refitting models; the coefficient heatmap uses robust, clamped fill limits
#'   so a single extreme value cannot dominate the scale, and each outcome
#'   column is annotated at the top with its omnibus p-value (ordinary models)
#'   or cross-validated deviance explained (penalized models) to discourage
#'   interpreting coefficients from a model that is not significant overall.
#'
#'   Multi-category outcomes add explicit `OutcomeLevel`, `ReferenceLevel`,
#'   `Contrast`, `ComparisonLabel`, and `OutcomeMode` fields. Unordered factors
#'   use nominal multinomial models by default. Ordered factors use
#'   proportional-odds models; their odds ratios describe movement toward a
#'   higher category, conditional on the predictors. One-vs-rest models are
#'   available for level-specific scientific questions, but their overlapping
#'   comparisons should be interpreted with multiplicity in mind. Binary
#'   subsets change both the analysis population and estimand. The resolved
#'   strategy, reference, engine, class counts, and concise scientific advice
#'   are recorded under `Metadata$Outcomes` and `Metadata$ModelingAdvice`.
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param OutcomeVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @param PredictorVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param Covars \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' result <- MultivariableRegressionTable(
#'   Labelled,
#'   outcome_vars = "AXL",
#'   predictor_vars = c("Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin"),
#'   covariates = "age"
#' )
#'
#' # Display the regression coefficient matrix plot
#' result$Plots$RegressionMatrix
#'
#' # Nominal outcome: each non-reference level versus the named reference.
#' ExampleData <- Labelled
#' ExampleData$Race <- factor(rep(c("Asian", "Black", "White"), length.out = nrow(ExampleData)))
#' nominal_result <- MultivariableRegressionTable(
#'   ExampleData,
#'   outcome_vars = "Race",
#'   predictor_vars = c("Adiponectin", "Alpha_1_Antitrypsin"),
#'   Method = "lasso",
#'   reference_levels = c(Race = "White")
#' )
#'
#' # Ordered outcome: one proportional-odds effect per predictor.
#' ExampleData$Severity <- ordered(
#'   rep(c("Mild", "Moderate", "Severe"), length.out = nrow(ExampleData)),
#'   levels = c("Mild", "Moderate", "Severe")
#' )
#' ordinal_result <- MultivariableRegressionTable(
#'   ExampleData,
#'   outcome_vars = "Severity",
#'   predictor_vars = c("Adiponectin", "Alpha_1_Antitrypsin"),
#'   Method = "lm"
#' )
#' }
#' @export
MultivariableRegressionTable <- function(data,
    outcome_vars,
    predictor_vars,
    covariates = NULL,
    Standardize = TRUE,
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    FDR = TRUE,
    FDRAlpha = 0.05,
    Method = c("lm", "ridge", "lasso", "elasticnet"),
    CVFolds = 10,
    Lambda = c("lambda.min", "lambda.1se"),
    Seed = 123,
    MissingDataStrategy = c("drop_sparse_impute", "impute", "complete_cases", "drop_sparse_complete_cases"),
    MaxMissingPredictor = 0.3,
    ImputeMethod = c("median_mode"),
    MinCompleteCases = NULL,
    outcome_modes = "auto",
    reference_levels = NULL,
    binary_subsets = NULL,
    Data = lifecycle::deprecated(),
    OutcomeVars = lifecycle::deprecated(),
    PredictorVars = lifecycle::deprecated(),
    Covars = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "MultivariableRegressionTable(Data)", "MultivariableRegressionTable(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(OutcomeVars)) {
    lifecycle::deprecate_warn("19.15.0", "MultivariableRegressionTable(OutcomeVars)", "MultivariableRegressionTable(outcome_vars)")
    outcome_vars <- OutcomeVars
  }
  if (!missing(outcome_vars)) OutcomeVars <- outcome_vars
  if (lifecycle::is_present(PredictorVars)) {
    lifecycle::deprecate_warn("19.15.0", "MultivariableRegressionTable(PredictorVars)", "MultivariableRegressionTable(predictor_vars)")
    predictor_vars <- PredictorVars
  }
  if (!missing(predictor_vars)) PredictorVars <- predictor_vars
  if (lifecycle::is_present(Covars)) {
    lifecycle::deprecate_warn("19.15.0", "MultivariableRegressionTable(Covars)", "MultivariableRegressionTable(covariates)")
    covariates <- Covars
  }
  Covars <- covariates


  Method <- match.arg(Method)
  Lambda <- match.arg(Lambda)
  MissingDataStrategy <- match.arg(MissingDataStrategy)
  ImputeMethod <- match.arg(ImputeMethod)

  valid_outcome_modes <- c("auto", "multinomial", "ordinal", "one_vs_rest", "binary_subset", "skip")
  if (!is.character(outcome_modes) || length(outcome_modes) == 0 ||
      any(!outcome_modes %in% valid_outcome_modes)) {
    stop("outcome_modes must contain only: ", paste(valid_outcome_modes, collapse = ", "), ".")
  }
  if (length(outcome_modes) > 1 && (is.null(names(outcome_modes)) || any(names(outcome_modes) == ""))) {
    stop("Multiple outcome_modes values must be named by outcome variable.")
  }
  if (!is.null(reference_levels) &&
      (!is.character(reference_levels) || is.null(names(reference_levels)) || any(names(reference_levels) == ""))) {
    stop("reference_levels must be NULL or a named character vector.")
  }
  if (!is.null(binary_subsets) && (!is.list(binary_subsets) || is.null(names(binary_subsets)))) {
    stop("binary_subsets must be NULL or a named list.")
  }

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
  if (!is.logical(Relabel) || length(Relabel) != 1) {
    stop("Relabel must be TRUE or FALSE.")
  }
  if (!is.logical(FDR) || length(FDR) != 1) {
    stop("FDR must be TRUE or FALSE.")
  }
  if (!is.numeric(FDRAlpha) || length(FDRAlpha) != 1 || is.na(FDRAlpha) || FDRAlpha <= 0 || FDRAlpha >= 1) {
    stop("FDRAlpha must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(CVFolds) || length(CVFolds) != 1 || is.na(CVFolds) || CVFolds < 2) {
    stop("CVFolds must be a single numeric value of at least 2.")
  }
  if (!is.numeric(Seed) || length(Seed) != 1 || is.na(Seed)) {
    stop("Seed must be a single numeric value.")
  }
  if (!is.numeric(MaxMissingPredictor) ||
      length(MaxMissingPredictor) != 1 ||
      is.na(MaxMissingPredictor) ||
      MaxMissingPredictor < 0 ||
      MaxMissingPredictor > 1) {
    stop("MaxMissingPredictor must be a single numeric value between 0 and 1.")
  }
  if (!is.null(MinCompleteCases) &&
      (!is.numeric(MinCompleteCases) ||
       length(MinCompleteCases) != 1 ||
       is.na(MinCompleteCases) ||
       MinCompleteCases < 1)) {
    stop("MinCompleteCases must be NULL or a single positive number.")
  }

  all_model_vars <- unique(c(OutcomeVars, PredictorVars, Covars))
  missing_vars <- setdiff(all_model_vars, names(Data))
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in Data: ", paste(missing_vars, collapse = ", "))
  }
  ordinal_treatment <- match.arg(TreatOrdinalAs, c("Categorical", "Continuous", "Both", "Exclude"))
  if (ordinal_treatment == "Both") {
    stop("TreatOrdinalAs = 'Both' is not meaningful for MultivariableRegressionTable().", call. = FALSE)
  }
  ordinal_reference <- ConvertOrdinalToNumeric(
    Data, all_model_vars, TreatOrdinalAs = "Categorical", ReturnMetadata = TRUE
  )
  if (ordinal_treatment == "Exclude" &&
      length(ordinal_reference$ordinal_variables)) {
    stop("TreatOrdinalAs = 'Exclude' cannot be used when ordinal model variables are explicitly supplied.", call. = FALSE)
  }
  Data <- ConvertOrdinalToNumeric(
    Data, all_model_vars, TreatOrdinalAs = ordinal_treatment, ReturnMetadata = TRUE
  )$data
  unknown_mode_outcomes <- setdiff(names(outcome_modes), OutcomeVars)
  unknown_reference_outcomes <- setdiff(names(reference_levels), OutcomeVars)
  unknown_subset_outcomes <- setdiff(names(binary_subsets), OutcomeVars)
  if (length(unknown_mode_outcomes) > 0) stop("outcome_modes names not found in outcome_vars: ", paste(unknown_mode_outcomes, collapse = ", "))
  if (length(unknown_reference_outcomes) > 0) stop("reference_levels names not found in outcome_vars: ", paste(unknown_reference_outcomes, collapse = ", "))
  if (length(unknown_subset_outcomes) > 0) stop("binary_subsets names not found in outcome_vars: ", paste(unknown_subset_outcomes, collapse = ", "))

  multicategory_outcomes <- OutcomeVars[vapply(OutcomeVars, function(outcome) {
    x <- Data[[outcome]]
    is.factor(x) && nlevels(droplevels(x[!is.na(x)])) > 2
  }, logical(1))]
  explicitly_multicategory <- intersect(
    names(outcome_modes)[outcome_modes != "auto"],
    OutcomeVars
  )
  explicit_references <- intersect(names(reference_levels), OutcomeVars)
  if (length(union(union(multicategory_outcomes, explicitly_multicategory), explicit_references)) > 0) {
    return(ScidrMulticategoryRegressionDispatch(
      Data = Data,
      OutcomeVars = OutcomeVars,
      PredictorVars = PredictorVars,
      Covars = Covars,
      Standardize = Standardize,
      Relabel = Relabel,
      FDR = FDR,
      FDRAlpha = FDRAlpha,
      Method = Method,
      CVFolds = CVFolds,
      Lambda = Lambda,
      Seed = Seed,
      MissingDataStrategy = MissingDataStrategy,
      MaxMissingPredictor = MaxMissingPredictor,
      ImputeMethod = ImputeMethod,
      MinCompleteCases = MinCompleteCases,
      outcome_modes = outcome_modes,
      reference_levels = reference_levels,
      binary_subsets = binary_subsets
    ))
  }

  if (Method != "lm" && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for ridge, lasso, and elasticnet methods.")
  }

  CVFolds <- as.integer(CVFolds)
  Seed <- as.integer(Seed)
  if (!is.null(MinCompleteCases)) {
    MinCompleteCases <- as.integer(MinCompleteCases)
  }
  model_terms <- unique(c(PredictorVars, Covars))
  label_lookup <- ScidrRegressionLabels(Data, unique(c(OutcomeVars, model_terms)), Relabel)
  outcome_families <- stats::setNames(
    vapply(OutcomeVars, function(outcome) ScidrOutcomeFamily(Data[[outcome]], outcome), character(1)),
    OutcomeVars
  )
  global_missingness <- ScidrPredictorMissingnessSummary(Data, model_terms)
  globally_dropped_terms <- if (MissingDataStrategy %in% c("drop_sparse_impute", "drop_sparse_complete_cases")) {
    global_missingness$Variable[global_missingness$MissingProportion > MaxMissingPredictor]
  } else {
    character(0)
  }
  retained_model_terms <- setdiff(model_terms, globally_dropped_terms)

  if (length(intersect(PredictorVars, retained_model_terms)) == 0) {
    stop(
      "All PredictorVars were dropped for missingness. ",
      "MaxMissingPredictor = ", MaxMissingPredictor, ". ",
      "Dropped predictors: ", paste(intersect(PredictorVars, globally_dropped_terms), collapse = ", "), ". ",
      "Increase MaxMissingPredictor, use MissingDataStrategy = 'impute', or reduce PredictorVars."
    )
  }

  multicollinearity <- ScidrRegressionMulticollinearity(Data, retained_model_terms)

  model_list <- list()
  large_tables <- list()
  prediction_tables <- list()
  diagnostics_tables <- list()
  missing_tables <- list()
  tuning_list <- list()

  for (outcome_index in seq_along(OutcomeVars)) {
    outcome <- OutcomeVars[[outcome_index]]
    outcome_family <- outcome_families[[outcome]]
    missing_info <- ScidrPrepareRegressionModelData(
      Data = Data,
      outcome = outcome,
      model_terms = model_terms,
      retained_model_terms = retained_model_terms,
      dropped_terms = globally_dropped_terms,
      missing_data_strategy = MissingDataStrategy,
      impute_method = ImputeMethod,
      max_missing_predictor = MaxMissingPredictor
    )
    df_model <- missing_info$ModelData
    complete_rows <- missing_info$ModelRows
    model_terms_outcome <- missing_info$ModelTerms
    predictor_vars_outcome <- intersect(PredictorVars, model_terms_outcome)
    missing_removed <- missing_info$MissingRemoved
    percent_removed <- missing_info$PercentRemoved
    min_rows_required <- if (is.null(MinCompleteCases)) length(model_terms_outcome) + 2 else MinCompleteCases

    if (nrow(df_model) < min_rows_required) {
      stop(ScidrNotEnoughRowsMessage(outcome, missing_info, min_rows_required))
    }

    if (outcome_family == "logistic") {
      outcome_info <- ScidrPrepareBinaryOutcome(df_model[[outcome]], outcome)
      df_model[[outcome]] <- outcome_info$Value
      remaining_levels <- unique(df_model[[outcome]][!is.na(df_model[[outcome]])])
      if (length(remaining_levels) < 2) {
        stop(
          "Outcome ", outcome, " has only one level after missing-data handling: ",
          paste(remaining_levels, collapse = ", "), ". ",
          "Use a less restrictive missing-data strategy or check the outcome coding."
        )
      }
      if (!requireNamespace("pROC", quietly = TRUE)) {
        stop("Package 'pROC' is required for logistic regression diagnostics.")
      }
    } else {
      outcome_info <- list(Levels = NA_character_)
    }

    fit_result <- if (Method == "lm") {
      ScidrFitOrdinaryRegression(
        df_model = df_model,
        outcome = outcome,
        model_terms = model_terms_outcome,
        predictor_vars = predictor_vars_outcome,
        outcome_family = outcome_family,
        standardize_estimates = Standardize
      )
    } else {
      ScidrFitPenalizedRegression(
        df_model = df_model,
        outcome = outcome,
        model_terms = model_terms_outcome,
        predictor_vars = predictor_vars_outcome,
        outcome_family = outcome_family,
        method = Method,
        cv_folds = CVFolds,
        lambda_choice = Lambda,
        seed = Seed + outcome_index
      )
    }

    fit_result$TermTable$OutcomeIndex <- outcome_index
    fit_result$TermTable$PredictorIndex <- match(fit_result$TermTable$Predictor, PredictorVars)
    fit_result$TermTable$Outcome <- outcome
    fit_result$TermTable$OutcomeLabel <- unname(label_lookup[[outcome]])
    fit_result$TermTable$OutcomeFamily <- outcome_family
    fit_result$TermTable$PredictorLabel <- unname(label_lookup[fit_result$TermTable$Predictor])
    fit_result$TermTable$RegressionMethod <- Method
    fit_result$TermTable$MissingRemoved <- missing_removed
    fit_result$TermTable$PercentRemoved <- percent_removed
    fit_result$TermTable$MissingDataStrategy <- MissingDataStrategy
    fit_result$TermTable$Imputed <- fit_result$TermTable$Predictor %in% missing_info$ImputedVariables
    fit_result$TermTable$DroppedForMissingness <- FALSE

    predictions <- ScidrRegressionPredictions(
      fit_result = fit_result,
      df_model = df_model,
      outcome = outcome,
      outcome_family = outcome_family,
      full_row_count = nrow(Data),
      complete_rows = complete_rows
    )

    diagnostics <- ScidrRegressionDiagnostics(
      fit_result = fit_result,
      df_model = df_model,
      outcome = outcome,
      outcome_family = outcome_family,
      method = Method,
      predictions_complete = predictions[predictions$.ModelComplete, , drop = FALSE]
    )

    diagnostics$Outcome <- outcome
    diagnostics$OutcomeLabel <- unname(label_lookup[[outcome]])
    diagnostics$OutcomeFamily <- outcome_family
    diagnostics$RegressionMethod <- Method
    diagnostics$MissingRemoved <- missing_removed
    diagnostics$PercentRemoved <- percent_removed
    diagnostics$MissingDataStrategy <- MissingDataStrategy
    diagnostics$DroppedPredictors <- I(list(intersect(PredictorVars, missing_info$DroppedVariables)))
    diagnostics$ImputedPredictors <- I(list(intersect(PredictorVars, missing_info$ImputedVariables)))
    diagnostics$DroppedPredictorCount <- length(intersect(PredictorVars, missing_info$DroppedVariables))
    diagnostics$ImputedPredictorCount <- length(intersect(PredictorVars, missing_info$ImputedVariables))

    fit_result$TermTable$SampleSize <- diagnostics$SampleSize[[1]]

    model_list[[outcome]] <- fit_result$Model
    large_tables[[outcome]] <- fit_result$TermTable
    prediction_tables[[outcome]] <- predictions
    diagnostics_tables[[outcome]] <- diagnostics
    missing_tables[[outcome]] <- data.frame(
      Outcome = outcome,
      MissingRemoved = missing_removed,
      PercentRemoved = percent_removed,
      OriginalN = nrow(Data),
      OutcomeNonMissingN = missing_info$OutcomeNonMissingN,
      CompleteCaseN = missing_info$CompleteCaseN,
      FinalN = nrow(df_model),
      MissingDataStrategy = MissingDataStrategy,
      DroppedVariables = I(list(missing_info$DroppedVariables)),
      ImputedVariables = I(list(missing_info$ImputedVariables)),
      stringsAsFactors = FALSE
    )
    tuning_list[[outcome]] <- fit_result$Tuning
  }

  large_table <- dplyr::bind_rows(large_tables)
  large_table <- ScidrAnnotateTableDefaults(large_table)
  if (Method == "lm" && FDR) {
    adjust_rows <- !is.na(large_table$PValue)
    large_table$FDR <- NA_real_
    large_table$FDR[adjust_rows] <- stats::p.adjust(large_table$PValue[adjust_rows], method = "fdr")
  } else {
    large_table$FDR <- NA_real_
  }

  star_source <- if (FDR) large_table$FDR else large_table$PValue
  large_table$Stars <- ScidrPValueStars(star_source)
  large_table$HoverText <- ScidrRegressionHoverText(large_table)

  reporting_rows <- large_table$Predictor %in% PredictorVars
  formatted_reporting <- large_table[reporting_rows, , drop = FALSE]
  formatted_table <- ScidrMultivariableGtTable(formatted_reporting, formatted = TRUE)
  regression_matrix <- large_table[reporting_rows, c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel",
    "ReferenceLevel", "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel",
    "PredictorLabel", "OutcomeFamily", "Estimate", "Effect", "EffectType",
    "StandardizedBeta", "PValue", "FDR", "Stars", "Selected",
    "VariableImportance", "VariableImportanceType", "HoverText"
  ), drop = FALSE]
  variable_importance_matrix <- large_table[reporting_rows, c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel",
    "ReferenceLevel", "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel",
    "PredictorLabel", "OutcomeFamily", "VariableImportance",
    "VariableImportanceType", "Selected", "HoverText"
  ), drop = FALSE]

  predictions <- dplyr::bind_rows(prediction_tables)
  predictions$.ModelComplete <- NULL

  diagnostics <- dplyr::bind_rows(diagnostics_tables)
  diagnostics <- ScidrAnnotateTableDefaults(diagnostics)
  model_summary <- diagnostics[, c(
    "Outcome", "OutcomeLabel", "OutcomeFamily", "RegressionMethod", "SampleSize",
    "MissingRemoved", "PercentRemoved", "PredictorCount", "Converged",
    "SeparationDetected", "AliasedTermCount",
    "ModelStat", "ModelStatType", "ModelPValue",
    "R2", "AdjustedR2", "AUC", "McFaddenR2", "RMSE", "AIC", "BIC",
    "DevianceExplained", "DroppedPredictorCount", "ImputedPredictorCount"
  ), drop = FALSE]
  model_summary <- ScidrAnnotateTableDefaults(model_summary)
  model_summary$MaximumVIF <- multicollinearity$MaximumVIF
  model_summary$MaximumCorrelation <- multicollinearity$MaximumCorrelation

  metadata <- list(
    AnalysisSettings = list(
      RegressionMethod = Method,
      OutcomeVars = OutcomeVars,
      PredictorVars = PredictorVars,
      Covars = Covars,
      Standardize = Standardize,
      FDR = FDR,
      FDRAlpha = FDRAlpha,
      CVFolds = CVFolds,
      Lambda = Lambda,
      Seed = Seed,
      MissingDataStrategy = MissingDataStrategy,
      MaxMissingPredictor = MaxMissingPredictor,
      ImputeMethod = ImputeMethod,
      MinCompleteCases = MinCompleteCases,
      ModelFamilies = outcome_families,
      Tuning = tuning_list
    ),
    Missingness = list(
      Outcomes = dplyr::bind_rows(missing_tables),
      Predictors = global_missingness
    ),
    Session = utils::sessionInfo(),
    PackageVersion = as.character(utils::packageVersion("SciDataReportR")),
    FunctionCall = match.call()
  )
  metadata$Outcomes <- dplyr::bind_rows(lapply(OutcomeVars, function(outcome) {
    x <- Data[[outcome]]
    retained_levels <- if (is.factor(x)) levels(droplevels(x[!is.na(x)])) else if (is.logical(x)) c("FALSE", "TRUE") else character(0)
    mode <- if (outcome_families[[outcome]] == "linear") "linear" else "binary"
    ScidrOutcomeMetadata(
      outcome, x, list(Requested = "auto"), mode,
      if (length(retained_levels) > 0) retained_levels[[1]] else NA_character_,
      retained_levels,
      if (Method == "lm") if (mode == "linear") "stats::lm" else "stats::glm" else "glmnet::cv.glmnet",
      Data
    )
  }))
  metadata$ModelingAdvice <- ScidrMulticategoryModelingAdvice()

  plots <- ScidrRegressionPlots(
    regression_matrix = regression_matrix,
    variable_importance_matrix = variable_importance_matrix,
    predictions = predictions,
    model_summary = model_summary
  )

  return(list(
    Models = model_list,
    FormattedTable = formatted_table,
    LargeTable = as.data.frame(large_table),
    RegressionMatrix = as.data.frame(regression_matrix),
    VariableImportanceMatrix = as.data.frame(variable_importance_matrix),
    Predictions = as.data.frame(predictions),
    Diagnostics = as.data.frame(diagnostics),
    ModelSummary = as.data.frame(model_summary),
    Multicollinearity = multicollinearity,
    Plots = plots,
    Metadata = metadata
  ))
}

ScidrResolveOutcomeMode <- function(outcome, x, outcome_modes) {
  requested <- if (!is.null(names(outcome_modes)) && outcome %in% names(outcome_modes)) {
    unname(outcome_modes[[outcome]])
  } else if (length(outcome_modes) == 1 && is.null(names(outcome_modes))) {
    outcome_modes[[1]]
  } else {
    "auto"
  }
  retained <- x[!is.na(x)]
  if (is.factor(retained)) retained <- droplevels(retained)
  n_levels <- if (is.factor(retained)) nlevels(retained) else length(unique(retained))
  if (requested != "auto") {
    return(list(Requested = requested, Resolved = requested))
  }
  resolved <- if (is.numeric(retained)) {
    "linear"
  } else if (is.logical(retained) || n_levels == 2) {
    "binary"
  } else if (is.ordered(retained)) {
    "ordinal"
  } else {
    "multinomial"
  }
  list(Requested = requested, Resolved = resolved)
}

ScidrMulticategoryRegressionDispatch <- function(Data,
                                                  OutcomeVars,
                                                  PredictorVars,
                                                  Covars,
                                                  Standardize,
                                                  Relabel,
                                                  FDR,
                                                  FDRAlpha,
                                                  Method,
                                                  CVFolds,
                                                  Lambda,
                                                  Seed,
                                                  MissingDataStrategy,
                                                  MaxMissingPredictor,
                                                  ImputeMethod,
                                                  MinCompleteCases,
                                                  outcome_modes,
                                                  reference_levels,
                                                  binary_subsets) {
  pieces <- list()
  outcome_metadata <- list()
  skipped <- character(0)

  for (outcome_index in seq_along(OutcomeVars)) {
    outcome <- OutcomeVars[[outcome_index]]
    x <- Data[[outcome]]
    mode_info <- ScidrResolveOutcomeMode(outcome, x, outcome_modes)
    mode <- mode_info$Resolved
    original_levels <- if (is.factor(x)) levels(x) else if (is.logical(x)) c("FALSE", "TRUE") else character(0)
    retained_levels <- if (is.factor(x)) levels(droplevels(x[!is.na(x)])) else original_levels
    reference <- if (!is.null(reference_levels) && outcome %in% names(reference_levels)) {
      unname(reference_levels[[outcome]])
    } else if (length(retained_levels) > 0) {
      retained_levels[[1]]
    } else {
      NA_character_
    }
    if (!is.na(reference) && is.factor(x) && !reference %in% retained_levels) {
      stop("Reference level '", reference, "' was not retained for outcome ", outcome, ".")
    }
    if (outcome %in% names(reference_levels) && mode %in% c("ordinal", "one_vs_rest")) {
      stop("reference_levels is not applicable to outcome '", outcome, "' in mode '", mode, "'.")
    }

    if (mode == "skip") {
      warning("Skipping multi-category outcome '", outcome, "' because outcome_modes requested 'skip'.", call. = FALSE)
      skipped <- c(skipped, outcome)
      outcome_metadata[[outcome]] <- ScidrOutcomeMetadata(
        outcome, x, mode_info, mode, reference, retained_levels, "none", Data
      )
      next
    }

    if (mode %in% c("linear", "binary")) {
      piece_data <- Data
      if (mode == "binary" && is.factor(x)) {
        piece_data[[outcome]] <- stats::relevel(droplevels(x), ref = reference)
        retained_levels <- levels(piece_data[[outcome]])
      }
      piece <- MultivariableRegressionTable(
        data = piece_data,
        outcome_vars = outcome,
        predictor_vars = PredictorVars,
        covariates = Covars,
        Standardize = Standardize,
        Relabel = Relabel,
        FDR = FALSE,
        FDRAlpha = FDRAlpha,
        Method = Method,
        CVFolds = CVFolds,
        Lambda = Lambda,
        Seed = Seed + outcome_index,
        MissingDataStrategy = MissingDataStrategy,
        MaxMissingPredictor = MaxMissingPredictor,
        ImputeMethod = ImputeMethod,
        MinCompleteCases = MinCompleteCases,
        outcome_modes = "auto"
      )
      if (mode == "binary") {
        event <- retained_levels[[2]]
        piece <- ScidrAnnotateRegressionPiece(piece, outcome, mode, event, reference,
          paste0(event, " vs ", reference), piece$LargeTable$OutcomeLabel[[1]])
      } else {
        piece <- ScidrAnnotateRegressionPiece(piece, outcome, mode, NA_character_, NA_character_,
          "Continuous", piece$LargeTable$OutcomeLabel[[1]])
      }
      pieces[[outcome]] <- piece
      outcome_metadata[[outcome]] <- ScidrOutcomeMetadata(
        outcome, x, mode_info, mode, reference, retained_levels,
        if (mode == "linear") "stats::lm" else if (Method == "lm") "stats::glm" else "glmnet::cv.glmnet",
        Data
      )
      next
    }

    if (!is.factor(x) || length(retained_levels) < 3) {
      stop("Outcome mode '", mode, "' requires a factor with at least three retained levels for outcome ", outcome, ".")
    }

    if (mode == "binary_subset") {
      subset_levels <- if (!is.null(binary_subsets) && outcome %in% names(binary_subsets)) binary_subsets[[outcome]] else NULL
      if (is.null(subset_levels) || length(subset_levels) != 2 || anyDuplicated(subset_levels) ||
          any(!subset_levels %in% retained_levels)) {
        stop("binary_subsets[['", outcome, "']] must contain exactly two distinct retained levels, ordered as reference then event.")
      }
      df_subset <- Data
      df_subset[[outcome]] <- factor(as.character(x), levels = subset_levels)
      piece <- MultivariableRegressionTable(
        data = df_subset, outcome_vars = outcome, predictor_vars = PredictorVars,
        covariates = Covars, Standardize = Standardize, Relabel = Relabel,
        FDR = FALSE, FDRAlpha = FDRAlpha, Method = Method, CVFolds = CVFolds,
        Lambda = Lambda, Seed = Seed + outcome_index,
        MissingDataStrategy = MissingDataStrategy, MaxMissingPredictor = MaxMissingPredictor,
        ImputeMethod = ImputeMethod, MinCompleteCases = MinCompleteCases,
        outcome_modes = "auto"
      )
      comparison <- paste0(subset_levels[[2]], " vs ", subset_levels[[1]])
      piece <- ScidrAnnotateRegressionPiece(piece, outcome, mode, subset_levels[[2]], subset_levels[[1]],
        comparison, piece$LargeTable$OutcomeLabel[[1]])
      pieces[[outcome]] <- piece
      outcome_metadata[[outcome]] <- ScidrOutcomeMetadata(
        outcome, x, mode_info, mode, subset_levels[[1]], subset_levels,
        if (Method == "lm") "stats::glm" else "glmnet::cv.glmnet", df_subset
      )
      next
    }

    if (mode == "one_vs_rest") {
      comparison_pieces <- list()
      for (level_index in seq_along(retained_levels)) {
        level <- retained_levels[[level_index]]
        df_binary <- Data
        df_binary[[outcome]] <- ifelse(is.na(x), NA, as.character(x) == level)
        piece <- MultivariableRegressionTable(
          data = df_binary, outcome_vars = outcome, predictor_vars = PredictorVars,
          covariates = Covars, Standardize = Standardize, Relabel = Relabel,
          FDR = FALSE, FDRAlpha = FDRAlpha, Method = Method, CVFolds = CVFolds,
          Lambda = Lambda, Seed = Seed + outcome_index + level_index,
          MissingDataStrategy = MissingDataStrategy, MaxMissingPredictor = MaxMissingPredictor,
          ImputeMethod = ImputeMethod, MinCompleteCases = MinCompleteCases,
          outcome_modes = "auto"
        )
        piece <- ScidrAnnotateRegressionPiece(piece, outcome, mode, level, "All other levels",
          paste0(level, " vs all others"), ScidrRegressionLabels(Data, outcome, Relabel)[[outcome]])
        comparison_pieces[[level]] <- piece
      }
      pieces[[outcome]] <- ScidrCombineRegressionPieces(comparison_pieces, FDR = FALSE, FDRAlpha = FDRAlpha)
      outcome_metadata[[outcome]] <- ScidrOutcomeMetadata(
        outcome, x, mode_info, mode, "All other levels", retained_levels,
        if (Method == "lm") "stats::glm (one-vs-rest)" else "glmnet::cv.glmnet (one-vs-rest)", Data
      )
      next
    }

    piece <- ScidrFitMulticategoryOutcome(
      Data = Data, outcome = outcome, PredictorVars = PredictorVars, Covars = Covars,
      Standardize = Standardize, Relabel = Relabel, Method = Method,
      CVFolds = CVFolds, Lambda = Lambda, Seed = Seed + outcome_index,
      MissingDataStrategy = MissingDataStrategy, MaxMissingPredictor = MaxMissingPredictor,
      ImputeMethod = ImputeMethod, MinCompleteCases = MinCompleteCases,
      mode = mode, reference = reference
    )
    pieces[[outcome]] <- piece
    outcome_metadata[[outcome]] <- ScidrOutcomeMetadata(
      outcome, x, mode_info, mode, if (mode == "ordinal") NA_character_ else reference, retained_levels,
      if (mode == "multinomial") {
        if (Method == "lm") "nnet::multinom" else "glmnet::cv.glmnet"
      } else if (Method == "lm") "MASS::polr" else "ordinalNet::ordinalNetTune",
      Data
    )
  }

  if (length(pieces) == 0) {
    stop("No outcomes remained after applying outcome_modes.")
  }
  out <- ScidrCombineRegressionPieces(pieces, FDR = FDR, FDRAlpha = FDRAlpha)
  out$Metadata$AnalysisSettings$OutcomeModes <- outcome_modes
  out$Metadata$AnalysisSettings$ReferenceLevels <- reference_levels
  out$Metadata$AnalysisSettings$BinarySubsets <- binary_subsets
  out$Metadata$Outcomes <- dplyr::bind_rows(outcome_metadata)
  out$Metadata$ModelingAdvice <- ScidrMulticategoryModelingAdvice()
  out$Metadata$SkippedOutcomes <- skipped
  out
}

ScidrOutcomeMetadata <- function(outcome, x, mode_info, resolved, reference, retained_levels, engine, data) {
  counts <- table(factor(as.character(x), levels = retained_levels), useNA = "no")
  data.frame(
    Outcome = outcome,
    RequestedMode = mode_info$Requested,
    ResolvedMode = resolved,
    Ordered = is.ordered(x),
    ReferenceLevel = ifelse(length(reference) == 0, NA_character_, reference),
    ContrastDirection = if (resolved == "ordinal") "Higher category vs lower category" else if (resolved == "one_vs_rest") "Level vs all other levels" else "Outcome level vs reference level",
    Engine = engine,
    OriginalLevels = I(list(if (is.factor(x)) levels(x) else if (is.logical(x)) c("FALSE", "TRUE") else character(0))),
    RetainedLevels = I(list(retained_levels)),
    LevelCounts = I(list(counts)),
    NonMissingN = sum(!is.na(x)),
    stringsAsFactors = FALSE
  )
}

ScidrFitMulticategoryOutcome <- function(Data,
                                         outcome,
                                         PredictorVars,
                                         Covars,
                                         Standardize,
                                         Relabel,
                                         Method,
                                         CVFolds,
                                         Lambda,
                                         Seed,
                                         MissingDataStrategy,
                                         MaxMissingPredictor,
                                         ImputeMethod,
                                         MinCompleteCases,
                                         mode,
                                         reference) {
  model_terms <- unique(c(PredictorVars, Covars))
  missingness <- ScidrPredictorMissingnessSummary(Data, model_terms)
  dropped <- if (MissingDataStrategy %in% c("drop_sparse_impute", "drop_sparse_complete_cases")) {
    missingness$Variable[missingness$MissingProportion > MaxMissingPredictor]
  } else character(0)
  retained_terms <- setdiff(model_terms, dropped)
  if (length(intersect(PredictorVars, retained_terms)) == 0) {
    stop("All PredictorVars were dropped for missingness for outcome ", outcome, ".")
  }
  missing_info <- ScidrPrepareRegressionModelData(
    Data, outcome, model_terms, retained_terms, dropped,
    MissingDataStrategy, ImputeMethod, MaxMissingPredictor
  )
  df_model <- missing_info$ModelData
  df_model[[outcome]] <- droplevels(df_model[[outcome]])
  retained_levels <- levels(df_model[[outcome]])
  if (length(retained_levels) < 3) {
    stop("Outcome ", outcome, " has fewer than three levels after missing-data handling.")
  }
  if (mode == "multinomial" && !reference %in% retained_levels) {
    stop("Reference level '", reference, "' was not retained for outcome ", outcome, ".")
  }
  if (mode == "multinomial") {
    df_model[[outcome]] <- stats::relevel(df_model[[outcome]], ref = reference)
  }
  retained_levels <- levels(df_model[[outcome]])
  min_rows_required <- if (is.null(MinCompleteCases)) length(retained_terms) + 2 else MinCompleteCases
  if (nrow(df_model) < min_rows_required) {
    stop(ScidrNotEnoughRowsMessage(outcome, missing_info, min_rows_required))
  }
  class_counts <- table(df_model[[outcome]])
  if (any(class_counts < 2)) {
    stop("Outcome ", outcome, " has a retained class with fewer than two observations: ",
      paste(names(class_counts), class_counts, sep = "=", collapse = ", "), ".")
  }
  if (Method != "lm" && min(class_counts) < CVFolds) {
    stop("CVFolds = ", CVFolds, " exceeds the smallest class size (", min(class_counts),
      ") for outcome ", outcome, ". Reduce CVFolds or combine sparse levels.")
  }

  fit <- if (mode == "multinomial") {
    if (Method == "lm") {
      ScidrFitOrdinaryMultinomial(df_model, outcome, retained_terms, PredictorVars, Standardize, reference)
    } else {
      ScidrFitPenalizedMultinomial(df_model, outcome, retained_terms, PredictorVars,
        Method, CVFolds, Lambda, Seed, reference)
    }
  } else if (Method == "lm") {
    ScidrFitOrdinaryOrdinal(df_model, outcome, retained_terms, PredictorVars, Standardize)
  } else {
    ScidrFitPenalizedOrdinal(df_model, outcome, retained_terms, PredictorVars,
      Method, CVFolds, Lambda, Seed)
  }

  label_lookup <- ScidrRegressionLabels(Data, unique(c(outcome, retained_terms)), Relabel)
  term_table <- fit$TermTable
  term_table$Outcome <- outcome
  term_table$OutcomeLabel <- unname(label_lookup[[outcome]])
  term_table$PredictorLabel <- unname(label_lookup[term_table$Predictor])
  term_table$PredictorIndex <- match(term_table$Predictor, PredictorVars)
  term_table$OutcomeIndex <- 1L
  term_table$OutcomeFamily <- mode
  term_table$OutcomeMode <- mode
  term_table$RegressionMethod <- Method
  term_table$MissingRemoved <- missing_info$MissingRemoved
  term_table$PercentRemoved <- missing_info$PercentRemoved
  term_table$MissingDataStrategy <- MissingDataStrategy
  term_table$Imputed <- term_table$Predictor %in% missing_info$ImputedVariables
  term_table$DroppedForMissingness <- FALSE
  term_table$SampleSize <- nrow(df_model)
  term_table$ComparisonLabel <- if (mode == "ordinal") {
    paste0(term_table$OutcomeLabel, ": higher vs lower")
  } else {
    paste0(term_table$OutcomeLabel, ": ", term_table$Contrast)
  }
  term_table$ComparisonIndex <- match(term_table$ComparisonLabel, unique(term_table$ComparisonLabel))
  term_table$FDR <- NA_real_
  term_table$Stars <- NA_character_
  term_table$HoverText <- ScidrRegressionHoverText(term_table)

  predictions <- ScidrMulticategoryPredictions(
    fit$Probabilities, df_model[[outcome]], outcome, mode,
    nrow(Data), missing_info$ModelRows
  )
  diagnostics <- ScidrMulticategoryDiagnostics(fit, df_model[[outcome]], predictions, outcome, mode, Method)
  diagnostics$OutcomeLabel <- unname(label_lookup[[outcome]])
  diagnostics$ComparisonLabel <- if (mode == "ordinal") {
    paste0(diagnostics$OutcomeLabel, ": higher vs lower")
  } else {
    paste0(diagnostics$OutcomeLabel, ": ", diagnostics$Contrast)
  }
  diagnostics$OutcomeMode <- mode
  diagnostics$OutcomeFamily <- mode
  diagnostics$RegressionMethod <- Method
  diagnostics$MissingRemoved <- missing_info$MissingRemoved
  diagnostics$PercentRemoved <- missing_info$PercentRemoved
  diagnostics$MissingDataStrategy <- MissingDataStrategy
  diagnostics$DroppedPredictors <- I(rep(list(intersect(PredictorVars, dropped)), nrow(diagnostics)))
  diagnostics$ImputedPredictors <- I(rep(list(intersect(PredictorVars, missing_info$ImputedVariables)), nrow(diagnostics)))
  diagnostics$DroppedPredictorCount <- length(intersect(PredictorVars, dropped))
  diagnostics$ImputedPredictorCount <- length(intersect(PredictorVars, missing_info$ImputedVariables))

  model_summary <- diagnostics[, c(
    "Outcome", "OutcomeLabel", "OutcomeLevel", "ReferenceLevel", "Contrast", "ComparisonLabel",
    "OutcomeMode", "OutcomeFamily", "RegressionMethod", "SampleSize", "MissingRemoved",
    "PercentRemoved", "PredictorCount", "Converged", "SeparationDetected", "AliasedTermCount",
    "ModelStat", "ModelStatType", "ModelPValue", "R2", "AdjustedR2", "AUC", "McFaddenR2",
    "RMSE", "AIC", "BIC", "DevianceExplained", "DroppedPredictorCount", "ImputedPredictorCount"
  ), drop = FALSE]

  reporting <- term_table$Predictor %in% PredictorVars
  regression_cols <- c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel", "ReferenceLevel",
    "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel", "PredictorLabel",
    "OutcomeFamily", "Estimate", "Effect", "EffectType", "StandardizedBeta", "PValue", "FDR",
    "Stars", "Selected", "VariableImportance", "VariableImportanceType", "HoverText"
  )
  importance_cols <- c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel", "ReferenceLevel",
    "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel", "PredictorLabel",
    "OutcomeFamily", "VariableImportance", "VariableImportanceType", "Selected", "HoverText"
  )
  multicollinearity <- ScidrRegressionMulticollinearity(Data, retained_terms)
  metadata <- list(
    AnalysisSettings = list(
      RegressionMethod = Method, OutcomeVars = outcome, PredictorVars = PredictorVars,
      Covars = Covars, Standardize = Standardize, FDR = FALSE, CVFolds = CVFolds,
      Lambda = Lambda, Seed = Seed, MissingDataStrategy = MissingDataStrategy,
      MaxMissingPredictor = MaxMissingPredictor, ImputeMethod = ImputeMethod,
      MinCompleteCases = MinCompleteCases, Tuning = setNames(list(fit$Tuning), outcome)
    ),
    Missingness = list(
      Outcomes = data.frame(
        Outcome = outcome, MissingRemoved = missing_info$MissingRemoved,
        PercentRemoved = missing_info$PercentRemoved, OriginalN = nrow(Data),
        OutcomeNonMissingN = missing_info$OutcomeNonMissingN,
        CompleteCaseN = missing_info$CompleteCaseN, FinalN = nrow(df_model),
        MissingDataStrategy = MissingDataStrategy,
        DroppedVariables = I(list(dropped)), ImputedVariables = I(list(missing_info$ImputedVariables))
      ),
      Predictors = missingness
    ),
    Session = utils::sessionInfo(), FunctionCall = match.call()
  )
  regression_matrix <- term_table[reporting, regression_cols, drop = FALSE]
  importance_matrix <- term_table[reporting, importance_cols, drop = FALSE]
  list(
    Models = setNames(list(fit$Model), outcome),
    FormattedTable = ScidrMultivariableGtTable(term_table[reporting, , drop = FALSE], TRUE),
    LargeTable = as.data.frame(term_table), RegressionMatrix = as.data.frame(regression_matrix),
    VariableImportanceMatrix = as.data.frame(importance_matrix), Predictions = as.data.frame(predictions),
    Diagnostics = as.data.frame(diagnostics), ModelSummary = as.data.frame(model_summary),
    Multicollinearity = multicollinearity,
    Plots = ScidrRegressionPlots(regression_matrix, importance_matrix, predictions, model_summary),
    Metadata = metadata
  )
}

ScidrMulticategoryTermRow <- function(term, estimate, std_estimate, se, p_value,
                                      outcome_level, reference_level, contrast,
                                      selected, lambda = NA_real_, alpha = NA_real_,
                                      converged = TRUE, warnings = "") {
  data.frame(
    Predictor = term, Estimate = estimate, RawEstimate = estimate,
    Effect = exp(estimate), EffectType = "Odds Ratio", StandardizedBeta = std_estimate,
    StandardError = se, TestStatistic = ifelse(is.na(se) || se == 0, NA_real_, estimate / se),
    DegreesFreedom = NA_real_, LowerCI = ifelse(is.na(se), NA_real_, exp(estimate - 1.96 * se)),
    UpperCI = ifelse(is.na(se), NA_real_, exp(estimate + 1.96 * se)), PValue = p_value,
    FDR = NA_real_, VariableImportance = abs(std_estimate),
    VariableImportanceType = "Absolute Standardized Beta", Selected = selected,
    Aliased = is.na(estimate), Lambda = lambda, Alpha = alpha, Converged = converged,
    Warnings = warnings, OutcomeLevel = outcome_level, ReferenceLevel = reference_level,
    Contrast = contrast, stringsAsFactors = FALSE
  )
}

ScidrQuoteFormulaNames <- function(var_names) {
  if (length(var_names) == 0) return(character(0))
  vapply(var_names, function(var_name) deparse1(as.name(var_name), backtick = TRUE), character(1), USE.NAMES = FALSE)
}

ScidrRegressionFormula <- function(model_terms, outcome = NULL) {
  term_symbols <- lapply(model_terms, as.name)
  rhs <- if (length(term_symbols) == 0) {
    1
  } else {
    Reduce(function(left, right) call("+", left, right), term_symbols)
  }
  formula_call <- if (is.null(outcome)) {
    call("~", rhs)
  } else {
    call("~", as.name(outcome), rhs)
  }
  stats::as.formula(formula_call, env = parent.frame())
}

ScidrModelMatrixInfo <- function(df_model, model_terms) {
  formula <- ScidrRegressionFormula(model_terms)
  x_full <- stats::model.matrix(formula, data = df_model)
  x <- x_full[, -1, drop = FALSE]
  assign <- attr(x_full, "assign")[-1]
  term_labels <- model_terms
  list(
    Formula = formula, X = x, Assign = assign, TermLabels = term_labels,
    ColumnTerms = term_labels[assign]
  )
}

ScidrFitOrdinaryMultinomial <- function(df_model, outcome, model_terms, predictor_vars, Standardize, reference) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' is required for ordinary multinomial regression.")
  formula <- ScidrRegressionFormula(model_terms, outcome)
  fit_data <- if (Standardize) ScidrScaleContinuousColumns(df_model, outcome = outcome) else df_model
  std_data <- ScidrScaleContinuousColumns(df_model, outcome = outcome)
  warnings_vec <- character(0)
  fit_model <- function(data) withCallingHandlers(
    nnet::multinom(formula, data = data, trace = FALSE, Hess = TRUE),
    warning = function(w) { warnings_vec <<- c(warnings_vec, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  model <- fit_model(fit_data)
  std_model <- fit_model(std_data)
  model_summary <- summary(model)
  std_summary <- summary(std_model)
  coefs <- model_summary$coefficients
  ses <- model_summary$standard.errors
  std_coefs <- std_summary$coefficients
  classes <- rownames(coefs)
  rows <- list()
  for (class in classes) {
    for (term in model_terms) {
      coef_name <- ScidrMatchingCoefficientName(colnames(coefs), term)
      estimate <- if (!is.na(coef_name)) coefs[class, coef_name] else NA_real_
      se <- if (!is.na(coef_name)) ses[class, coef_name] else NA_real_
      std_estimate <- if (!is.na(coef_name)) std_coefs[class, coef_name] else NA_real_
      p_value <- if (!is.na(se) && se > 0) 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE) else NA_real_
      rows[[paste(class, term, sep = "::")]] <- ScidrMulticategoryTermRow(
        term, estimate, std_estimate, se, p_value, class, reference,
        paste0(class, " vs ", reference), NA, converged = model$convergence == 0,
        warnings = paste(unique(warnings_vec), collapse = "; ")
      )
    }
  }
  probabilities <- stats::predict(model, type = "probs")
  null_model <- nnet::multinom(ScidrRegressionFormula(character(0), outcome), data = fit_data, trace = FALSE)
  lr <- 2 * (as.numeric(stats::logLik(model)) - as.numeric(stats::logLik(null_model)))
  lr_df <- attr(stats::logLik(model), "df") - attr(stats::logLik(null_model), "df")
  list(
    Model = model, StandardizedModel = std_model, TermTable = dplyr::bind_rows(rows),
    Probabilities = probabilities, Warnings = warnings_vec, Converged = model$convergence == 0,
    ModelStat = lr, ModelStatType = "LR chi-square",
    ModelPValue = stats::pchisq(lr, lr_df, lower.tail = FALSE),
    AIC = stats::AIC(model), BIC = stats::BIC(model), DevianceExplained = NA_real_, Tuning = list()
  )
}

ScidrFitOrdinaryOrdinal <- function(df_model, outcome, model_terms, predictor_vars, Standardize) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required for proportional-odds regression.")
  formula <- ScidrRegressionFormula(model_terms, outcome)
  fit_data <- if (Standardize) ScidrScaleContinuousColumns(df_model, outcome = outcome) else df_model
  std_data <- ScidrScaleContinuousColumns(df_model, outcome = outcome)
  warnings_vec <- character(0)
  fit_model <- function(data) withCallingHandlers(
    MASS::polr(formula, data = data, Hess = TRUE, method = "logistic"),
    warning = function(w) { warnings_vec <<- c(warnings_vec, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  model <- fit_model(fit_data)
  std_model <- fit_model(std_data)
  coef_table <- coef(summary(model))
  std_coefs <- stats::coef(std_model)
  rows <- lapply(model_terms, function(term) {
    coef_name <- ScidrMatchingCoefficientName(rownames(coef_table), term)
    estimate <- if (!is.na(coef_name)) coef_table[coef_name, "Value"] else NA_real_
    se <- if (!is.na(coef_name)) coef_table[coef_name, "Std. Error"] else NA_real_
    std_estimate <- if (!is.na(coef_name)) std_coefs[[coef_name]] else NA_real_
    p_value <- if (!is.na(se) && se > 0) 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE) else NA_real_
    ScidrMulticategoryTermRow(
      term, estimate, std_estimate, se, p_value, "Higher category", "Lower category",
      "higher vs lower", NA, converged = model$convergence == 0,
      warnings = paste(unique(warnings_vec), collapse = "; ")
    )
  })
  probabilities <- stats::predict(model, type = "probs")
  null_model <- MASS::polr(ScidrRegressionFormula(character(0), outcome), data = fit_data, Hess = TRUE)
  lr <- 2 * (as.numeric(stats::logLik(model)) - as.numeric(stats::logLik(null_model)))
  lr_df <- attr(stats::logLik(model), "df") - attr(stats::logLik(null_model), "df")
  list(
    Model = model, StandardizedModel = std_model, TermTable = dplyr::bind_rows(rows),
    Probabilities = probabilities, Warnings = warnings_vec, Converged = model$convergence == 0,
    Thresholds = model$zeta, ModelStat = lr, ModelStatType = "LR chi-square",
    ModelPValue = stats::pchisq(lr, lr_df, lower.tail = FALSE),
    AIC = stats::AIC(model), BIC = stats::BIC(model), DevianceExplained = NA_real_,
    Tuning = list(Thresholds = model$zeta, ProportionalOdds = TRUE)
  )
}

ScidrFitPenalizedMultinomial <- function(df_model, outcome, model_terms, predictor_vars,
                                         method, cv_folds, lambda_choice, seed, reference) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' is required for penalized multinomial regression.")
  mm <- ScidrModelMatrixInfo(df_model, model_terms)
  x <- mm$X
  x_std <- scale(x)
  x_std[, attr(x_std, "scaled:scale") == 0] <- 0
  covariates <- setdiff(model_terms, predictor_vars)
  penalty_factor <- ifelse(mm$ColumnTerms %in% covariates, 0, 1)
  names(penalty_factor) <- colnames(x)
  alpha_grid <- if (method == "ridge") 0 else if (method == "lasso") 1 else seq(0, 1, 0.1)
  set.seed(seed)
  foldid <- sample(rep(seq_len(cv_folds), length.out = nrow(df_model)))
  fits <- list()
  std_fits <- list()
  errors <- data.frame(Alpha = numeric(), LambdaMin = numeric(), Lambda1SE = numeric(), CVMMin = numeric())
  for (alpha in alpha_grid) {
    fit <- glmnet::cv.glmnet(x, df_model[[outcome]], family = "multinomial",
      type.multinomial = "grouped", alpha = alpha, foldid = foldid,
      standardize = TRUE, penalty.factor = penalty_factor)
    std_fit <- glmnet::cv.glmnet(x_std, df_model[[outcome]], family = "multinomial",
      type.multinomial = "grouped", alpha = alpha, foldid = foldid,
      lambda = fit$lambda, standardize = FALSE, penalty.factor = penalty_factor)
    fits[[as.character(alpha)]] <- fit
    std_fits[[as.character(alpha)]] <- std_fit
    errors <- rbind(errors, data.frame(
      Alpha = alpha, LambdaMin = fit$lambda.min, Lambda1SE = fit$lambda.1se,
      CVMMin = min(fit$cvm, na.rm = TRUE)
    ))
  }
  best_alpha <- errors$Alpha[which.min(errors$CVMMin)]
  fit <- fits[[as.character(best_alpha)]]
  std_fit <- std_fits[[as.character(best_alpha)]]
  selected_lambda <- if (lambda_choice == "lambda.min") fit$lambda.min else fit$lambda.1se
  coefs <- stats::coef(fit, s = selected_lambda)
  std_coefs <- stats::coef(std_fit, s = selected_lambda)
  classes <- names(coefs)
  rows <- list()
  for (class in setdiff(classes, reference)) {
    class_coef <- as.matrix(coefs[[class]])[, 1]
    ref_coef <- as.matrix(coefs[[reference]])[, 1]
    class_std <- as.matrix(std_coefs[[class]])[, 1]
    ref_std <- as.matrix(std_coefs[[reference]])[, 1]
    for (term in model_terms) {
      coef_name <- ScidrMatchingCoefficientName(names(class_coef), term)
      estimate <- if (!is.na(coef_name)) class_coef[[coef_name]] - ref_coef[[coef_name]] else 0
      std_estimate <- if (!is.na(coef_name)) class_std[[coef_name]] - ref_std[[coef_name]] else 0
      rows[[paste(class, term, sep = "::")]] <- ScidrMulticategoryTermRow(
        term, estimate, std_estimate, NA_real_, NA_real_, class, reference,
        paste0(class, " vs ", reference), abs(estimate) > 1e-6,
        selected_lambda, best_alpha
      )
    }
  }
  probabilities <- stats::predict(fit, newx = x, s = selected_lambda, type = "response")[, , 1]
  list(
    Model = fit, StandardizedModel = std_fit, TermTable = dplyr::bind_rows(rows),
    Probabilities = probabilities, Warnings = character(0), Converged = TRUE,
    ModelStat = NA_real_, ModelStatType = NA_character_, ModelPValue = NA_real_,
    AIC = NA_real_, BIC = NA_real_,
    DevianceExplained = tryCatch(max(fit$glmnet.fit$dev.ratio, na.rm = TRUE), error = function(e) NA_real_),
    Tuning = list(
      AlphaGrid = alpha_grid, CrossValidationErrors = errors, BestAlpha = best_alpha,
      BestLambda = selected_lambda, PenaltyFactors = penalty_factor,
      UnpenalizedTerms = covariates, TypeMultinomial = "grouped"
    )
  )
}

ScidrFitPenalizedOrdinal <- function(df_model, outcome, model_terms, predictor_vars,
                                     method, cv_folds, lambda_choice, seed) {
  if (!requireNamespace("ordinalNet", quietly = TRUE)) {
    stop(
      "Package 'ordinalNet' is required for penalized proportional-odds regression. ",
      "Install it or explicitly use outcome_modes = c(", outcome, " = 'multinomial') ",
      "to fit a nominal penalized model instead."
    )
  }
  mm <- ScidrModelMatrixInfo(df_model, model_terms)
  x <- mm$X
  x_std <- scale(x)
  x_std[, attr(x_std, "scaled:scale") == 0] <- 0
  covariates <- setdiff(model_terms, predictor_vars)
  penalty_factor <- ifelse(mm$ColumnTerms %in% covariates, 0, 1)
  names(penalty_factor) <- colnames(x)
  alpha_grid <- if (method == "ridge") 0 else if (method == "lasso") 1 else seq(0, 1, 0.1)
  set.seed(seed)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(df_model)))
  folds <- split(seq_len(nrow(df_model)), fold_id)
  tuning <- list()
  tuning_rows <- list()
  for (alpha in alpha_grid) {
    tune <- ordinalNet::ordinalNetTune(
      x = x, y = df_model[[outcome]], folds = folds, nFolds = cv_folds,
      printProgress = FALSE, warn = FALSE, alpha = alpha,
      family = "cumulative", link = "logit", parallelTerms = TRUE,
      nonparallelTerms = FALSE, standardize = TRUE,
      penaltyFactors = penalty_factor
    )
    mean_loglik <- rowMeans(tune$loglik, na.rm = TRUE)
    se_loglik <- apply(tune$loglik, 1, stats::sd, na.rm = TRUE) / sqrt(cv_folds)
    tuning[[as.character(alpha)]] <- tune
    tuning_rows[[as.character(alpha)]] <- data.frame(
      Alpha = alpha, LambdaIndex = seq_along(tune$lambdaVals), Lambda = tune$lambdaVals,
      MeanLogLikelihood = mean_loglik, SELogLikelihood = se_loglik,
      MeanMisclassification = rowMeans(tune$misclass, na.rm = TRUE),
      MeanDevianceExplained = rowMeans(tune$devPct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  cv_table <- dplyr::bind_rows(tuning_rows)
  best_row <- which.max(cv_table$MeanLogLikelihood)
  best_alpha <- cv_table$Alpha[[best_row]]
  alpha_table <- cv_table[cv_table$Alpha == best_alpha, , drop = FALSE]
  best_within_alpha <- which.max(alpha_table$MeanLogLikelihood)
  if (lambda_choice == "lambda.1se") {
    cutoff <- alpha_table$MeanLogLikelihood[[best_within_alpha]] - alpha_table$SELogLikelihood[[best_within_alpha]]
    eligible <- which(alpha_table$MeanLogLikelihood >= cutoff)
    selected_row <- eligible[[which.max(alpha_table$Lambda[eligible])]]
  } else {
    selected_row <- best_within_alpha
  }
  selected_lambda <- alpha_table$Lambda[[selected_row]]
  selected_index <- alpha_table$LambdaIndex[[selected_row]]
  lambda_values <- tuning[[as.character(best_alpha)]]$lambdaVals
  model <- ordinalNet::ordinalNet(
    x = x, y = df_model[[outcome]], alpha = best_alpha, standardize = TRUE,
    penaltyFactors = penalty_factor, family = "cumulative", link = "logit",
    parallelTerms = TRUE, nonparallelTerms = FALSE, lambdaVals = lambda_values,
    warn = FALSE, keepTrainingData = TRUE
  )
  std_model <- ordinalNet::ordinalNet(
    x = x_std, y = df_model[[outcome]], alpha = best_alpha, standardize = FALSE,
    penaltyFactors = penalty_factor, family = "cumulative", link = "logit",
    parallelTerms = TRUE, nonparallelTerms = FALSE, lambdaVals = lambda_values,
    warn = FALSE, keepTrainingData = TRUE
  )
  coef_values <- stats::coef(model, whichLambda = selected_index)
  std_values <- stats::coef(std_model, whichLambda = selected_index)
  rows <- lapply(model_terms, function(term) {
    coef_name <- ScidrMatchingCoefficientName(names(coef_values), term)
    estimate <- if (!is.na(coef_name)) unname(coef_values[[coef_name]]) else 0
    std_estimate <- if (!is.na(coef_name) && coef_name %in% names(std_values)) unname(std_values[[coef_name]]) else 0
    ScidrMulticategoryTermRow(
      term, estimate, std_estimate, NA_real_, NA_real_, "Higher category", "Lower category",
      "higher vs lower", abs(estimate) > 1e-6, selected_lambda, best_alpha
    )
  })
  probabilities <- stats::predict(model, newx = x, whichLambda = selected_index, type = "response")
  list(
    Model = model, StandardizedModel = std_model, TermTable = dplyr::bind_rows(rows),
    Probabilities = probabilities, Warnings = character(0), Converged = TRUE,
    ModelStat = NA_real_, ModelStatType = NA_character_, ModelPValue = NA_real_,
    AIC = NA_real_, BIC = NA_real_,
    DevianceExplained = alpha_table$MeanDevianceExplained[[selected_row]],
    Tuning = list(
      AlphaGrid = alpha_grid, CrossValidationErrors = cv_table, BestAlpha = best_alpha,
      BestLambda = selected_lambda, BestLambdaIndex = selected_index,
      PenaltyFactors = penalty_factor, UnpenalizedTerms = covariates,
      Family = "cumulative", Link = "logit", ParallelTerms = TRUE
    )
  )
}

ScidrMulticategoryPredictions <- function(probabilities, observed, outcome, mode,
                                          full_row_count, complete_rows) {
  probabilities <- as.matrix(probabilities)
  class_levels <- colnames(probabilities)
  if (is.null(class_levels)) class_levels <- levels(observed)
  predicted_index <- max.col(probabilities, ties.method = "first")
  predicted_labels <- class_levels[predicted_index]
  model_rows <- which(complete_rows)
  rows <- lapply(seq_along(class_levels), function(class_index) {
    out <- data.frame(
      Outcome = rep(outcome, full_row_count),
      OutcomeLevel = class_levels[[class_index]],
      ReferenceLevel = if (mode == "ordinal") "Lower category" else levels(observed)[[1]],
      Contrast = if (mode == "ordinal") "higher vs lower" else paste0(class_levels[[class_index]], " probability"),
      Observed = rep(NA_real_, full_row_count), Predicted = rep(NA_real_, full_row_count),
      PredictedProbability = rep(NA_real_, full_row_count), PredictedClass = rep(NA_integer_, full_row_count),
      ObservedClass = rep(NA_character_, full_row_count), PredictedClassLabel = rep(NA_character_, full_row_count),
      ClassificationThreshold = rep(NA_real_, full_row_count), Residual = rep(NA_real_, full_row_count),
      StandardizedResidual = rep(NA_real_, full_row_count), Leverage = rep(NA_real_, full_row_count),
      CooksDistance = rep(NA_real_, full_row_count), stringsAsFactors = FALSE
    )
    observed_index <- as.integer(observed)
    out$Observed[model_rows] <- as.integer(observed_index == class_index)
    out$Predicted[model_rows] <- as.integer(predicted_index == class_index)
    out$PredictedProbability[model_rows] <- probabilities[, class_index]
    out$PredictedClass[model_rows] <- predicted_index
    out$ObservedClass[model_rows] <- as.character(observed)
    out$PredictedClassLabel[model_rows] <- predicted_labels
    out$Residual[model_rows] <- out$Observed[model_rows] - probabilities[, class_index]
    out
  })
  dplyr::bind_rows(rows)
}

ScidrMulticategoryDiagnostics <- function(fit, observed, predictions, outcome, mode, method) {
  class_levels <- levels(observed)
  predicted_labels <- predictions[
    predictions$OutcomeLevel == class_levels[[1]],
    c("ObservedClass", "PredictedClassLabel"),
    drop = FALSE
  ]
  predicted_labels <- predicted_labels[!is.na(predicted_labels$ObservedClass), , drop = FALSE]
  confusion <- table(
    factor(predicted_labels$ObservedClass, levels = class_levels),
    factor(predicted_labels$PredictedClassLabel, levels = class_levels),
    dnn = c("Observed", "Predicted")
  )
  accuracy <- sum(diag(confusion)) / sum(confusion)
  recalls <- diag(confusion) / rowSums(confusion)
  balanced_accuracy <- mean(recalls, na.rm = TRUE)
  auc <- NA_real_
  if (requireNamespace("pROC", quietly = TRUE) && nlevels(observed) > 2) {
    probs <- fit$Probabilities
    auc <- tryCatch(as.numeric(pROC::auc(pROC::multiclass.roc(observed, probs, quiet = TRUE))), error = function(e) NA_real_)
  }
  comparison_levels <- if (mode == "ordinal") "Higher category" else setdiff(class_levels, class_levels[[1]])
  reference_levels <- if (mode == "ordinal") "Lower category" else rep(class_levels[[1]], length(comparison_levels))
  contrasts <- if (mode == "ordinal") "higher vs lower" else paste0(comparison_levels, " vs ", reference_levels)
  data.frame(
    Outcome = outcome, OutcomeLevel = comparison_levels, ReferenceLevel = reference_levels,
    Contrast = contrasts, R2 = NA_real_, AdjustedR2 = NA_real_, RMSE = NA_real_, ResidualSD = NA_real_,
    AIC = fit$AIC, BIC = fit$BIC, AUC = auc, AUCLowerCI = NA_real_, AUCUpperCI = NA_real_,
    OptimalThreshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
    PositivePredictiveValue = NA_real_, NegativePredictiveValue = NA_real_,
    BalancedAccuracy = balanced_accuracy, Accuracy = accuracy,
    ConfusionMatrix = I(rep(list(confusion), length(comparison_levels))), McFaddenR2 = NA_real_,
    CrossValidationError = NA_real_, SelectedLambda = ifelse(method == "lm", NA_real_, fit$Tuning$BestLambda),
    SelectedAlpha = ifelse(method == "lm", NA_real_, fit$Tuning$BestAlpha),
    RetainedPredictorCount = sum(fit$TermTable$Selected %in% TRUE, na.rm = TRUE),
    RetainedPredictors = I(rep(list(unique(fit$TermTable$Predictor[fit$TermTable$Selected %in% TRUE])), length(comparison_levels))),
    DevianceExplained = fit$DevianceExplained, SampleSize = length(observed),
    PredictorCount = length(unique(fit$TermTable$Predictor)), Converged = fit$Converged,
    SeparationDetected = FALSE, AliasedTermCount = sum(fit$TermTable$Aliased %in% TRUE, na.rm = TRUE),
    AliasedTerms = I(rep(list(unique(fit$TermTable$Predictor[fit$TermTable$Aliased %in% TRUE])), length(comparison_levels))),
    ModelPValue = fit$ModelPValue, ModelStat = fit$ModelStat, ModelStatType = fit$ModelStatType,
    Warnings = paste(unique(fit$Warnings), collapse = "; "), stringsAsFactors = FALSE
  )
}

ScidrMulticategoryModelingAdvice <- function() {
  c(
    nominal = "Use multinomial regression for unordered categories when comparisons share one coherent outcome model.",
    ordinal = "Proportional-odds models use ordering efficiently, but the proportional-odds assumption should be assessed scientifically and diagnostically.",
    one_vs_rest = "One-vs-rest models answer level-versus-all-others questions but produce overlapping comparisons and require multiplicity-aware interpretation.",
    binary_subset = "Binary subsets change the analysis population and estimand; excluded levels and retained sample sizes should be reported.",
    skip = "Skip is appropriate when a multi-category outcome is outside the analysis question; skipped outcomes are recorded in metadata."
  )
}

ScidrAnnotateRegressionPiece <- function(piece, outcome, mode, outcome_level, reference_level, contrast, outcome_label) {
  comparison_label <- if (mode %in% c("linear", "binary")) outcome_label else paste0(outcome_label, ": ", contrast)
  annotate <- function(tbl) {
    if (is.null(tbl) || nrow(tbl) == 0) return(tbl)
    tbl$Outcome <- outcome
    tbl$OutcomeLabel <- outcome_label
    tbl$OutcomeLevel <- outcome_level
    tbl$ReferenceLevel <- reference_level
    tbl$Contrast <- contrast
    tbl$ComparisonLabel <- comparison_label
    tbl$ComparisonIndex <- 1L
    tbl$OutcomeMode <- mode
    tbl
  }
  piece$LargeTable <- annotate(piece$LargeTable)
  piece$RegressionMatrix <- annotate(piece$RegressionMatrix)
  piece$VariableImportanceMatrix <- annotate(piece$VariableImportanceMatrix)
  piece$Predictions <- annotate(piece$Predictions)
  piece$Diagnostics <- annotate(piece$Diagnostics)
  piece$ModelSummary <- annotate(piece$ModelSummary)
  piece
}

ScidrCombineRegressionPieces <- function(pieces, FDR, FDRAlpha) {
  large_table <- dplyr::bind_rows(lapply(pieces, `[[`, "LargeTable"))
  if (!"ComparisonLabel" %in% names(large_table)) {
    large_table <- ScidrAnnotateTableDefaults(large_table)
  }
  comparison_levels <- unique(large_table$ComparisonLabel)
  outcome_levels <- unique(large_table$Outcome)
  large_table$OutcomeIndex <- match(large_table$Outcome, outcome_levels)
  large_table$ComparisonIndex <- match(large_table$ComparisonLabel, comparison_levels)
  if (FDR) {
    adjust_rows <- !is.na(large_table$PValue)
    large_table$FDR <- NA_real_
    large_table$FDR[adjust_rows] <- stats::p.adjust(large_table$PValue[adjust_rows], method = "fdr")
  } else {
    large_table$FDR <- NA_real_
  }
  large_table$Stars <- ScidrPValueStars(if (FDR) large_table$FDR else large_table$PValue)
  large_table$HoverText <- ScidrRegressionHoverText(large_table)
  reporting_rows <- large_table$Predictor %in% unique(unlist(lapply(pieces, function(x) x$Metadata$AnalysisSettings$PredictorVars)))
  regression_cols <- c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel", "ReferenceLevel",
    "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel", "PredictorLabel",
    "OutcomeFamily", "Estimate", "Effect", "EffectType", "StandardizedBeta", "PValue", "FDR",
    "Stars", "Selected", "VariableImportance", "VariableImportanceType", "HoverText"
  )
  importance_cols <- c(
    "OutcomeIndex", "ComparisonIndex", "PredictorIndex", "Outcome", "OutcomeLevel", "ReferenceLevel",
    "Contrast", "ComparisonLabel", "OutcomeMode", "Predictor", "OutcomeLabel", "PredictorLabel",
    "OutcomeFamily", "VariableImportance", "VariableImportanceType", "Selected", "HoverText"
  )
  regression_matrix <- large_table[reporting_rows, regression_cols, drop = FALSE]
  importance_matrix <- large_table[reporting_rows, importance_cols, drop = FALSE]
  predictions <- dplyr::bind_rows(lapply(pieces, `[[`, "Predictions"))
  diagnostics <- dplyr::bind_rows(lapply(pieces, `[[`, "Diagnostics"))
  model_summary <- dplyr::bind_rows(lapply(pieces, `[[`, "ModelSummary"))
  if (!"ComparisonLabel" %in% names(model_summary)) model_summary <- ScidrAnnotateTableDefaults(model_summary)
  model_summary$OutcomeIndex <- match(model_summary$Outcome, outcome_levels)
  model_summary$ComparisonIndex <- match(model_summary$ComparisonLabel, comparison_levels)
  plots <- ScidrRegressionPlots(regression_matrix, importance_matrix, predictions, model_summary)
  metadata_first <- pieces[[1]]$Metadata
  metadata_first$AnalysisSettings$FDR <- FDR
  metadata_first$AnalysisSettings$FDRAlpha <- FDRAlpha
  list(
    Models = lapply(pieces, function(piece) {
      models <- piece$Models
      if (length(models) == 1) return(models[[1]])
      lapply(models, function(model) if (is.list(model) && length(model) == 1) model[[1]] else model)
    }),
    FormattedTable = ScidrMultivariableGtTable(large_table[reporting_rows, , drop = FALSE], formatted = TRUE),
    LargeTable = as.data.frame(large_table),
    RegressionMatrix = as.data.frame(regression_matrix),
    VariableImportanceMatrix = as.data.frame(importance_matrix),
    Predictions = as.data.frame(predictions),
    Diagnostics = as.data.frame(diagnostics),
    ModelSummary = as.data.frame(model_summary),
    Multicollinearity = pieces[[1]]$Multicollinearity,
    Plots = plots,
    Metadata = metadata_first
  )
}

ScidrAnnotateTableDefaults <- function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl)
  if (!"OutcomeLevel" %in% names(tbl)) tbl$OutcomeLevel <- NA_character_
  if (!"ReferenceLevel" %in% names(tbl)) tbl$ReferenceLevel <- NA_character_
  if (!"Contrast" %in% names(tbl)) tbl$Contrast <- ifelse(tbl$OutcomeFamily == "linear", "Continuous", "Event vs reference")
  if (!"ComparisonLabel" %in% names(tbl)) tbl$ComparisonLabel <- tbl$OutcomeLabel
  if (!"ComparisonIndex" %in% names(tbl)) tbl$ComparisonIndex <- match(tbl$ComparisonLabel, unique(tbl$ComparisonLabel))
  if (!"OutcomeMode" %in% names(tbl)) tbl$OutcomeMode <- tbl$OutcomeFamily
  tbl
}

ScidrMultivariableGtTable <- function(reporting, formatted = TRUE) {
  # Report-facing table mirroring ScidrUnivariateGtTable(): the primary
  # estimate is the effect (odds ratio for logistic, beta for linear) on the
  # same scale as its confidence interval, combined into a single cell, with
  # significant p-values shown in bold. Penalized fits (lasso/ridge/elasticnet)
  # carry no p-values or CIs, so those cells are left blank and nothing is
  # bolded; the estimate cell then shows the penalized coefficient alone.
  if (nrow(reporting) == 0) {
    return(gt::gt(reporting))
  }

  # Significance follows the same source as the plotted stars: FDR when it was
  # calculated (ordinary models with FDR = TRUE), otherwise the raw p-value.
  sig_source <- ifelse(!is.na(reporting$FDR), reporting$FDR, reporting$PValue)

  table_data <- reporting %>%
    dplyr::mutate(
      Significant = !is.na(sig_source) & sig_source < 0.05,
      Estimate_CI = dplyr::case_when(
        is.na(.data$Effect) ~ NA_character_,
        is.na(.data$LowerCI) | is.na(.data$UpperCI) ~
          formatC(.data$Effect, digits = 3, format = "fg"),
        TRUE ~ paste0(
          formatC(.data$Effect, digits = 3, format = "fg"),
          " (",
          formatC(.data$LowerCI, digits = 3, format = "fg"),
          ", ",
          formatC(.data$UpperCI, digits = 3, format = "fg"),
          ")"
        )
      ),
      P = dplyr::case_when(
        is.na(.data$PValue) ~ NA_character_,
        .data$PValue < 0.001 ~ "<0.001",
        TRUE ~ formatC(.data$PValue, digits = 2, format = "fg")
      )
    )

  if (formatted) {
    group_column <- if ("ComparisonLabel" %in% names(table_data)) "ComparisonLabel" else "OutcomeLabel"
    table_data <- table_data %>%
      dplyr::select(
        dplyr::all_of(group_column),
        PredictorLabel,
        EffectType,
        SampleSize,
        Estimate_CI,
        P,
        Significant
      )

    out <- gt::gt(table_data, groupname_col = group_column) %>%
      gt::cols_label(
        PredictorLabel = "Variable",
        EffectType = "Effect",
        SampleSize = "N",
        Estimate_CI = "Estimate (95% CI)",
        P = "p-value"
      ) %>%
      gt::cols_hide(columns = "Significant") %>%
      gt::sub_missing(columns = c("Estimate_CI", "P"), missing_text = "") %>%
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_body(
          columns = "P",
          rows = Significant
        )
      )
  } else {
    table_data <- table_data %>%
      dplyr::select(
        Outcome,
        OutcomeLabel,
        PredictorLabel,
        EffectType,
        SampleSize,
        Effect,
        StandardizedBeta,
        LowerCI,
        UpperCI,
        dplyr::all_of("PValue"),
        FDR
      )

    out <- gt::gt(table_data, groupname_col = "OutcomeLabel") %>%
      gt::cols_label(
        Outcome = "Outcome",
        PredictorLabel = "Variable",
        EffectType = "Effect",
        SampleSize = "N",
        Effect = "Estimate",
        StandardizedBeta = "Std. beta",
        LowerCI = "95% CI Low",
        UpperCI = "95% CI High",
        PValue = "p-value",
        FDR = "FDR"
      ) %>%
      gt::fmt_number(
        columns = c("Effect", "StandardizedBeta", "LowerCI", "UpperCI"),
        decimals = 3
      ) %>%
      gt::fmt_number(columns = c("PValue", "FDR"), decimals = 3)
  }

  out
}

ScidrOutcomeFamily <- function(x, outcome) {
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
      return(if (is.ordered(non_missing)) "ordinal" else "multinomial")
    }
  }
  stop("Outcome ", outcome, " must be numeric, logical, or a two-level factor.")
}

ScidrPredictorMissingnessSummary <- function(Data, model_terms) {
  if (length(model_terms) == 0) {
    return(data.frame(
      Variable = character(0),
      MissingCount = integer(0),
      MissingProportion = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Variable = model_terms,
    MissingCount = vapply(model_terms, function(var) sum(is.na(Data[[var]])), integer(1)),
    MissingProportion = vapply(model_terms, function(var) mean(is.na(Data[[var]])), numeric(1)),
    stringsAsFactors = FALSE
  )
}

ScidrPrepareRegressionModelData <- function(Data,
                                            outcome,
                                            model_terms,
                                            retained_model_terms,
                                            dropped_terms,
                                            missing_data_strategy,
                                            impute_method,
                                            max_missing_predictor) {
  use_imputation <- missing_data_strategy %in% c("drop_sparse_impute", "impute")
  model_terms_outcome <- retained_model_terms
  model_data <- Data[, unique(c(outcome, model_terms_outcome)), drop = FALSE]
  outcome_nonmissing_rows <- !is.na(model_data[[outcome]])
  outcome_nonmissing_n <- sum(outcome_nonmissing_rows)
  complete_case_n <- sum(stats::complete.cases(model_data))

  if (use_imputation) {
    model_data <- model_data[outcome_nonmissing_rows, , drop = FALSE]
    imputed_variables <- character(0)
    for (term in model_terms_outcome) {
      missing_idx <- is.na(model_data[[term]])
      if (any(missing_idx)) {
        model_data[[term]][missing_idx] <- ScidrImputeValue(model_data[[term]], impute_method)
        imputed_variables <- c(imputed_variables, term)
      }
    }
    model_rows <- rep(FALSE, nrow(Data))
    model_rows[which(outcome_nonmissing_rows)] <- TRUE
  } else {
    complete_rows <- stats::complete.cases(model_data)
    model_data <- model_data[complete_rows, , drop = FALSE]
    model_rows <- rep(FALSE, nrow(Data))
    model_rows[which(complete_rows)] <- TRUE
    imputed_variables <- character(0)
  }

  list(
    ModelData = model_data,
    ModelRows = model_rows,
    ModelTerms = model_terms_outcome,
    DroppedVariables = dropped_terms,
    ImputedVariables = unique(imputed_variables),
    MissingRemoved = nrow(Data) - nrow(model_data),
    PercentRemoved = (nrow(Data) - nrow(model_data)) / nrow(Data),
    OriginalN = nrow(Data),
    OutcomeNonMissingN = outcome_nonmissing_n,
    CompleteCaseN = complete_case_n,
    MaxMissingPredictor = max_missing_predictor
  )
}

ScidrImputeValue <- function(x, impute_method) {
  if (all(is.na(x))) {
    stop("Cannot impute a predictor or covariate with all values missing.")
  }

  if (is.numeric(x)) {
    return(stats::median(x, na.rm = TRUE))
  }

  observed <- x[!is.na(x)]
  mode_value <- names(sort(table(observed), decreasing = TRUE))[1]

  if (is.factor(x)) {
    return(factor(mode_value, levels = levels(x)))
  }
  if (is.logical(x)) {
    return(mode_value == "TRUE")
  }
  mode_value
}

ScidrNotEnoughRowsMessage <- function(outcome, missing_info, min_rows_required) {
  paste0(
    "Outcome ", outcome, " does not have enough rows after missing-data handling. ",
    "Original N = ", missing_info$OriginalN, "; ",
    "outcome non-missing N = ", missing_info$OutcomeNonMissingN, "; ",
    "complete-case N before imputation = ", missing_info$CompleteCaseN, "; ",
    "final model N = ", nrow(missing_info$ModelData), "; ",
    "minimum required N = ", min_rows_required, ". ",
    "Dropped for missingness: ",
    ifelse(length(missing_info$DroppedVariables) == 0, "none", paste(missing_info$DroppedVariables, collapse = ", ")),
    ". Try MissingDataStrategy = 'drop_sparse_impute', increasing MaxMissingPredictor, reducing predictors, or lowering MinCompleteCases if statistically appropriate."
  )
}

ScidrPrepareBinaryOutcome <- function(x, outcome) {
  if (is.logical(x)) {
    return(list(Value = as.integer(x), Levels = c("FALSE", "TRUE")))
  }
  x <- as.factor(x)
  if (length(levels(x)) != 2) {
    stop("Outcome ", outcome, " must have exactly two levels for logistic regression.")
  }
  list(Value = as.integer(x == levels(x)[2]), Levels = levels(x))
}

ScidrRegressionLabels <- function(Data, vars, Relabel) {
  labels <- stats::setNames(vars, vars)
  if (Relabel) {
    attr_labels <- sjlabelled::get_label(Data[, vars, drop = FALSE], def.value = vars)
    attr_labels <- as.character(attr_labels)
    names(attr_labels) <- vars
    missing_labels <- is.na(attr_labels) | attr_labels == ""
    attr_labels[missing_labels] <- vars[missing_labels]
    labels <- attr_labels
  }
  as.list(labels)
}

ScidrScaleContinuousColumns <- function(df, outcome = NULL, standardize_outcome = FALSE) {
  for (nm in names(df)) {
    if (!standardize_outcome && identical(nm, outcome)) {
      next
    }
    if (is.numeric(df[[nm]]) && length(unique(df[[nm]][!is.na(df[[nm]])])) > 2) {
      df[[nm]] <- as.numeric(scale(df[[nm]]))
    }
  }
  df
}

ScidrFitOrdinaryRegression <- function(df_model,
                                       outcome,
                                       model_terms,
                                       predictor_vars,
                                       outcome_family,
                                       standardize_estimates) {
  formula <- ScidrRegressionFormula(model_terms, outcome)
  fit_data <- if (standardize_estimates) {
    ScidrScaleContinuousColumns(df_model, outcome = outcome, standardize_outcome = outcome_family == "linear")
  } else {
    df_model
  }

  warnings_vec <- character(0)
  model <- withCallingHandlers(
    if (outcome_family == "linear") {
      stats::lm(formula, data = fit_data)
    } else {
      stats::glm(formula, data = fit_data, family = stats::binomial())
    },
    warning = function(w) {
      warnings_vec <<- c(warnings_vec, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  standardized_data <- ScidrScaleContinuousColumns(
    df_model,
    outcome = outcome,
    standardize_outcome = outcome_family == "linear"
  )
  standardized_model <- withCallingHandlers(
    if (outcome_family == "linear") {
      stats::lm(formula, data = standardized_data)
    } else {
      stats::glm(formula, data = standardized_data, family = stats::binomial())
    },
    warning = function(w) {
      warnings_vec <<- c(warnings_vec, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  coef_table <- as.data.frame(summary(model)$coefficients)
  coef_table$Term <- rownames(coef_table)
  names(coef_table)[1:min(4, ncol(coef_table))] <- c("Estimate", "StandardError", "TestStatistic", "PValue")[1:min(4, ncol(coef_table))]

  std_coef <- stats::coef(standardized_model)
  ci <- tryCatch(stats::confint(model), error = function(e) NULL)
  base_model <- if (outcome_family == "linear") {
    stats::lm(ScidrRegressionFormula(character(0), outcome), data = fit_data)
  } else {
    stats::glm(ScidrRegressionFormula(character(0), outcome), data = fit_data, family = stats::binomial())
  }

  # Rank deficiency: lm()/glm() set aliased (perfectly collinear) coefficients
  # to NA and drop them from summary(), so they would otherwise appear as
  # unexplained blank cells. Surface them explicitly.
  full_coef <- stats::coef(model)
  aliased_terms <- names(full_coef)[is.na(full_coef)]
  if (length(aliased_terms) > 0) {
    warnings_vec <- c(
      warnings_vec,
      paste0(
        "Perfectly collinear (aliased) terms dropped by the model: ",
        paste(aliased_terms, collapse = ", "),
        ". Check for redundant predictors (e.g. change scores alongside their components)."
      )
    )
  }

  # Separation / non-convergence: for logistic fits, quasi-complete separation
  # makes coefficients and standard errors explode (fitted probabilities pinned
  # at 0/1). The estimates are not interpretable, so flag the fit and blank its
  # numeric outputs rather than let them poison downstream tables and plots.
  separation_detected <- FALSE
  if (outcome_family == "logistic") {
    fitted_probs <- stats::fitted(model)
    extreme_fitted <- any(fitted_probs < 1e-5 | fitted_probs > 1 - 1e-5, na.rm = TRUE)
    std_slopes <- std_coef[names(std_coef) != "(Intercept)"]
    huge_coef <- any(abs(std_slopes) > 10, na.rm = TRUE)
    warn_sep <- any(grepl(
      "fitted probabilities numerically 0 or 1|did not converge",
      warnings_vec,
      ignore.case = TRUE
    ))
    separation_detected <- extreme_fitted || huge_coef || warn_sep || !isTRUE(model$converged)
    if (separation_detected) {
      warnings_vec <- c(
        warnings_vec,
        paste0(
          "Separation / non-convergence detected for outcome '", outcome,
          "': logistic coefficients are unreliable and have been blanked. ",
          "Consider a penalized method (lasso/ridge) or fewer predictors."
        )
      )
    }
  }
  converged <- if (outcome_family == "logistic") {
    isTRUE(model$converged) && !separation_detected
  } else {
    TRUE
  }

  rows <- lapply(model_terms, function(term) {
    row_name <- ScidrMatchingCoefficientName(coef_table$Term, term)
    is_aliased <- any(startsWith(aliased_terms, term))
    estimate <- if (!is.na(row_name)) coef_table$Estimate[coef_table$Term == row_name][1] else NA_real_
    p_value <- if (!is.na(row_name)) coef_table$PValue[coef_table$Term == row_name][1] else NA_real_
    se <- if (!is.na(row_name)) coef_table$StandardError[coef_table$Term == row_name][1] else NA_real_
    test_stat <- if (!is.na(row_name)) coef_table$TestStatistic[coef_table$Term == row_name][1] else NA_real_
    lower_ci <- if (!is.null(ci) && row_name %in% rownames(ci)) ci[row_name, 1] else NA_real_
    upper_ci <- if (!is.null(ci) && row_name %in% rownames(ci)) ci[row_name, 2] else NA_real_
    std_beta <- if (!is.na(row_name) && row_name %in% names(std_coef)) unname(std_coef[[row_name]]) else NA_real_

    # A separated logistic fit yields no trustworthy estimates for any term, so
    # blank the whole model's inferential outputs rather than plot noise (and
    # skip the reduced-model refit that would only add more separation warnings).
    if (separation_detected) {
      estimate <- NA_real_
      p_value <- NA_real_
      se <- NA_real_
      test_stat <- NA_real_
      lower_ci <- NA_real_
      upper_ci <- NA_real_
      std_beta <- NA_real_
      variable_importance <- NA_real_
    } else {
      variable_importance <- ScidrOrdinaryVariableImportance(model, base_model, term, outcome_family)
    }

    data.frame(
      Predictor = term,
      Estimate = estimate,
      RawEstimate = estimate,
      Effect = if (outcome_family == "logistic") exp(estimate) else estimate,
      EffectType = if (outcome_family == "logistic") "Odds Ratio" else "Beta",
      StandardizedBeta = std_beta,
      StandardError = se,
      TestStatistic = test_stat,
      DegreesFreedom = if (outcome_family == "linear") stats::df.residual(model) else NA_real_,
      LowerCI = if (outcome_family == "logistic") exp(lower_ci) else lower_ci,
      UpperCI = if (outcome_family == "logistic") exp(upper_ci) else upper_ci,
      PValue = p_value,
      FDR = NA_real_,
      VariableImportance = variable_importance,
      VariableImportanceType = if (outcome_family == "logistic") "Likelihood Ratio Chi-square" else "Partial R2",
      Selected = NA,
      Aliased = is_aliased,
      Lambda = NA_real_,
      Alpha = NA_real_,
      Converged = converged,
      Warnings = paste(unique(warnings_vec), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  list(
    Model = model,
    StandardizedModel = standardized_model,
    TermTable = dplyr::bind_rows(rows),
    Warnings = warnings_vec,
    SeparationDetected = separation_detected,
    AliasedTerms = aliased_terms,
    Converged = converged,
    Tuning = list()
  )
}

ScidrFitPenalizedRegression <- function(df_model,
                                        outcome,
                                        model_terms,
                                        predictor_vars,
                                        outcome_family,
                                        method,
                                        cv_folds,
                                        lambda_choice,
                                        seed) {
  model_formula <- ScidrRegressionFormula(model_terms)
  x_full <- stats::model.matrix(model_formula, data = df_model)
  column_terms <- attr(x_full, "assign")[-1]
  x_raw <- x_full[, -1, drop = FALSE]
  x_std <- scale(x_raw)
  x_std[, attr(x_std, "scaled:scale") == 0] <- 0
  x_std <- as.matrix(x_std)
  y <- df_model[[outcome]]
  family <- if (outcome_family == "logistic") "binomial" else "gaussian"
  alpha_grid <- if (method == "ridge") 0 else if (method == "lasso") 1 else seq(0, 1, 0.1)

  # Covariates are mandatory adjustments: exempt them from the penalty so
  # ridge/lasso/elasticnet never shrink or drop them. penalty.factor is
  # per model-matrix column; map columns back to their originating term.
  term_labels <- model_terms
  covariate_terms <- setdiff(model_terms, predictor_vars)
  penalty_factor <- ifelse(term_labels[column_terms] %in% covariate_terms, 0, 1)
  names(penalty_factor) <- colnames(x_raw)

  set.seed(seed)
  foldid <- sample(rep(seq_len(cv_folds), length.out = nrow(df_model)))

  warnings_vec <- character(0)
  cv_fits <- list()
  cv_std_fits <- list()
  cv_errors <- data.frame(Alpha = numeric(0), LambdaMin = numeric(0), Lambda1SE = numeric(0), CVMMin = numeric(0))
  for (alpha in alpha_grid) {
    fit <- withCallingHandlers(
      glmnet::cv.glmnet(
        x = x_raw,
        y = y,
        family = family,
        alpha = alpha,
        foldid = foldid,
        standardize = TRUE,
        penalty.factor = penalty_factor
      ),
      warning = function(w) {
        warnings_vec <<- c(warnings_vec, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    std_fit <- withCallingHandlers(
      glmnet::cv.glmnet(
        x = x_std,
        y = y,
        family = family,
        alpha = alpha,
        foldid = foldid,
        lambda = fit$lambda,
        standardize = FALSE,
        penalty.factor = penalty_factor
      ),
      warning = function(w) {
        warnings_vec <<- c(warnings_vec, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    cv_fits[[as.character(alpha)]] <- fit
    cv_std_fits[[as.character(alpha)]] <- std_fit
    cv_errors <- rbind(
      cv_errors,
      data.frame(
        Alpha = alpha,
        LambdaMin = fit$lambda.min,
        Lambda1SE = fit$lambda.1se,
        CVMMin = min(fit$cvm, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }

  best_alpha <- cv_errors$Alpha[which.min(cv_errors$CVMMin)]
  cv_fit <- cv_fits[[as.character(best_alpha)]]
  cv_std_fit <- cv_std_fits[[as.character(best_alpha)]]
  selected_lambda <- if (lambda_choice == "lambda.min") cv_fit$lambda.min else cv_fit$lambda.1se
  coef_mat <- as.matrix(stats::coef(cv_fit, s = selected_lambda))
  coef_values <- coef_mat[, 1]
  std_coef_mat <- as.matrix(stats::coef(cv_std_fit, s = selected_lambda))
  std_coef_values <- std_coef_mat[, 1]
  linear_predictor <- as.numeric(stats::predict(cv_fit, newx = x_raw, s = selected_lambda, type = "link"))
  predicted_response <- as.numeric(stats::predict(cv_fit, newx = x_raw, s = selected_lambda, type = "response"))
  retained <- names(coef_values)[names(coef_values) != "(Intercept)" & abs(coef_values) > 0]

  rows <- lapply(model_terms, function(term) {
    coef_names <- names(coef_values)
    row_name <- ScidrMatchingCoefficientName(coef_names, term)
    estimate <- if (!is.na(row_name)) unname(coef_values[[row_name]]) else 0
    std_estimate <- if (!is.na(row_name) && row_name %in% names(std_coef_values)) unname(std_coef_values[[row_name]]) else 0
    data.frame(
      Predictor = term,
      Estimate = estimate,
      RawEstimate = estimate,
      Effect = if (outcome_family == "logistic") exp(estimate) else estimate,
      EffectType = if (outcome_family == "logistic") "Odds Ratio" else "Beta",
      StandardizedBeta = std_estimate,
      StandardError = NA_real_,
      TestStatistic = NA_real_,
      DegreesFreedom = NA_real_,
      LowerCI = NA_real_,
      UpperCI = NA_real_,
      PValue = NA_real_,
      FDR = NA_real_,
      VariableImportance = abs(std_estimate),
      VariableImportanceType = "Absolute Standardized Beta",
      Selected = !is.na(row_name) && abs(estimate) > 0,
      Aliased = FALSE,
      Lambda = selected_lambda,
      Alpha = best_alpha,
      Converged = TRUE,
      Warnings = paste(unique(warnings_vec), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  list(
    Model = cv_fit,
    StandardizedModel = cv_std_fit,
    TermTable = dplyr::bind_rows(rows),
    Warnings = warnings_vec,
    SeparationDetected = FALSE,
    AliasedTerms = character(0),
    Converged = TRUE,
    Tuning = list(
      AlphaGrid = alpha_grid,
      CrossValidationErrors = cv_errors,
      BestAlpha = best_alpha,
      BestLambda = selected_lambda,
      PenaltyFactors = penalty_factor,
      UnpenalizedTerms = covariate_terms
    ),
    PenalizedPredictions = list(
      LinearPredictor = linear_predictor,
      PredictedResponse = predicted_response,
      SelectedLambda = selected_lambda,
      SelectedAlpha = best_alpha,
      RetainedPredictors = retained,
      CrossValidationError = min(cv_fit$cvm, na.rm = TRUE),
      DevianceExplained = tryCatch(max(cv_fit$glmnet.fit$dev.ratio, na.rm = TRUE), error = function(e) NA_real_)
    )
  )
}

ScidrMatchingCoefficientName <- function(coef_names, term) {
  quoted_term <- ScidrQuoteFormulaNames(term)
  candidates <- coef_names[
    coef_names %in% c(term, quoted_term) |
      startsWith(coef_names, term) |
      startsWith(coef_names, quoted_term)
  ]
  candidates <- setdiff(candidates, "(Intercept)")
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  candidates[[1]]
}

ScidrOrdinaryVariableImportance <- function(model, base_model, term, outcome_family) {
  full_loglik <- as.numeric(stats::logLik(model))
  reduced_formula <- stats::update.formula(
    stats::formula(model),
    paste(". ~ . -", ScidrQuoteFormulaNames(term))
  )
  reduced_model <- tryCatch(
    if (outcome_family == "linear") stats::lm(reduced_formula, data = model$model) else stats::glm(reduced_formula, data = model$model, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(reduced_model)) {
    return(NA_real_)
  }
  if (outcome_family == "linear") {
    rss_full <- sum(stats::residuals(model)^2)
    rss_reduced <- sum(stats::residuals(reduced_model)^2)
    return((rss_reduced - rss_full) / rss_reduced)
  }
  2 * (full_loglik - as.numeric(stats::logLik(reduced_model)))
}

ScidrRegressionPredictions <- function(fit_result,
                                       df_model,
                                       outcome,
                                       outcome_family,
                                       full_row_count,
                                       complete_rows) {
  out <- data.frame(
    Outcome = rep(outcome, full_row_count),
    Observed = rep(NA_real_, full_row_count),
    Predicted = rep(NA_real_, full_row_count),
    PredictedProbability = rep(NA_real_, full_row_count),
    PredictedClass = rep(NA, full_row_count),
    ClassificationThreshold = rep(NA_real_, full_row_count),
    Residual = rep(NA_real_, full_row_count),
    StandardizedResidual = rep(NA_real_, full_row_count),
    Leverage = rep(NA_real_, full_row_count),
    CooksDistance = rep(NA_real_, full_row_count),
    .ModelComplete = complete_rows,
    stringsAsFactors = FALSE
  )

  if (inherits(fit_result$Model, "cv.glmnet")) {
    pred_response <- fit_result$PenalizedPredictions$PredictedResponse
    out$Observed[complete_rows] <- df_model[[outcome]]
    if (outcome_family == "linear") {
      out$Predicted[complete_rows] <- pred_response
      out$Residual[complete_rows] <- df_model[[outcome]] - pred_response
      out$StandardizedResidual[complete_rows] <- as.numeric(scale(out$Residual[complete_rows]))
    } else {
      threshold <- ScidrOptimalThreshold(df_model[[outcome]], pred_response)$Threshold
      out$PredictedProbability[complete_rows] <- pred_response
      out$PredictedClass[complete_rows] <- as.integer(pred_response >= threshold)
      out$Predicted[complete_rows] <- out$PredictedClass[complete_rows]
      out$ClassificationThreshold[complete_rows] <- threshold
    }
    return(out)
  }

  model <- fit_result$Model
  out$Observed[complete_rows] <- df_model[[outcome]]
  if (outcome_family == "linear") {
    fitted_values <- as.numeric(stats::predict(model, type = "response"))
    residuals <- df_model[[outcome]] - fitted_values
    out$Predicted[complete_rows] <- fitted_values
    out$Residual[complete_rows] <- residuals
    out$StandardizedResidual[complete_rows] <- as.numeric(stats::rstandard(model))
    out$Leverage[complete_rows] <- as.numeric(stats::hatvalues(model))
    out$CooksDistance[complete_rows] <- as.numeric(stats::cooks.distance(model))
  } else {
    probs <- as.numeric(stats::predict(model, type = "response"))
    threshold <- ScidrOptimalThreshold(df_model[[outcome]], probs)$Threshold
    pred_class <- as.integer(probs >= threshold)
    out$PredictedProbability[complete_rows] <- probs
    out$PredictedClass[complete_rows] <- pred_class
    out$Predicted[complete_rows] <- pred_class
    out$ClassificationThreshold[complete_rows] <- threshold
    out$Residual[complete_rows] <- df_model[[outcome]] - probs
    out$StandardizedResidual[complete_rows] <- tryCatch(as.numeric(stats::rstandard(model)), error = function(e) NA_real_)
    out$Leverage[complete_rows] <- tryCatch(as.numeric(stats::hatvalues(model)), error = function(e) NA_real_)
    out$CooksDistance[complete_rows] <- tryCatch(as.numeric(stats::cooks.distance(model)), error = function(e) NA_real_)
  }
  out
}

ScidrRegressionDiagnostics <- function(fit_result,
                                       df_model,
                                       outcome,
                                       outcome_family,
                                       method,
                                       predictions_complete) {
  base <- data.frame(
    R2 = NA_real_, AdjustedR2 = NA_real_, RMSE = NA_real_, ResidualSD = NA_real_,
    AIC = NA_real_, BIC = NA_real_, AUC = NA_real_, AUCLowerCI = NA_real_,
    AUCUpperCI = NA_real_, OptimalThreshold = NA_real_, Sensitivity = NA_real_,
    Specificity = NA_real_, PositivePredictiveValue = NA_real_,
    NegativePredictiveValue = NA_real_, BalancedAccuracy = NA_real_,
    Accuracy = NA_real_, ConfusionMatrix = I(list(matrix(NA, nrow = 2, ncol = 2))),
    McFaddenR2 = NA_real_, CrossValidationError = NA_real_,
    SelectedLambda = NA_real_, SelectedAlpha = NA_real_,
    RetainedPredictorCount = NA_integer_, RetainedPredictors = I(list(character(0))),
    DevianceExplained = NA_real_, SampleSize = nrow(df_model),
    PredictorCount = length(setdiff(names(df_model), outcome)),
    Converged = isTRUE(fit_result$Converged),
    SeparationDetected = isTRUE(fit_result$SeparationDetected),
    AliasedTermCount = length(fit_result$AliasedTerms),
    AliasedTerms = I(list(fit_result$AliasedTerms)),
    ModelPValue = NA_real_, ModelStat = NA_real_, ModelStatType = NA_character_,
    Warnings = paste(unique(fit_result$Warnings), collapse = "; "),
    stringsAsFactors = FALSE
  )

  if (inherits(fit_result$Model, "cv.glmnet")) {
    retained <- fit_result$PenalizedPredictions$RetainedPredictors
    base$CrossValidationError <- fit_result$PenalizedPredictions$CrossValidationError
    base$SelectedLambda <- fit_result$PenalizedPredictions$SelectedLambda
    base$SelectedAlpha <- fit_result$PenalizedPredictions$SelectedAlpha
    base$RetainedPredictorCount <- length(retained)
    base$RetainedPredictors <- I(list(retained))
    base$DevianceExplained <- fit_result$PenalizedPredictions$DevianceExplained
    # Penalized models have no valid classical omnibus test; ModelPValue stays NA
    # and cross-validated deviance explained is the column-level fit summary.
  } else if (outcome_family == "linear") {
    model_summary <- summary(fit_result$Model)
    residuals <- predictions_complete$Residual
    base$R2 <- unname(model_summary$r.squared)
    base$AdjustedR2 <- unname(model_summary$adj.r.squared)
    base$RMSE <- sqrt(mean(residuals^2, na.rm = TRUE))
    base$ResidualSD <- stats::sd(residuals, na.rm = TRUE)
    base$AIC <- stats::AIC(fit_result$Model)
    base$BIC <- stats::BIC(fit_result$Model)
    # Omnibus F-test for the whole model.
    fstat <- model_summary$fstatistic
    if (!is.null(fstat)) {
      base$ModelStat <- unname(fstat[1])
      base$ModelStatType <- "F"
      base$ModelPValue <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    }
  } else {
    probs <- predictions_complete$PredictedProbability
    observed <- predictions_complete$Observed
    threshold_info <- ScidrOptimalThreshold(observed, probs)
    threshold <- threshold_info$Threshold
    pred_class <- as.integer(probs >= threshold)
    metrics <- ScidrConfusionMetrics(observed, pred_class)
    null_model <- stats::glm(ScidrRegressionFormula(character(0), outcome), data = df_model, family = stats::binomial())
    base$AUC <- threshold_info$AUC
    base$AUCLowerCI <- threshold_info$AUCLowerCI
    base$AUCUpperCI <- threshold_info$AUCUpperCI
    base$OptimalThreshold <- threshold
    base$Sensitivity <- metrics$Sensitivity
    base$Specificity <- metrics$Specificity
    base$PositivePredictiveValue <- metrics$PositivePredictiveValue
    base$NegativePredictiveValue <- metrics$NegativePredictiveValue
    base$BalancedAccuracy <- metrics$BalancedAccuracy
    base$Accuracy <- metrics$Accuracy
    base$ConfusionMatrix <- I(list(metrics$ConfusionMatrix))
    base$McFaddenR2 <- 1 - as.numeric(stats::logLik(fit_result$Model)) / as.numeric(stats::logLik(null_model))
    base$AIC <- stats::AIC(fit_result$Model)
    base$BIC <- stats::BIC(fit_result$Model)
    base$Converged <- isTRUE(fit_result$Model$converged) && !isTRUE(fit_result$SeparationDetected)
    # Omnibus likelihood-ratio test of the full model against the null. Suppress
    # for separated fits, where the deviance drop is an artefact, not evidence.
    model_obj <- fit_result$Model
    lr_chisq <- model_obj$null.deviance - model_obj$deviance
    lr_df <- model_obj$df.null - model_obj$df.residual
    if (!isTRUE(fit_result$SeparationDetected) && is.finite(lr_chisq) && lr_df > 0) {
      base$ModelStat <- lr_chisq
      base$ModelStatType <- "LR chi-square"
      base$ModelPValue <- stats::pchisq(lr_chisq, lr_df, lower.tail = FALSE)
    }
  }
  base
}

ScidrOptimalThreshold <- function(observed, probability) {
  roc_obj <- pROC::roc(observed, probability, quiet = TRUE, direction = "<")
  coords <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) c(NA_real_, NA_real_, NA_real_))
  data.frame(
    Threshold = as.numeric(coords$threshold[[1]]),
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUCLowerCI = ci[[1]],
    AUCUpperCI = ci[[3]],
    stringsAsFactors = FALSE
  )
}

ScidrConfusionMetrics <- function(observed, predicted) {
  observed <- as.integer(observed)
  predicted <- as.integer(predicted)
  tp <- sum(observed == 1 & predicted == 1, na.rm = TRUE)
  tn <- sum(observed == 0 & predicted == 0, na.rm = TRUE)
  fp <- sum(observed == 0 & predicted == 1, na.rm = TRUE)
  fn <- sum(observed == 1 & predicted == 0, na.rm = TRUE)
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
  data.frame(
    Sensitivity = sensitivity,
    Specificity = specificity,
    PositivePredictiveValue = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    NegativePredictiveValue = ifelse(tn + fn > 0, tn / (tn + fn), NA_real_),
    BalancedAccuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    Accuracy = (tp + tn) / (tp + tn + fp + fn),
    ConfusionMatrix = I(list(matrix(c(tn, fp, fn, tp), nrow = 2, byrow = TRUE, dimnames = list(Observed = c("0", "1"), Predicted = c("0", "1"))))),
    stringsAsFactors = FALSE
  )
}

ScidrRegressionMulticollinearity <- function(Data, model_terms) {
  model_data <- Data[, model_terms, drop = FALSE]
  model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]
  x <- stats::model.matrix(ScidrRegressionFormula(model_terms), data = model_data)[, -1, drop = FALSE]
  empty <- list(
    CorrelationMatrix = matrix(numeric(0), nrow = 0, ncol = 0),
    CorrelationDataFrame = data.frame(),
    VIF = data.frame(Variable = character(0), VIF = numeric(0), stringsAsFactors = FALSE),
    Tolerance = data.frame(Variable = character(0), Tolerance = numeric(0), stringsAsFactors = FALSE),
    ConditionIndex = NA_real_,
    MaximumCorrelation = NA_real_,
    MaximumVIF = NA_real_,
    HighCorrelationPairs = data.frame(Variable1 = character(0), Variable2 = character(0), Correlation = numeric(0), stringsAsFactors = FALSE),
    HighVIFVariables = character(0),
    Recommendations = character(0)
  )
  if (ncol(x) == 0 || nrow(x) < 3) {
    return(empty)
  }
  if (ncol(x) == 1) {
    return(list(
      CorrelationMatrix = matrix(NA_real_, nrow = 1, ncol = 1, dimnames = list(colnames(x), colnames(x))),
      CorrelationDataFrame = data.frame(Variable1 = character(0), Variable2 = character(0), Correlation = numeric(0), stringsAsFactors = FALSE),
      VIF = data.frame(Variable = colnames(x), VIF = 1, stringsAsFactors = FALSE),
      Tolerance = data.frame(Variable = colnames(x), Tolerance = 1, stringsAsFactors = FALSE),
      ConditionIndex = 1,
      MaximumCorrelation = NA_real_,
      MaximumVIF = 1,
      HighCorrelationPairs = data.frame(Variable1 = character(0), Variable2 = character(0), Correlation = numeric(0), stringsAsFactors = FALSE),
      HighVIFVariables = character(0),
      Recommendations = character(0)
    ))
  }
  cor_mat <- stats::cor(x, use = "pairwise.complete.obs")
  diag(cor_mat) <- NA_real_
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("Variable1", "Variable2", "Correlation")
  cor_df <- cor_df[!is.na(cor_df$Correlation), , drop = FALSE]
  high_pairs <- cor_df[abs(cor_df$Correlation) > 0.80 & as.character(cor_df$Variable1) < as.character(cor_df$Variable2), , drop = FALSE]
  vif_values <- vapply(seq_len(ncol(x)), function(i) {
    other <- x[, -i, drop = FALSE]
    if (ncol(other) == 0) {
      return(1)
    }
    r2 <- summary(stats::lm(x[, i] ~ other))$r.squared
    1 / (1 - r2)
  }, numeric(1))
  names(vif_values) <- colnames(x)
  eigen_values <- eigen(stats::cor(x, use = "pairwise.complete.obs"), symmetric = TRUE, only.values = TRUE)$values
  condition_index <- sqrt(max(eigen_values, na.rm = TRUE) / min(eigen_values[eigen_values > 0], na.rm = TRUE))
  recommendations <- character(0)
  if (any(abs(cor_df$Correlation) > 0.80, na.rm = TRUE)) {
    recommendations <- c(recommendations, "Remove predictors", "Consider PCA")
  }
  if (any(vif_values > 5, na.rm = TRUE)) {
    recommendations <- c(recommendations, "Consider ridge regression", "Consider lasso", "Consider elastic net")
  }
  list(
    CorrelationMatrix = cor_mat,
    CorrelationDataFrame = cor_df,
    VIF = data.frame(Variable = names(vif_values), VIF = unname(vif_values), stringsAsFactors = FALSE),
    Tolerance = data.frame(Variable = names(vif_values), Tolerance = 1 / unname(vif_values), stringsAsFactors = FALSE),
    ConditionIndex = condition_index,
    MaximumCorrelation = max(abs(cor_df$Correlation), na.rm = TRUE),
    MaximumVIF = max(vif_values, na.rm = TRUE),
    HighCorrelationPairs = high_pairs,
    HighVIFVariables = names(vif_values)[vif_values > 5],
    Recommendations = unique(recommendations)
  )
}

ScidrPValueStars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

ScidrRegressionHoverText <- function(tbl) {
  importance_label <- ifelse(is.na(tbl$VariableImportanceType), "Variable Importance", tbl$VariableImportanceType)
  contrast_text <- if ("Contrast" %in% names(tbl)) paste0("<b>Contrast:</b> ", tbl$Contrast, "<br>") else ""
  reference_text <- if ("ReferenceLevel" %in% names(tbl)) paste0("<b>Reference:</b> ", tbl$ReferenceLevel, "<br>") else ""
  paste0(
    "<b>Outcome:</b> ", tbl$OutcomeLabel, "<br>",
    contrast_text,
    reference_text,
    "<b>Predictor:</b> ", tbl$PredictorLabel, "<br>",
    "<b>Estimate:</b> ", signif(tbl$Estimate, 3), "<br>",
    "<b>Effect:</b> ", signif(tbl$Effect, 3), " (", tbl$EffectType, ")<br>",
    "<b>Standardized Beta:</b> ", signif(tbl$StandardizedBeta, 3), "<br>",
    "<b>95% CI:</b> ", signif(tbl$LowerCI, 3), " to ", signif(tbl$UpperCI, 3), "<br>",
    "<b>P-value:</b> ", signif(tbl$PValue, 3), "<br>",
    "<b>FDR:</b> ", signif(tbl$FDR, 3), "<br>",
    "<b>", importance_label, ":</b> ", signif(tbl$VariableImportance, 3), "<br>",
    "<b>Sample size:</b> ", tbl$SampleSize
  )
}

ScidrRobustFillLimits <- function(values, probs = 0.99, floor = 0.1) {
  vals <- values[is.finite(values)]
  if (length(vals) == 0) {
    return(NULL)
  }
  m <- stats::quantile(abs(vals), probs = probs, na.rm = TRUE, names = FALSE)
  if (!is.finite(m) || m <= 0) {
    m <- max(abs(vals), na.rm = TRUE)
  }
  if (!is.finite(m) || m <= 0) {
    return(NULL)
  }
  m <- max(m, floor)
  c(-m, m)
}

ScidrColumnAnnotations <- function(model_summary) {
  empty <- data.frame(OutcomeLabel = character(0), Label = character(0), stringsAsFactors = FALSE)
  if (is.null(model_summary) || nrow(model_summary) == 0) {
    return(empty)
  }
  format_p <- function(p) {
    if (is.na(p)) {
      return(NA_character_)
    }
    if (p < 0.001) {
      return("p < 0.001")
    }
    paste0("p = ", formatC(p, digits = 2, format = "fg"))
  }
  labels <- vapply(seq_len(nrow(model_summary)), function(i) {
    row <- model_summary[i, , drop = FALSE]
    penalized <- !identical(as.character(row$RegressionMethod), "lm")
    if (penalized) {
      # No valid omnibus p-value for penalized fits: report CV deviance explained.
      dev <- row$DevianceExplained
      if (is.na(dev)) {
        return("")
      }
      return(paste0("Dev. expl. = ", formatC(dev, digits = 2, format = "f")))
    }
    if (isTRUE(row$SeparationDetected)) {
      return("unstable fit")
    }
    p <- row$ModelPValue
    if (is.na(p)) {
      return("")
    }
    format_p(p)
  }, character(1))
  data.frame(
    OutcomeLabel = if ("ComparisonLabel" %in% names(model_summary)) as.character(model_summary$ComparisonLabel) else as.character(model_summary$OutcomeLabel),
    Label = labels,
    stringsAsFactors = FALSE
  )
}

ScidrAddColumnAnnotations <- function(p, annotations) {
  annotations <- annotations[!is.na(annotations$Label) & annotations$Label != "", , drop = FALSE]
  if (nrow(annotations) == 0) {
    return(p)
  }
  p +
    ggplot2::geom_text(
      data = annotations,
      mapping = ggplot2::aes(x = .data$OutcomeLabel, y = Inf, label = .data$Label),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 3,
      na.rm = TRUE
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 22, r = 6, b = 6, l = 6))
}

ScidrRegressionPlots <- function(regression_matrix,
                                 variable_importance_matrix,
                                 predictions,
                                 model_summary) {
  plots <- list()
  column_annotations <- ScidrColumnAnnotations(model_summary)

  if (nrow(regression_matrix) > 0) {
    plot_data <- regression_matrix
    if ("ComparisonLabel" %in% names(plot_data)) plot_data$OutcomeLabel <- plot_data$ComparisonLabel
    plot_data$PredictorLabel <- factor(
      plot_data$PredictorLabel,
      levels = rev(unique(plot_data$PredictorLabel[order(plot_data$PredictorIndex)]))
    )
    plot_data$OutcomeLabel <- factor(
      plot_data$OutcomeLabel,
      levels = unique(plot_data$OutcomeLabel[order(plot_data$ComparisonIndex)])
    )

    # Robust symmetric limits so a single extreme (e.g. an unstable coefficient
    # that slipped through) can't dominate the diverging scale and wash out the
    # rest of the panel. Values outside the limits are clamped, not dropped.
    fill_limits <- ScidrRobustFillLimits(plot_data$StandardizedBeta)

    plots$RegressionMatrix <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$OutcomeLabel,
        y = .data$PredictorLabel,
        fill = .data$StandardizedBeta,
        text = .data$HoverText
      )
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        na.value = "grey90",
        limits = fill_limits,
        oob = function(x, range = fill_limits) {
          if (is.null(range)) {
            return(x)
          }
          x[x < range[1]] <- range[1]
          x[x > range[2]] <- range[2]
          x
        }
      ) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        fill = "Standardized beta"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    # Omnibus model summary at the top of each column: p-value for ordinary
    # models (so cells in a non-significant model are not over-interpreted),
    # or cross-validated deviance explained for penalized models.
    plots$RegressionMatrix <- ScidrAddColumnAnnotations(plots$RegressionMatrix, column_annotations)
  }

  if (nrow(variable_importance_matrix) > 0) {
    importance_data <- variable_importance_matrix
    if ("ComparisonLabel" %in% names(importance_data)) importance_data$OutcomeLabel <- importance_data$ComparisonLabel
    importance_data$PredictorLabel <- factor(
      importance_data$PredictorLabel,
      levels = rev(unique(importance_data$PredictorLabel[order(importance_data$PredictorIndex)]))
    )
    importance_data$OutcomeLabel <- factor(
      importance_data$OutcomeLabel,
      levels = unique(importance_data$OutcomeLabel[order(importance_data$ComparisonIndex)])
    )

    importance_upper <- ScidrRobustFillLimits(importance_data$VariableImportance)
    importance_limits <- if (is.null(importance_upper)) NULL else c(0, importance_upper[2])

    plots$VariableImportanceMatrix <- ggplot2::ggplot(
      importance_data,
      ggplot2::aes(
        x = .data$OutcomeLabel,
        y = .data$PredictorLabel,
        fill = .data$VariableImportance,
        text = .data$HoverText
      )
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradient(
        low = "grey95",
        high = "#1B7837",
        na.value = "grey90",
        limits = importance_limits,
        oob = function(x, range = importance_limits) {
          if (is.null(range)) {
            return(x)
          }
          x[x < range[1]] <- range[1]
          x[x > range[2]] <- range[2]
          x
        }
      ) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        fill = unique(stats::na.omit(importance_data$VariableImportanceType))[1]
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    plots$VariableImportanceMatrix <- ScidrAddColumnAnnotations(plots$VariableImportanceMatrix, column_annotations)
  }

  if (nrow(model_summary) > 0) {
    diagnostics_long <- ScidrRegressionDiagnosticsPlotData(model_summary)
    if (nrow(diagnostics_long) > 0) {
      plots$RegressionDiagnostics <- ggplot2::ggplot(
        diagnostics_long,
        ggplot2::aes(
          x = .data$OutcomeLabel,
          y = .data$Value,
          fill = .data$Metric
        )
      ) +
        ggplot2::geom_col(position = "dodge", width = 0.7) +
        ggplot2::facet_wrap(~Metric, scales = "free_y") +
        ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }
  }

  if (nrow(predictions) > 0) {
    if (any(!is.na(predictions$Predicted) & !is.na(predictions$Observed))) {
      plots$ObservedVsPredicted <- ggplot2::ggplot(
        predictions,
        ggplot2::aes(x = .data$Observed, y = .data$Predicted)
      ) +
        ggplot2::geom_point(alpha = 0.7, na.rm = TRUE) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::facet_wrap(~Outcome, scales = "free") +
        ggplot2::labs(x = "Observed", y = "Predicted") +
        ggplot2::theme_minimal()
    }

    if (any(!is.na(predictions$Residual))) {
      plots$Residuals <- ggplot2::ggplot(
        predictions,
        ggplot2::aes(x = .data$Predicted, y = .data$Residual)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::geom_point(alpha = 0.7, na.rm = TRUE) +
        ggplot2::facet_wrap(~Outcome, scales = "free") +
        ggplot2::labs(x = "Predicted", y = "Residual") +
        ggplot2::theme_minimal()
    }

    if (any(!is.na(predictions$PredictedProbability))) {
      plots$Calibration <- ggplot2::ggplot(
        predictions,
        ggplot2::aes(x = .data$PredictedProbability, y = .data$Observed)
      ) +
        ggplot2::geom_point(alpha = 0.35, na.rm = TRUE) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, na.rm = TRUE) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::facet_wrap(~Outcome) +
        ggplot2::labs(x = "Predicted probability", y = "Observed") +
        ggplot2::theme_minimal()
    }

    if (any(!is.na(predictions$Leverage)) || any(!is.na(predictions$CooksDistance))) {
      plots$Influence <- ggplot2::ggplot(
        predictions,
        ggplot2::aes(x = .data$Leverage, y = .data$CooksDistance)
      ) +
        ggplot2::geom_point(alpha = 0.7, na.rm = TRUE) +
        ggplot2::facet_wrap(~Outcome, scales = "free") +
        ggplot2::labs(x = "Leverage", y = "Cook's distance") +
        ggplot2::theme_minimal()
    }
  }

  plots
}

ScidrRegressionDiagnosticsPlotData <- function(model_summary) {
  metric_candidates <- c("R2", "AdjustedR2", "AUC", "McFaddenR2", "RMSE", "AIC", "BIC")
  metric_candidates <- intersect(metric_candidates, names(model_summary))
  if (length(metric_candidates) == 0) {
    return(data.frame(Outcome = character(0), OutcomeLabel = character(0), Metric = character(0), Value = numeric(0)))
  }

  rows <- lapply(metric_candidates, function(metric) {
    data.frame(
      Outcome = model_summary$Outcome,
      OutcomeLabel = if ("ComparisonLabel" %in% names(model_summary)) model_summary$ComparisonLabel else model_summary$OutcomeLabel,
      Metric = metric,
      Value = suppressWarnings(as.numeric(model_summary[[metric]])),
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  out <- out[!is.na(out$Value), , drop = FALSE]
  out
}
