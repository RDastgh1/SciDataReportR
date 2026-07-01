#' Multivariable regression table
#'
#' Fit one multivariable regression model per outcome and return a stable,
#' label-aware regression object for downstream tables, diagnostics, and plots.
#'
#' @param Data Data frame containing outcomes, predictors, and covariates.
#' @param OutcomeVars Character vector of outcome variable names.
#' @param PredictorVars Character vector of predictor variable names.
#' @param Covars Optional character vector of covariate variable names.
#' @param Standardize Logical. If `TRUE`, ordinary models are fit on
#'   standardized continuous variables for the primary estimate. Standardized
#'   coefficients are always calculated separately regardless of this setting.
#' @param Relabel Logical. If `TRUE`, use variable labels from `sjlabelled`
#'   when available.
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
#'
#' @return A named list with stable components: `Models`, `FormattedTable`,
#'   `LargeTable`, `RegressionMatrix`, `VariableImportanceMatrix`,
#'   `Predictions`, `Diagnostics`, `ModelSummary`, `Multicollinearity`,
#'   `Plots`, and `Metadata`.
#' @export
MultivariableRegressionTable <- function(
    Data,
    OutcomeVars,
    PredictorVars,
    Covars = NULL,
    Standardize = TRUE,
    Relabel = TRUE,
    FDR = TRUE,
    FDRAlpha = 0.05,
    Method = c("lm", "ridge", "lasso", "elasticnet"),
    CVFolds = 10,
    Lambda = c("lambda.min", "lambda.1se"),
    Seed = 123,
    MissingDataStrategy = c("drop_sparse_impute", "impute", "complete_cases", "drop_sparse_complete_cases"),
    MaxMissingPredictor = 0.30,
    ImputeMethod = c("median_mode"),
    MinCompleteCases = NULL
) {

  Method <- match.arg(Method)
  Lambda <- match.arg(Lambda)
  MissingDataStrategy <- match.arg(MissingDataStrategy)
  ImputeMethod <- match.arg(ImputeMethod)

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

  if (any(outcome_families == "logistic") && !requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required for logistic regression diagnostics.")
  }

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
  formatted_table <- large_table[reporting_rows, c(
    "Outcome", "Predictor", "OutcomeFamily", "Estimate", "Effect", "EffectType",
    "StandardizedBeta", "LowerCI", "UpperCI", "PValue", "FDR",
    "VariableImportance", "VariableImportanceType", "SampleSize", "Selected",
    "Lambda", "Alpha"
  ), drop = FALSE]
  regression_matrix <- large_table[reporting_rows, c(
    "OutcomeIndex", "PredictorIndex", "Outcome", "Predictor", "OutcomeLabel",
    "PredictorLabel", "OutcomeFamily", "Estimate", "Effect", "EffectType",
    "StandardizedBeta", "PValue", "FDR", "Stars", "Selected",
    "VariableImportance", "VariableImportanceType", "HoverText"
  ), drop = FALSE]
  variable_importance_matrix <- large_table[reporting_rows, c(
    "OutcomeIndex", "PredictorIndex", "Outcome", "Predictor", "OutcomeLabel",
    "PredictorLabel", "OutcomeFamily", "VariableImportance",
    "VariableImportanceType", "Selected", "HoverText"
  ), drop = FALSE]

  predictions <- dplyr::bind_rows(prediction_tables)
  predictions$.ModelComplete <- NULL

  diagnostics <- dplyr::bind_rows(diagnostics_tables)
  model_summary <- diagnostics[, c(
    "Outcome", "OutcomeLabel", "OutcomeFamily", "RegressionMethod", "SampleSize",
    "MissingRemoved", "PercentRemoved", "PredictorCount", "Converged",
    "R2", "AdjustedR2", "AUC", "McFaddenR2", "RMSE", "AIC", "BIC",
    "DroppedPredictorCount", "ImputedPredictorCount"
  ), drop = FALSE]
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

  return(list(
    Models = model_list,
    FormattedTable = as.data.frame(formatted_table),
    LargeTable = as.data.frame(large_table),
    RegressionMatrix = as.data.frame(regression_matrix),
    VariableImportanceMatrix = as.data.frame(variable_importance_matrix),
    Predictions = as.data.frame(predictions),
    Diagnostics = as.data.frame(diagnostics),
    ModelSummary = as.data.frame(model_summary),
    Multicollinearity = multicollinearity,
    Plots = list(),
    Metadata = metadata
  ))
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
      stop(
        "Outcome ", outcome, " has ", n_levels, " levels. ",
        "Multinomial regression is not yet supported; logistic regression requires a two-level factor or logical outcome."
      )
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
  formula <- stats::reformulate(model_terms, response = outcome)
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
    warning = function(w) invokeRestart("muffleWarning")
  )

  coef_table <- as.data.frame(summary(model)$coefficients)
  coef_table$Term <- rownames(coef_table)
  names(coef_table)[1:min(4, ncol(coef_table))] <- c("Estimate", "StandardError", "TestStatistic", "PValue")[1:min(4, ncol(coef_table))]

  std_coef <- stats::coef(standardized_model)
  ci <- tryCatch(stats::confint(model), error = function(e) NULL)
  base_model <- if (outcome_family == "linear") {
    stats::lm(stats::reformulate("1", response = outcome), data = fit_data)
  } else {
    stats::glm(stats::reformulate("1", response = outcome), data = fit_data, family = stats::binomial())
  }

  rows <- lapply(model_terms, function(term) {
    row_name <- ScidrMatchingCoefficientName(coef_table$Term, term)
    estimate <- if (!is.na(row_name)) coef_table$Estimate[coef_table$Term == row_name][1] else NA_real_
    p_value <- if (!is.na(row_name)) coef_table$PValue[coef_table$Term == row_name][1] else NA_real_
    se <- if (!is.na(row_name)) coef_table$StandardError[coef_table$Term == row_name][1] else NA_real_
    test_stat <- if (!is.na(row_name)) coef_table$TestStatistic[coef_table$Term == row_name][1] else NA_real_
    lower_ci <- if (!is.null(ci) && row_name %in% rownames(ci)) ci[row_name, 1] else NA_real_
    upper_ci <- if (!is.null(ci) && row_name %in% rownames(ci)) ci[row_name, 2] else NA_real_
    std_beta <- if (!is.na(row_name) && row_name %in% names(std_coef)) unname(std_coef[[row_name]]) else NA_real_
    variable_importance <- ScidrOrdinaryVariableImportance(model, base_model, term, outcome_family)
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
      Lambda = NA_real_,
      Alpha = NA_real_,
      Converged = if (outcome_family == "logistic") isTRUE(model$converged) else TRUE,
      Warnings = paste(unique(warnings_vec), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  list(
    Model = model,
    StandardizedModel = standardized_model,
    TermTable = dplyr::bind_rows(rows),
    Warnings = warnings_vec,
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
  x_raw <- stats::model.matrix(stats::reformulate(model_terms), data = df_model)[, -1, drop = FALSE]
  x_std <- scale(x_raw)
  x_std[, attr(x_std, "scaled:scale") == 0] <- 0
  x_std <- as.matrix(x_std)
  y <- df_model[[outcome]]
  family <- if (outcome_family == "logistic") "binomial" else "gaussian"
  alpha_grid <- if (method == "ridge") 0 else if (method == "lasso") 1 else seq(0, 1, 0.1)

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
        standardize = TRUE
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
        standardize = FALSE
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
    Tuning = list(
      AlphaGrid = alpha_grid,
      CrossValidationErrors = cv_errors,
      BestAlpha = best_alpha,
      BestLambda = selected_lambda
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
  candidates <- coef_names[coef_names == term | startsWith(coef_names, paste0(term))]
  candidates <- setdiff(candidates, "(Intercept)")
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  candidates[[1]]
}

ScidrOrdinaryVariableImportance <- function(model, base_model, term, outcome_family) {
  full_loglik <- as.numeric(stats::logLik(model))
  reduced_formula <- stats::update.formula(stats::formula(model), paste(". ~ . -", term))
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
    Converged = TRUE, Warnings = paste(unique(fit_result$Warnings), collapse = "; "),
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
  } else if (outcome_family == "linear") {
    model_summary <- summary(fit_result$Model)
    residuals <- predictions_complete$Residual
    base$R2 <- unname(model_summary$r.squared)
    base$AdjustedR2 <- unname(model_summary$adj.r.squared)
    base$RMSE <- sqrt(mean(residuals^2, na.rm = TRUE))
    base$ResidualSD <- stats::sd(residuals, na.rm = TRUE)
    base$AIC <- stats::AIC(fit_result$Model)
    base$BIC <- stats::BIC(fit_result$Model)
  } else {
    probs <- predictions_complete$PredictedProbability
    observed <- predictions_complete$Observed
    threshold_info <- ScidrOptimalThreshold(observed, probs)
    threshold <- threshold_info$Threshold
    pred_class <- as.integer(probs >= threshold)
    metrics <- ScidrConfusionMetrics(observed, pred_class)
    null_model <- stats::glm(stats::reformulate("1", response = outcome), data = df_model, family = stats::binomial())
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
    base$Converged <- isTRUE(fit_result$Model$converged)
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
  x <- stats::model.matrix(stats::reformulate(model_terms), data = model_data)[, -1, drop = FALSE]
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
  paste0(
    "<b>Outcome:</b> ", tbl$OutcomeLabel, "<br>",
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
