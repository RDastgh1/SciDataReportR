#' Calculate diagnostic likelihood ratios
#'
#' Calculates diagnostic likelihood ratios for categorical diagnostic results
#' and binary outcomes. A likelihood ratio is the probability of a result among
#' outcome-positive participants divided by its probability among outcome-negative
#' participants; this is a diagnostic accuracy measure, not a nested-model test.
#'
#' Binary predictors produce sensitivity, specificity, LR+, and LR-. Predictors
#' with more than two levels return one likelihood ratio for each result level.
#' Numeric predictors with more than two observed values must be categorized first.
#'
#' @param data A data frame.
#' @param outcome_vars Character vector of binary outcome variable names.
#' @param predictor_vars Character vector of categorical or binary diagnostic predictor names.
#' @param positive_level Optional outcome-positive level specification: `NULL`, one
#'   level for all outcomes, or a named vector keyed by outcome variable.
#' @param predictor_positive_level Optional positive-result specification for binary
#'   predictors: `NULL`, one level for all predictors, or a named vector.
#' @param stratify_by Optional character vector of variables defining strata.
#' @param confidence_level Confidence level for log-scale LR confidence intervals.
#' @param continuity_correction Optional positive value added to all four calculation
#'   cells only if a zero cell occurs. `NULL` preserves zero and infinite estimates.
#' @param Relabel Logical; use attached variable labels when available.
#'
#' @return A list containing `Results`, `BinarySummary`, compact and expanded `gt`
#'   tables (`FormattedTable`, `LargeTable`, `BinaryFormattedTable`, and
#'   `BinaryLargeTable`), and `Metadata`.
#'
#' @seealso [PlotDiagnosticLRHeatmap()] to visualize `Results`.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' df_Labelled$DiagnosisBinary <- factor(df_Labelled$Diagnosis,
#'   levels = c("Control", "Impaired"))
#' lr <- DiagnosticLikelihoodRatioTable(df_Labelled, "DiagnosisBinary",
#'   c("sex", "Genotype"))
#' lr$FormattedTable
#' PlotDiagnosticLRHeatmap(lr)
#' @export
DiagnosticLikelihoodRatioTable <- function(data, outcome_vars, predictor_vars,
    positive_level = NULL, predictor_positive_level = NULL, stratify_by = NULL,
    confidence_level = 0.95, continuity_correction = NULL, Relabel = TRUE) {

  # Validate inputs

  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!is.character(outcome_vars) || !length(outcome_vars)) stop("outcome_vars must be a non-empty character vector.", call. = FALSE)
  if (!is.character(predictor_vars) || !length(predictor_vars)) stop("predictor_vars must be a non-empty character vector.", call. = FALSE)
  if (!is.null(stratify_by) && !is.character(stratify_by)) stop("stratify_by must be NULL or a character vector.", call. = FALSE)
  if (!is.numeric(confidence_level) || length(confidence_level) != 1 || is.na(confidence_level) || confidence_level <= 0 || confidence_level >= 1) stop("confidence_level must be one number strictly between 0 and 1.", call. = FALSE)
  if (!is.null(continuity_correction) && (!is.numeric(continuity_correction) || length(continuity_correction) != 1 || is.na(continuity_correction) || continuity_correction <= 0)) stop("continuity_correction must be NULL or one positive number.", call. = FALSE)
  if (!is.logical(Relabel) || length(Relabel) != 1 || is.na(Relabel)) stop("Relabel must be TRUE or FALSE.", call. = FALSE)
  all_vars <- unique(c(outcome_vars, predictor_vars, stratify_by))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars)) stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  if (!requireNamespace("gt", quietly = TRUE)) stop("Package 'gt' is required by DiagnosticLikelihoodRatioTable().", call. = FALSE)
  continuous_predictors <- predictor_vars[vapply(predictor_vars, function(variable) {
    values <- stats::na.omit(data[[variable]])
    is.numeric(values) && length(unique(values)) > 2
  }, logical(1))]
  if (length(continuous_predictors)) stop("DiagnosticLikelihoodRatioTable() does not automatically dichotomize continuous predictors. Categorize these variables before analysis: ", paste(continuous_predictors, collapse = ", "), call. = FALSE)

  # Prepare data

  display_labels <- ScidrDisplayLabels(data, all_vars, Relabel)
  outcome_levels <- lapply(outcome_vars, function(variable) ScidrDiagnosticResolveLevels(data[[variable]], variable, positive_level, "Outcome"))
  names(outcome_levels) <- outcome_vars
  predictor_levels <- lapply(predictor_vars, function(variable) ScidrDiagnosticObservedLevels(data[[variable]]))
  names(predictor_levels) <- predictor_vars
  binary_levels <- lapply(predictor_vars, function(variable) {
    if (length(predictor_levels[[variable]]) != 2) return(NULL)
    ScidrDiagnosticResolveLevels(data[[variable]], variable, predictor_positive_level, "Predictor")
  })
  names(binary_levels) <- predictor_vars
  for (outcome in outcome_vars) message(display_labels[[outcome]], ": '", outcome_levels[[outcome]]$positive, "' treated as outcome-positive.")
  if (is.null(stratify_by)) {
    strata <- tibble::tibble(Stratum = "All participants")
  } else {
    strata <- data %>% dplyr::select(dplyr::all_of(stratify_by)) %>%
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x))) %>% dplyr::distinct()
    strata$Stratum <- apply(as.data.frame(strata), 1, function(values) paste(paste0(unname(display_labels[stratify_by]), " = ", values), collapse = " | "))
  }

  # Calculate likelihood ratios

  results <- list()
  binary_summary <- list()
  warning_messages <- character()
  z_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  for (outcome in outcome_vars) for (predictor in predictor_vars) for (stratum_index in seq_len(nrow(strata))) {
    outcome_info <- outcome_levels[[outcome]]
    required_vars <- unique(c(outcome, predictor, stratify_by))
    df_Analysis <- data %>% dplyr::select(dplyr::all_of(required_vars))
    stratum_row <- strata[stratum_index, , drop = FALSE]
    if (!is.null(stratify_by)) for (strat_var in stratify_by) {
      df_Analysis <- df_Analysis[!is.na(df_Analysis[[strat_var]]) & as.character(df_Analysis[[strat_var]]) == as.character(stratum_row[[strat_var]]), , drop = FALSE]
    }
    df_Analysis <- df_Analysis %>% dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[predictor]])) %>%
      dplyr::mutate(.OutcomePositive = as.character(.data[[outcome]]) == outcome_info$positive, .Result = as.character(.data[[predictor]]))
    stratum_label <- stratum_row$Stratum[[1]]
    observed_levels <- ScidrDiagnosticObservedLevels(df_Analysis[[predictor]])
    if (length(observed_levels) < 2) {
      warning_messages <- c(warning_messages, paste0("Predictor '", predictor, "' has fewer than two observed levels within stratum '", stratum_label, "' and was skipped."))
      next
    }
    n_case <- sum(df_Analysis$.OutcomePositive)
    n_control <- sum(!df_Analysis$.OutcomePositive)
    if (!n_case || !n_control) {
      warning_messages <- c(warning_messages, paste0("Outcome '", outcome, "' does not contain both levels within stratum '", stratum_label, "' and was skipped."))
      next
    }
    comparison_results <- purrr::map_dfr(observed_levels, function(result_level) {
      lr <- ScidrDiagnosticLRFromCounts(sum(df_Analysis$.OutcomePositive & df_Analysis$.Result == result_level), sum(!df_Analysis$.OutcomePositive & df_Analysis$.Result == result_level), n_case, n_control, z_value, continuity_correction)
      out <- tibble::tibble(Outcome = outcome, OutcomeLabel = display_labels[[outcome]], OutcomePositiveLevel = outcome_info$positive, Predictor = predictor, PredictorLabel = display_labels[[predictor]], PredictorLevels = length(observed_levels), ResultLevel = result_level, Stratum = stratum_label, N = nrow(df_Analysis), NCase = n_case, NControl = n_control)
      if (!is.null(stratify_by)) out <- dplyr::bind_cols(out, stratum_row[stratify_by])
      dplyr::bind_cols(out, lr)
    })
    results[[length(results) + 1L]] <- comparison_results
    if (length(observed_levels) == 2) {
      predictor_info <- binary_levels[[predictor]]
      if (is.null(predictor_info) || !all(c(predictor_info$positive, predictor_info$negative) %in% observed_levels)) {
        warning_messages <- c(warning_messages, paste0("Binary predictor '", predictor, "' has a requested level absent within stratum '", stratum_label, "' and was skipped."))
      } else {
        positive_row <- comparison_results %>% dplyr::filter(.data$ResultLevel == predictor_info$positive)
        negative_row <- comparison_results %>% dplyr::filter(.data$ResultLevel == predictor_info$negative)
        out <- tibble::tibble(Outcome = outcome, OutcomeLabel = display_labels[[outcome]], OutcomePositiveLevel = outcome_info$positive, Predictor = predictor, PredictorLabel = display_labels[[predictor]], PositiveTestLevel = predictor_info$positive, NegativeTestLevel = predictor_info$negative, Stratum = stratum_label, N = nrow(df_Analysis), NCase = n_case, NControl = n_control, Sensitivity = positive_row$CaseProbability, Specificity = negative_row$ControlProbability, LRPositive = positive_row$LikelihoodRatio, LRPositiveLowerCI = positive_row$LRLowerCI, LRPositiveUpperCI = positive_row$LRUpperCI, LRNegative = negative_row$LikelihoodRatio, LRNegativeLowerCI = negative_row$LRLowerCI, LRNegativeUpperCI = negative_row$LRUpperCI, ZeroCell = positive_row$ZeroCell | negative_row$ZeroCell, Corrected = positive_row$Corrected | negative_row$Corrected)
        if (!is.null(stratify_by)) out <- dplyr::bind_cols(out, stratum_row[stratify_by])
        binary_summary[[length(binary_summary) + 1L]] <- out
      }
    }
  }
  if (length(warning_messages)) warning(paste(unique(warning_messages), collapse = "\n"), call. = FALSE)
  if (!length(results)) stop("No diagnostic likelihood ratios could be calculated from the supplied data.", call. = FALSE)
  Results <- purrr::list_rbind(results)
  BinarySummary <- if (length(binary_summary)) purrr::list_rbind(binary_summary) else tibble::tibble()

  # Build outputs

  FormattedTable <- ScidrDiagnosticLRGtTable(Results, formatted = TRUE)
  LargeTable <- ScidrDiagnosticLRGtTable(Results, formatted = FALSE)
  BinaryFormattedTable <- if (nrow(BinarySummary)) ScidrDiagnosticBinaryGtTable(BinarySummary, formatted = TRUE) else NULL
  BinaryLargeTable <- if (nrow(BinarySummary)) ScidrDiagnosticBinaryGtTable(BinarySummary, formatted = FALSE) else NULL
  Metadata <- list(OutcomeLevels = outcome_levels, PredictorLevels = predictor_levels, BinaryPredictorLevels = binary_levels, OutcomeVariables = outcome_vars, PredictorVariables = predictor_vars, StratifyBy = stratify_by, ConfidenceLevel = confidence_level, ContinuityCorrection = continuity_correction, Relabel = Relabel)

  # Return result

  list(FormattedTable = FormattedTable, LargeTable = LargeTable, BinaryFormattedTable = BinaryFormattedTable, BinaryLargeTable = BinaryLargeTable, Results = Results, BinarySummary = BinarySummary, Metadata = Metadata)
}

# Centralize the numerically sensitive 2 x 2 calculation so level-specific and
# binary summaries cannot diverge.
ScidrDiagnosticLRFromCounts <- function(case_with_result, control_with_result,
    n_case, n_control, z_value, continuity_correction = NULL) {
  counts <- c(case_with_result, control_with_result, n_case - case_with_result, n_control - control_with_result)
  zero_cell <- any(counts == 0)
  corrected <- zero_cell && !is.null(continuity_correction)
  calculation_counts <- if (corrected) counts + continuity_correction else counts
  a <- calculation_counts[[1]]; b <- calculation_counts[[2]]; c <- calculation_counts[[3]]; d <- calculation_counts[[4]]
  case_probability <- a / (a + c)
  control_probability <- b / (b + d)
  likelihood_ratio <- case_probability / control_probability
  log_lr_se <- sqrt((1 / a) - (1 / (a + c)) + (1 / b) - (1 / (b + d)))
  ci_available <- is.finite(likelihood_ratio) && likelihood_ratio > 0 && is.finite(log_lr_se)
  tibble::tibble(CaseWithResult = counts[[1]], ControlWithResult = counts[[2]], CaseWithoutResult = counts[[3]], ControlWithoutResult = counts[[4]], CaseProbability = case_probability, ControlProbability = control_probability, LikelihoodRatio = likelihood_ratio, LRLowerCI = if (ci_available) exp(log(likelihood_ratio) - z_value * log_lr_se) else NA_real_, LRUpperCI = if (ci_available) exp(log(likelihood_ratio) + z_value * log_lr_se) else NA_real_, Log2LR = if (is.na(likelihood_ratio)) NA_real_ else if (likelihood_ratio == 0) -Inf else log2(likelihood_ratio), ZeroCell = zero_cell, Corrected = corrected)
}

ScidrDiagnosticObservedLevels <- function(x) {
  values <- stats::na.omit(x)
  if (is.factor(x)) return(as.character(levels(droplevels(x))))
  if (is.logical(x)) return(as.character(c(FALSE, TRUE)[c(FALSE, TRUE) %in% unique(values)]))
  if (is.numeric(x)) return(as.character(sort(unique(values))))
  sort(unique(as.character(values)))
}

ScidrDiagnosticResolveLevels <- function(x, variable, specification, role) {
  observed <- ScidrDiagnosticObservedLevels(x)
  if (length(observed) != 2) stop(role, " variable '", variable, "' must have exactly two observed levels. Observed levels: ", paste(observed, collapse = ", "), call. = FALSE)
  requested <- if (!is.null(specification) && !is.null(names(specification)) && variable %in% names(specification)) as.character(specification[[variable]]) else if (!is.null(specification) && length(specification) == 1) as.character(specification) else observed[[2]]
  if (!requested %in% observed) stop("Requested positive level '", requested, "' was not found in ", role, " variable '", variable, "'. Observed levels: ", paste(observed, collapse = ", "), call. = FALSE)
  list(positive = requested, negative = setdiff(observed, requested), levels = observed)
}

ScidrDiagnosticFormatLR <- function(estimate, lower, upper) {
  if (is.na(estimate)) return(NA_character_)
  if (!is.finite(estimate)) return(ifelse(estimate > 0, "Inf", "0"))
  if (is.na(lower) || is.na(upper)) return(formatC(estimate, digits = 3, format = "fg"))
  paste0(formatC(estimate, digits = 3, format = "fg"), " (", formatC(lower, digits = 3, format = "fg"), ", ", formatC(upper, digits = 3, format = "fg"), ")")
}

ScidrDiagnosticLRGtTable <- function(results, formatted = TRUE) {
  table_data <- results %>% dplyr::mutate(`Diagnostic result` = paste0(.data$PredictorLabel, ": ", .data$ResultLevel), `Case %` = 100 * .data$CaseProbability, `Control %` = 100 * .data$ControlProbability, `Likelihood ratio (95% CI)` = mapply(ScidrDiagnosticFormatLR, .data$LikelihoodRatio, .data$LRLowerCI, .data$LRUpperCI))
  columns <- if (formatted) c("Stratum", "OutcomeLabel", "Diagnostic result", "N", "Case %", "Control %", "Likelihood ratio (95% CI)") else c("Stratum", "Outcome", "OutcomeLabel", "OutcomePositiveLevel", "Predictor", "PredictorLabel", "ResultLevel", "N", "NCase", "NControl", "CaseWithResult", "ControlWithResult", "CaseWithoutResult", "ControlWithoutResult", "CaseProbability", "ControlProbability", "LikelihoodRatio", "LRLowerCI", "LRUpperCI", "Log2LR", "ZeroCell", "Corrected")
  table_data <- table_data %>% dplyr::select(dplyr::all_of(columns))
  table <- gt::gt(table_data, groupname_col = if (dplyr::n_distinct(table_data$Stratum) > 1) "Stratum" else NULL) %>% gt::tab_header(title = if (formatted) "Diagnostic likelihood ratios" else "Diagnostic likelihood ratios: expanded results") %>% gt::opt_row_striping()
  if (formatted) table %>% gt::fmt_number(columns = c(`Case %`, `Control %`), decimals = 1) %>% gt::cols_label(OutcomeLabel = "Outcome", N = "N") else table %>% gt::fmt_number(columns = c(CaseProbability, ControlProbability, LikelihoodRatio, LRLowerCI, LRUpperCI, Log2LR), decimals = 3)
}

ScidrDiagnosticBinaryGtTable <- function(results, formatted = TRUE) {
  table_data <- results %>% dplyr::mutate(Sensitivity = 100 * .data$Sensitivity, Specificity = 100 * .data$Specificity, `LR+ (95% CI)` = mapply(ScidrDiagnosticFormatLR, .data$LRPositive, .data$LRPositiveLowerCI, .data$LRPositiveUpperCI), `LR- (95% CI)` = mapply(ScidrDiagnosticFormatLR, .data$LRNegative, .data$LRNegativeLowerCI, .data$LRNegativeUpperCI))
  columns <- if (formatted) c("Stratum", "OutcomeLabel", "PredictorLabel", "PositiveTestLevel", "NegativeTestLevel", "N", "Sensitivity", "Specificity", "LR+ (95% CI)", "LR- (95% CI)") else names(results)
  table_data <- table_data %>% dplyr::select(dplyr::all_of(columns))
  table <- gt::gt(table_data, groupname_col = if (dplyr::n_distinct(table_data$Stratum) > 1) "Stratum" else NULL) %>% gt::tab_header(title = if (formatted) "Binary diagnostic test performance" else "Binary diagnostic test performance: expanded results") %>% gt::opt_row_striping()
  if (formatted) table %>% gt::fmt_number(columns = c(Sensitivity, Specificity), decimals = 1) %>% gt::cols_label(OutcomeLabel = "Outcome", PredictorLabel = "Predictor", PositiveTestLevel = "Positive result", NegativeTestLevel = "Negative result", N = "N", Sensitivity = "Sensitivity (%)", Specificity = "Specificity (%)") else table
}
