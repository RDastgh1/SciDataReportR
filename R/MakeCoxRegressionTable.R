#' Cox Proportional Hazards Regression Table
#'
#' Fits one Cox proportional hazards model per event-predictor pair and returns
#' report-ready tables plus a tidy results dataframe compatible with
#' [PlotForestFromTable()].
#'
#' @details
#' Each model has the form `Surv(time_var, event) ~ predictor + covariates`.
#' Without `design`, models are fit with [survival::coxph()]. When `design` is
#' a `survey.design`, the corresponding design variables are used and models
#' are fit with [survey::svycoxph()]. `N` and `Events` are unweighted analytic
#' row counts in both cases.
#'
#' Numeric events must be coded 0/1. Logical events use `TRUE` as the event.
#' For two-level factors, the second factor level is the event. Character event
#' variables are deliberately rejected so the event ordering is never implicit.
#'
#' @param data Data frame containing the time, event, predictor, and covariate
#'   variables. When `design` is supplied, `data` provides the variable-label
#'   contract and the design variables provide the modeling data.
#' @param time_var Character scalar naming the follow-up-time variable.
#' @param event_vars Character vector naming event indicators.
#' @param predictor_vars Character vector naming predictors.
#' @param covariates Optional character vector naming adjustment variables.
#' @param design Optional `survey.design` whose variables contain all model
#'   variables. Models are fit with [survey::svycoxph()] when supplied.
#' @param reference_levels Optional named list or named character vector giving
#'   the reference level for categorical predictors. Unspecified predictors
#'   retain their first factor level.
#' @param Standardize Logical. If `TRUE`, standardize numeric predictors and
#'   covariates within each model's complete-case analysis set.
#' @param FDR Logical. If `TRUE`, add `FDR`, adjusted across all returned rows
#'   using [ApplyFDRCorrection()].
#' @param FDRAlpha Numeric threshold used for FDR-adjusted significance.
#' @param CheckPH Logical. If `TRUE`, calculate the global [survival::cox.zph()]
#'   p-value for ordinary Cox models. Survey models return `NA`.
#' @param ReturnModels Logical. If `TRUE`, retain fitted models in
#'   `ModelSummaries`.
#' @param Relabel Logical. If `TRUE`, use attached variable labels.
#' @param TreatOrdinalAs How ordered predictors and covariates are handled.
#'
#' @return A list containing:
#' * `FormattedTable`: a report-facing `gt` table.
#' * `LargeTable`: a detailed `gt` table.
#' * `Results`: one row per predictor coefficient, with the tidy regression
#'   fields used by [MakeUnivariateRegressionTable()] plus `Events`,
#'   `Concordance`, `PH_PValue`, and (when requested) `FDR`.
#' * `ModelSummaries`: fitted models when `ReturnModels = TRUE`, otherwise
#'   `NULL`.
#' * `Metadata`: event coding, modeling engine, and analysis settings.
#'
#' @examples
#' lung <- survival::lung
#' lung$event <- lung$status == 2
#' attr(lung$age, "label") <- "Age"
#'
#' cox_results <- MakeCoxRegressionTable(
#'   data = lung,
#'   time_var = "time",
#'   event_vars = "event",
#'   predictor_vars = c("age", "sex")
#' )
#' cox_results$FormattedTable
#' PlotForestFromTable(cox_results$Results)
#'
#' @seealso [MakeUnivariateRegressionTable()], [PlotForestFromTable()],
#'   [ApplyFDRCorrection()]
#' @export
MakeCoxRegressionTable <- function(data,
    time_var,
    event_vars,
    predictor_vars,
    covariates = NULL,
    design = NULL,
    reference_levels = NULL,
    Standardize = FALSE,
    FDR = TRUE,
    FDRAlpha = 0.05,
    CheckPH = TRUE,
    ReturnModels = FALSE,
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical") {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (!is.character(time_var) || length(time_var) != 1 || is.na(time_var) || !nzchar(time_var)) {
    stop("time_var must be a single variable name.", call. = FALSE)
  }
  if (!is.character(event_vars) || length(event_vars) == 0 || anyNA(event_vars)) {
    stop("event_vars must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(predictor_vars) || length(predictor_vars) == 0 || anyNA(predictor_vars)) {
    stop("predictor_vars must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates))) {
    stop("covariates must be NULL or a character vector.", call. = FALSE)
  }
  ScidrCoxValidateFlag(Standardize, "Standardize")
  ScidrCoxValidateFlag(FDR, "FDR")
  ScidrCoxValidateFlag(CheckPH, "CheckPH")
  ScidrCoxValidateFlag(ReturnModels, "ReturnModels")
  ScidrCoxValidateFlag(Relabel, "Relabel")
  if (!is.numeric(FDRAlpha) || length(FDRAlpha) != 1 || is.na(FDRAlpha) ||
      FDRAlpha <= 0 || FDRAlpha >= 1) {
    stop("FDRAlpha must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.null(design) && !inherits(design, "survey.design")) {
    stop("design must be NULL or a survey.design object.", call. = FALSE)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Package 'survival' is required by MakeCoxRegressionTable(). ",
      "Install it with install.packages('survival').",
      call. = FALSE
    )
  }
  if (!is.null(design) && !requireNamespace("survey", quietly = TRUE)) {
    stop(
      "Package 'survey' is required when design is supplied. ",
      "Install it with install.packages('survey').",
      call. = FALSE
    )
  }
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop(
      "Package 'gt' is required by MakeCoxRegressionTable(). ",
      "Install it with install.packages('gt').",
      call. = FALSE
    )
  }

  all_model_vars <- unique(c(time_var, event_vars, predictor_vars, covariates))
  missing_vars <- setdiff(all_model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (!is.null(design)) {
    missing_design_vars <- setdiff(all_model_vars, names(design$variables))
    if (length(missing_design_vars) > 0) {
      stop(
        "The following variables were not found in design$variables: ",
        paste(missing_design_vars, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!is.numeric(data[[time_var]])) {
    stop("time_var must name a numeric variable.", call. = FALSE)
  }
  if (!is.null(design) && !is.numeric(design$variables[[time_var]])) {
    stop("time_var must name a numeric variable in design$variables.", call. = FALSE)
  }

  TreatOrdinalAs <- match.arg(
    TreatOrdinalAs,
    c("Categorical", "Continuous", "Both", "Exclude")
  )
  if (TreatOrdinalAs == "Both") {
    stop("TreatOrdinalAs = 'Both' is not meaningful for MakeCoxRegressionTable().", call. = FALSE)
  }
  ordinal_vars <- unique(c(predictor_vars, covariates))
  ordinal_data <- if (is.null(design)) data else design$variables
  ordinal_reference <- ConvertOrdinalToNumeric(
    ordinal_data,
    ordinal_vars,
    TreatOrdinalAs = "Categorical",
    ReturnMetadata = TRUE
  )
  if (TreatOrdinalAs == "Exclude" && length(ordinal_reference$ordinal_variables)) {
    stop(
      "TreatOrdinalAs = 'Exclude' cannot be used when ordinal model variables are explicitly supplied.",
      call. = FALSE
    )
  }

  display_labels <- ScidrDisplayLabels(data, all_model_vars, Relabel)
  data <- ConvertOrdinalToNumeric(
    data,
    ordinal_vars,
    TreatOrdinalAs = TreatOrdinalAs,
    ReturnMetadata = TRUE
  )$data
  if (!is.null(design)) {
    design$variables <- ConvertOrdinalToNumeric(
      design$variables,
      ordinal_vars,
      TreatOrdinalAs = TreatOrdinalAs,
      ReturnMetadata = TRUE
    )$data
  }

  reference_levels <- ScidrCoxReferenceLevels(reference_levels, predictor_vars)
  for (predictor in names(reference_levels)) {
    data[[predictor]] <- ScidrCoxApplyReference(
      data[[predictor]],
      predictor,
      reference_levels[[predictor]]
    )
    if (!is.null(design)) {
      design$variables[[predictor]] <- ScidrCoxApplyReference(
        design$variables[[predictor]],
        predictor,
        reference_levels[[predictor]]
      )
    }
  }

  modeling_data <- if (is.null(design)) data else design$variables
  event_info <- lapply(event_vars, function(event_var) {
    ScidrCoxPrepareEvent(modeling_data[[event_var]], event_var)
  })
  names(event_info) <- event_vars
  for (event_var in event_vars) {
    if (is.null(design)) {
      data[[event_var]] <- event_info[[event_var]]$Value
    } else {
      design$variables[[event_var]] <- event_info[[event_var]]$Value
    }
  }

  results_list <- list()
  model_list <- list()
  engine <- if (is.null(design)) "coxph" else "svycoxph"

  for (event_var in event_vars) {
    event_models <- list()
    for (predictor in predictor_vars) {
      model_vars <- unique(c(time_var, event_var, predictor, covariates))
      tryCatch(expr = {
        if (is.null(design)) {
          complete_rows <- stats::complete.cases(data[, model_vars, drop = FALSE])
          df_Model <- data[complete_rows, model_vars, drop = FALSE]
          df_Model <- ScidrCoxStandardizeData(
            df_Model,
            c(predictor, covariates),
            Standardize
          )
          n_model <- nrow(df_Model)
          events_model <- sum(df_Model[[event_var]])
          ScidrCoxValidateAnalysisSet(n_model, events_model)
          formula <- ScidrCoxFormula(time_var, event_var, c(predictor, covariates))
          model <- survival::coxph(
            formula = formula,
            data = df_Model,
            model = ReturnModels,
            x = ReturnModels
          )
        } else {
          complete_rows <- stats::complete.cases(design$variables[, model_vars, drop = FALSE])
          design_model <- design[complete_rows, ]
          design_model$variables <- ScidrCoxStandardizeData(
            design_model$variables,
            c(predictor, covariates),
            Standardize
          )
          n_model <- nrow(design_model$variables)
          events_model <- sum(design_model$variables[[event_var]])
          ScidrCoxValidateAnalysisSet(n_model, events_model)
          formula <- ScidrCoxFormula(time_var, event_var, c(predictor, covariates))
          model <- survey::svycoxph(formula = formula, design = design_model)
        }

        coefficient_table <- as.data.frame(summary(model)$coefficients)
        coefficient_table$Term <- rownames(coefficient_table)
        rownames(coefficient_table) <- NULL
        predictor_term <- ScidrQuoteFormulaNames(predictor)
        predictor_indices <- model$assign[[predictor_term]]
        if (is.null(predictor_indices)) {
          predictor_indices <- integer(0)
        }
        coefficient_table <- coefficient_table[predictor_indices, , drop = FALSE]
        if (nrow(coefficient_table) == 0) {
          stop("the fitted model did not estimate a coefficient for the predictor")
        }

        raw_estimate <- coefficient_table[["coef"]]
        std_error <- coefficient_table[["se(coef)"]]
        estimate <- exp(raw_estimate)
        conf_low <- exp(raw_estimate - stats::qnorm(0.975) * std_error)
        conf_high <- exp(raw_estimate + stats::qnorm(0.975) * std_error)
        p_value <- coefficient_table[["Pr(>|z|)"]]
        predictor_label <- display_labels[[predictor]]
        source_data <- if (is.null(design)) df_Model else design_model$variables
        categorical_predictor <- is.factor(source_data[[predictor]]) || is.character(source_data[[predictor]])
        level <- rep(NA_character_, nrow(coefficient_table))
        term_label <- rep(predictor_label, nrow(coefficient_table))
        if (categorical_predictor) {
          level <- substring(coefficient_table$Term, nchar(predictor_term) + 1L)
          level[level == coefficient_table$Term | level == ""] <- NA_character_
          term_label <- ifelse(
            is.na(level),
            predictor_label,
            paste0(predictor_label, " : ", level)
          )
        }

        concordance <- unname(summary(model)$concordance[["C"]])
        ph_p_value <- NA_real_
        if (is.null(design) && CheckPH) {
          ph_table <- survival::cox.zph(model)$table
          ph_p_value <- unname(ph_table["GLOBAL", "p"])
        }

        results_list[[length(results_list) + 1L]] <- data.frame(
          Outcome = event_var,
          OutcomeLabel = display_labels[[event_var]],
          OutcomeFamily = "survival",
          EffectType = "HR",
          Predictor = predictor,
          PredictorLabel = predictor_label,
          Term = coefficient_table$Term,
          Level = level,
          TermLabel = term_label,
          N = n_model,
          Events = events_model,
          Estimate = estimate,
          StdError = std_error,
          ConfLow = conf_low,
          ConfHigh = conf_high,
          PValue = p_value,
          Significant = p_value < 0.05,
          ReferenceValue = 1,
          Concordance = concordance,
          PH_PValue = ph_p_value,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        if (ReturnModels) {
          event_models[[predictor]] <- model
        }
      }, error = function(e) {
        stop(
          "Error processing event '", event_var, "' and predictor '", predictor,
          "': ", conditionMessage(e),
          call. = FALSE
        )
      })
    }
    if (ReturnModels) {
      model_list[[event_var]] <- event_models
    }
  }

  results <- dplyr::bind_rows(results_list)
  if (FDR) {
    results$FDR <- ApplyFDRCorrection(results$PValue)
    results$Significant <- !is.na(results$FDR) & results$FDR < FDRAlpha
  }

  outcomes_metadata <- dplyr::bind_rows(lapply(event_vars, function(event_var) {
    data.frame(
      Outcome = event_var,
      OutcomeLabel = display_labels[[event_var]],
      OutcomeFamily = "survival",
      ReferenceLevel = event_info[[event_var]]$ReferenceLevel,
      EventLevel = event_info[[event_var]]$EventLevel,
      Method = engine,
      stringsAsFactors = FALSE
    )
  }))
  metadata <- list(
    Outcomes = outcomes_metadata,
    AnalysisSettings = list(
      Method = engine,
      TimeVar = time_var,
      EventVars = event_vars,
      PredictorVars = predictor_vars,
      Covars = covariates,
      ReferenceLevels = reference_levels,
      Standardize = Standardize,
      FDR = FDR,
      FDRAlpha = FDRAlpha,
      CheckPH = CheckPH,
      ReturnModels = ReturnModels,
      TreatOrdinalAs = TreatOrdinalAs
    )
  )

  list(
    FormattedTable = ScidrCoxGtTable(results, formatted = TRUE),
    LargeTable = ScidrCoxGtTable(results, formatted = FALSE),
    Results = results,
    ModelSummaries = if (ReturnModels) model_list else NULL,
    Metadata = metadata
  )
}

ScidrCoxValidateFlag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

ScidrCoxReferenceLevels <- function(reference_levels, predictor_vars) {
  if (is.null(reference_levels)) return(character(0))
  if (!is.list(reference_levels) && !is.character(reference_levels)) {
    stop("reference_levels must be NULL, a named list, or a named character vector.", call. = FALSE)
  }
  if (is.null(names(reference_levels)) || any(names(reference_levels) == "") || anyDuplicated(names(reference_levels))) {
    stop("reference_levels must have unique, non-empty predictor names.", call. = FALSE)
  }
  values <- unlist(reference_levels, recursive = FALSE, use.names = TRUE)
  if (length(values) != length(reference_levels) ||
      !is.character(values) || anyNA(values) || any(!nzchar(values))) {
    stop("Each reference_levels entry must contain one non-missing character level.", call. = FALSE)
  }
  unknown <- setdiff(names(values), predictor_vars)
  if (length(unknown) > 0) {
    stop("reference_levels names not found in predictor_vars: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  values
}

ScidrCoxApplyReference <- function(x, predictor, reference) {
  if (is.character(x)) x <- factor(x)
  if (!is.factor(x)) {
    stop("reference_levels can only be set for categorical predictor '", predictor, "'.", call. = FALSE)
  }
  if (!reference %in% levels(x)) {
    stop(
      "Reference level '", reference, "' was not found for predictor '", predictor, "'.",
      call. = FALSE
    )
  }
  stats::relevel(x, ref = reference)
}

ScidrCoxPrepareEvent <- function(x, event_var) {
  non_missing <- x[!is.na(x)]
  if (is.numeric(x)) {
    if (length(non_missing) == 0 || any(!non_missing %in% c(0, 1))) {
      stop("Event variable '", event_var, "' must be coded only 0/1.", call. = FALSE)
    }
    return(list(Value = as.numeric(x), ReferenceLevel = "0", EventLevel = "1"))
  }
  if (is.logical(x)) {
    if (length(non_missing) == 0) {
      stop("Event variable '", event_var, "' contains no observed values.", call. = FALSE)
    }
    return(list(Value = as.integer(x), ReferenceLevel = "FALSE", EventLevel = "TRUE"))
  }
  if (is.factor(x)) {
    observed <- droplevels(non_missing)
    if (nlevels(observed) != 2) {
      stop("Event variable '", event_var, "' must be a two-level factor.", call. = FALSE)
    }
    event_levels <- levels(observed)
    value <- rep(NA_real_, length(x))
    value[!is.na(x)] <- as.numeric(as.character(x[!is.na(x)]) == event_levels[[2]])
    return(list(
      Value = value,
      ReferenceLevel = event_levels[[1]],
      EventLevel = event_levels[[2]]
    ))
  }
  stop(
    "Event variable '", event_var, "' must be numeric 0/1, logical, or a two-level factor.",
    call. = FALSE
  )
}

ScidrCoxValidateAnalysisSet <- function(n_model, events_model) {
  if (n_model == 0) stop("no complete cases remain")
  if (events_model == 0) stop("no events remain after removing missing values")
  if (events_model == n_model) stop("no censored observations remain after removing missing values")
}

ScidrCoxStandardizeData <- function(data, variables, Standardize) {
  if (!Standardize) return(data)
  numeric_variables <- variables[vapply(data[variables], is.numeric, logical(1))]
  for (variable in numeric_variables) {
    variable_sd <- stats::sd(data[[variable]])
    if (!is.finite(variable_sd) || variable_sd == 0) {
      stop("numeric model variable '", variable, "' cannot be standardized because it has zero variance")
    }
    data[[variable]] <- as.numeric(scale(data[[variable]]))
  }
  data
}

ScidrCoxFormula <- function(time_var, event_var, model_terms) {
  quote_name <- function(x) deparse1(as.name(x), backtick = TRUE)
  response <- paste0(
    "survival::Surv(", quote_name(time_var), ", ", quote_name(event_var), ")"
  )
  stats::reformulate(
    termlabels = vapply(model_terms, quote_name, character(1)),
    response = response,
    env = parent.frame()
  )
}

ScidrCoxGtTable <- function(results, formatted = TRUE) {
  results <- results %>%
    dplyr::mutate(RowKey = paste(.data$Predictor, .data$Term, sep = "\r"))
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
      dplyr::filter(.data$Outcome == outcome)
    outcome_results <- outcome_results[match(row_data$RowKey, outcome_results$RowKey), , drop = FALSE]
    column_prefix <- paste0("Outcome", outcome_index)
    outcome_label <- results$OutcomeLabel[match(outcome, results$Outcome)]

    if (formatted) {
      effect_column <- paste0(column_prefix, "_EffectCI")
      p_column <- paste0(column_prefix, "_P")
      table_data[[effect_column]] <- dplyr::case_when(
        is.na(outcome_results$Estimate) ~ NA_character_,
        TRUE ~ paste0(
          formatC(outcome_results$Estimate, digits = 3, format = "fg"),
          " (", formatC(outcome_results$ConfLow, digits = 3, format = "fg"),
          ", ", formatC(outcome_results$ConfHigh, digits = 3, format = "fg"), ")",
          ifelse(outcome_results$Significant, ScidrPValueStars(
            if ("FDR" %in% names(outcome_results)) outcome_results$FDR else outcome_results$PValue,
            ns_label = ""
          ), "")
        )
      )
      table_data[[p_column]] <- dplyr::case_when(
        is.na(outcome_results$PValue) ~ NA_character_,
        outcome_results$PValue < 0.001 ~ "<0.001",
        TRUE ~ formatC(outcome_results$PValue, digits = 2, format = "fg")
      )
      outcome_columns <- c(effect_column, p_column)
      column_labels[[effect_column]] <- "HR (95% CI)"
      column_labels[[p_column]] <- "p-value"
      if ("FDR" %in% names(outcome_results)) {
        fdr_column <- paste0(column_prefix, "_FDR")
        table_data[[fdr_column]] <- outcome_results$FDR
        outcome_columns <- c(outcome_columns, fdr_column)
        column_labels[[fdr_column]] <- "FDR"
        numeric_columns <- c(numeric_columns, fdr_column)
        significant_rows[[fdr_column]] <- which(outcome_results$Significant)
      }
      significant_rows[[effect_column]] <- which(outcome_results$Significant)
      significant_rows[[p_column]] <- which(outcome_results$Significant)
    } else {
      large_columns <- c(
        "N", "Events", "Estimate", "StdError", "ConfLow", "ConfHigh",
        "PValue", if ("FDR" %in% names(outcome_results)) "FDR", "Concordance", "PH_PValue"
      )
      large_labels <- c(
        N = "N", Events = "Events", Estimate = "HR", StdError = "SE (log HR)",
        ConfLow = "95% CI Low", ConfHigh = "95% CI High", PValue = "p-value",
        FDR = "FDR", Concordance = "Concordance", PH_PValue = "Global PH p-value"
      )
      outcome_columns <- paste0(column_prefix, "_", large_columns)
      for (column_index in seq_along(large_columns)) {
        table_data[[outcome_columns[[column_index]]]] <- outcome_results[[large_columns[[column_index]]]]
      }
      column_labels <- c(
        column_labels,
        as.list(stats::setNames(large_labels[large_columns], outcome_columns))
      )
      numeric_columns <- c(numeric_columns, outcome_columns)
    }
    spanners[[column_prefix]] <- outcome_columns
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
  if (length(numeric_columns) > 0) {
    out <- gt::fmt_number(out, columns = numeric_columns, decimals = 3)
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
  }
  out
}
