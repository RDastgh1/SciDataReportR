#' Univariate Regression Table
#'
#' Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified).
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
#' @export
#'
MakeUnivariateRegressionTable <- function(data,
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
      tryCatch(expr = {
        f <- as.formula(paste(YVar, paste(c(xVar, Covars),
                                          collapse = "+"), sep = " ~ "))
        if (Standardize) {
          ModelData <- model_data_full %>% select(all_of(c(YVar, xVar, Covars)))
          numeric_cols <- sapply(ModelData, is.numeric)
          if (outcome_family == "logistic") {
            numeric_cols[YVar] <- FALSE
          }
          if (any(numeric_cols)) {
            ModelData[, numeric_cols] <- scale(ModelData[, numeric_cols, drop = FALSE])
          }
        }else {
          ModelData <- model_data_full
        }
        mod <- if (outcome_family == "logistic") {
          stats::glm(formula = f, data = ModelData, family = stats::binomial())
        } else {
          stats::lm(formula = f, data = ModelData)
        }
        if (ReturnModels) {
          mod_list[[xVar]] <- mod
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

        labels <- sjlabelled::get_label(Data, def.value = colnames(Data))
        labels <- labels[c(PredictorVars, OutcomeVars)]
        labels <- stats::setNames(as.character(labels), c(PredictorVars, OutcomeVars))
        predictor_label <- labels[[xVar]]
        effect_type <- ifelse(
          outcome_family == "logistic" && LogisticExponentiate,
          "Odds ratio",
          "Estimate"
        )

        model_terms <- attr(stats::terms(mod), "term.labels")
        keep_terms <- model_terms[model_terms %!in% Covars]
        keep_rows <- coefficient_table$Term %in% keep_terms
        if (length(keep_terms) > 0) {
          keep_rows <- keep_rows | startsWith(coefficient_table$Term, paste0(keep_terms, ""))
        }
        coefficient_table <- coefficient_table[keep_rows, , drop = FALSE]

        if (nrow(coefficient_table) == 0) {
          return(NULL)
        }

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
          level <- sub(paste0("^", xVar), "", coefficient_table$Term)
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
          Outcome = YVar,
          OutcomeLabel = sjlabelled::get_label(Data[[YVar]], def.value = YVar) %>% as.character(),
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
        stop(paste("Error processing", YVar, "and", xVar,
                   ": ", conditionMessage(e)))
      })
    }
    if (ReturnModels) {
      Wide_mod_list[[YVar]] <- mod_list
    }
    outcome_label <- sjlabelled::get_label(Data[[YVar]], def.value = YVar) %>% as.character()
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
  FinalTable <- ScidrUnivariateGtTable(results, formatted = FALSE)
  FinalFormattedTable <- ScidrUnivariateGtTable(results, formatted = TRUE)

  metadata <- list(
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

ScidrUnivariateGtTable <- function(results, formatted = TRUE) {
  table_data <- results %>%
    dplyr::mutate(
      Estimate_CI = paste0(
        formatC(.data$Estimate, digits = 3, format = "fg"),
        " (",
        formatC(.data$ConfLow, digits = 3, format = "fg"),
        ", ",
        formatC(.data$ConfHigh, digits = 3, format = "fg"),
        ")"
      ),
      P = dplyr::case_when(
        is.na(.data$PValue) ~ NA_character_,
        .data$PValue < 0.001 ~ "<0.001",
        TRUE ~ formatC(.data$PValue, digits = 2, format = "fg")
      )
    )

  if (formatted) {
    table_data <- table_data %>%
      dplyr::select(
        OutcomeLabel,
        TermLabel,
        EffectType,
        N,
        Estimate_CI,
        P,
        Significant
      )

    out <- gt::gt(table_data, groupname_col = "OutcomeLabel") %>%
      gt::cols_label(
        TermLabel = "Variable",
        EffectType = "Effect",
        N = "N",
        Estimate_CI = "Estimate (95% CI)",
        P = "p-value"
      ) %>%
      gt::cols_hide(columns = "Significant") %>%
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
        TermLabel,
        EffectType,
        N,
        Estimate,
        StdError,
        ConfLow,
        ConfHigh,
        PValue
      )

    out <- gt::gt(table_data, groupname_col = "OutcomeLabel") %>%
      gt::cols_label(
        Outcome = "Outcome",
        TermLabel = "Variable",
        EffectType = "Effect",
        N = "N",
        Estimate = "Estimate",
        StdError = "SE",
        ConfLow = "95% CI Low",
        ConfHigh = "95% CI High",
        PValue = "p-value"
      ) %>%
      gt::fmt_number(
        columns = c("Estimate", "StdError", "ConfLow", "ConfHigh"),
        decimals = 3
      ) %>%
      gt::fmt_number(columns = "PValue", decimals = 3)
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
