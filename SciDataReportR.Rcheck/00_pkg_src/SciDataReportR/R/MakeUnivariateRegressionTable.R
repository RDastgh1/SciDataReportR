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
#' @importFrom sjlabelled get_label set_label
#' @return A list containing:
#'   - FormattedTable: A merged table with formatted regression results
#'   - LargeTable: A merged table with unformatted regression results
#'   - Results: A tidy dataframe with one row per estimated term. Columns:
#'     `Outcome`, `OutcomeLabel`, `OutcomeFamily`, `EffectType`, `Predictor`,
#'     `PredictorLabel`, `Term`, `Level`, `TermLabel`, `N`, `Estimate`,
#'     `StdError`, `ConfLow`, `ConfHigh`, `PValue`, `Significant`, and
#'     `ReferenceValue`. This dataframe can be filtered and passed directly to
#'     [PlotForestFromTable()].
#'   - ModelSummaries: A list of fitted model objects
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

  all_model_vars <- unique(c(OutcomeVars, PredictorVars, Covars))
  missing_vars <- setdiff(all_model_vars, names(Data))
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in Data: ", paste(missing_vars, collapse = ", "))
  }

  Wide_tbl_list <- list()
  Wide_mod_list <- list()
  Wide_tblformatted_list <- list()
  results_list <- list()
  outcome_metadata <- list()

  for (yVarIndex in seq_along(OutcomeVars)) {
    YVar <- OutcomeVars[yVarIndex]
    outcome_family <- ScidrUnivariateOutcomeFamily(Data[[YVar]], YVar, Method)
    outcome_info <- ScidrUnivariatePrepareOutcome(Data[[YVar]], YVar, outcome_family)
    model_data_full <- Data
    model_data_full[[YVar]] <- outcome_info$Value

    tbl_list <- list()
    mod_list <- list()
    tblformatted_list <- list()
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
        mod_list[[xVar]] <- mod
        labels <- get_label(Data, def.value = colnames(Data))
        labels <- labels[c(PredictorVars, OutcomeVars)]
        label_list <- setNames(as.list(labels), c(PredictorVars,
                                                  OutcomeVars))
        modTableP <- tbl_regression(mod, pvalue_fun = ~style_pvalue(.x,
                                                                    digits = 2), label = label_list[xVar],
                                     exponentiate = outcome_family == "logistic" && LogisticExponentiate) %>%
          bold_p() %>% bold_labels() %>% italicize_levels()
        modTableP$table_body <- modTableP$table_body %>%
          filter(variable %!in% Covars)
        modTableP$table_body$outcome_family <- outcome_family
        modTableP$table_body$effect_type <- ifelse(
          outcome_family == "logistic" && LogisticExponentiate,
          "Odds ratio",
          "Estimate"
        )
        modTableP$table_body$var_label <- as.character(modTableP$table_body$var_label)
        tbl_list[[xVar]] <- modTableP
        modTableCombined <- modTableP %>% add_significance_stars() %>%
          modify_table_styling(columns = "estimate",
                               cols_merge_pattern = "{estimate} ({std.error}){stars}")
        tblformatted_list[[xVar]] <- modTableCombined
      }, error = function(e) {
        stop(paste("Error processing", YVar, "and", xVar,
                   ": ", e$message))
      })
    }
    Wide_tbl_list[[YVar]] <- tbl_stack(tbl_list) %>% remove_row_type(type = "reference")
    Wide_mod_list[[YVar]] <- mod_list
    Wide_tblformatted_list[[YVar]] <- tbl_stack(tblformatted_list) %>%
      remove_row_type(type = "reference")
    outcome_label <- sjlabelled::get_label(Data[[YVar]], def.value = YVar) %>% as.character()
    results_list[[YVar]] <- ScidrUnivariateResultsFromBody(
      Wide_tbl_list[[YVar]]$table_body,
      outcome = YVar,
      outcome_label = outcome_label
    )
    outcome_metadata[[YVar]] <- data.frame(
      Outcome = YVar,
      OutcomeLabel = outcome_label,
      OutcomeFamily = outcome_family,
      ReferenceLevel = outcome_info$ReferenceLevel,
      EventLevel = outcome_info$EventLevel,
      stringsAsFactors = FALSE
    )
  }
  s <- vapply(
    OutcomeVars,
    function(var) {
      # for each var: if it has an sjlabel, use it; otherwise fall back to the var name
      sjlabelled::get_label(Data[[var]], def.value = var) %>% as.character()
    },
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )

  FinalTable <- tbl_merge(Wide_tbl_list, tab_spanner = unname(s))
  FinalFormattedTable <- tbl_merge(Wide_tblformatted_list,
                                   tab_spanner = unname(s))

  metadata <- list(
    Outcomes = dplyr::bind_rows(outcome_metadata),
    AnalysisSettings = list(
      Method = Method,
      OutcomeVars = OutcomeVars,
      PredictorVars = PredictorVars,
      Covars = Covars,
      Standardize = Standardize,
      LogisticExponentiate = LogisticExponentiate
    )
  )

  return(list(FormattedTable = FinalFormattedTable, LargeTable = FinalTable,
              Results = dplyr::bind_rows(results_list),
              ModelSummaries = Wide_mod_list, Metadata = metadata))
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
    Data = lifecycle::deprecated(),
    OutcomeVars = lifecycle::deprecated(),
    PredictorVars = lifecycle::deprecated(),
    Covars = lifecycle::deprecated()) {
  lifecycle::deprecate_soft("20.5.0", "UnivariateRegressionTable()", "MakeUnivariateRegressionTable()")
  call <- match.call()
  call[[1L]] <- MakeUnivariateRegressionTable
  eval.parent(call)
}

ScidrUnivariateResultsFromBody <- function(body, outcome, outcome_label) {
  rows <- body[!is.na(body$estimate), , drop = FALSE]
  n_col <- intersect(c("N_obs", "N"), names(rows))
  data.frame(
    Outcome = outcome,
    OutcomeLabel = outcome_label,
    OutcomeFamily = rows$outcome_family,
    EffectType = rows$effect_type,
    Predictor = rows$variable,
    PredictorLabel = rows$var_label,
    Term = if ("term" %in% names(rows)) rows$term else NA_character_,
    Level = ifelse(rows$row_type == "level", rows$label, NA_character_),
    TermLabel = ifelse(
      !is.na(rows$label) & rows$label != rows$var_label,
      paste0(rows$var_label, " : ", rows$label),
      rows$var_label
    ),
    N = if (length(n_col) > 0) as.numeric(rows[[n_col[[1]]]]) else NA_real_,
    Estimate = rows$estimate,
    StdError = rows$std.error,
    ConfLow = rows$conf.low,
    ConfHigh = rows$conf.high,
    PValue = rows$p.value,
    Significant = rows$p.value < 0.05,
    ReferenceValue = ifelse(rows$effect_type == "Odds ratio", 1, 0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
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
