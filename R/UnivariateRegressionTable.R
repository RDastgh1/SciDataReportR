#' Univariate Regression Table
#'
#' Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified).
#'
#' @param Data Dataframe containing the variables
#' @param OutcomeVars Character vector of outcome variable names
#' @param PredictorVars Character vector of predictor variable names
#' @param Covars Character vector of covariate variable names (default: NULL)
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
#'   - ModelSummaries: A list of fitted model objects
#'   - Metadata: Outcome families and analysis settings
#' @export
#'
UnivariateRegressionTable <- function (
    Data,
    OutcomeVars,
    PredictorVars,
    Covars = NULL,
    Standardize = FALSE,
    Method = c("auto", "lm", "logistic"),
    LogisticExponentiate = TRUE
)
{
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
    outcome_metadata[[YVar]] <- data.frame(
      Outcome = YVar,
      OutcomeLabel = sjlabelled::get_label(Data[[YVar]], def.value = YVar) %>% as.character(),
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
              ModelSummaries = Wide_mod_list, Metadata = metadata))
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
