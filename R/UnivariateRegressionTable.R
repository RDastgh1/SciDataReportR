#' Univariate Regression Table
#'
#' Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified).
#'
#' @param Data Dataframe containing the variables
#' @param OutcomeVars Character vector of outcome variable names
#' @param PredictorVars Character vector of predictor variable names
#' @param Covars Character vector of covariate variable names (default: NULL)
#' @param Standardize Logical indicating whether to standardize numeric variables (default: FALSE)
#' @importFrom sjlabelled get_label set_label
#' @return A list containing:
#'   - FormattedTable: A merged table with formatted regression results
#'   - LargeTable: A merged table with unformatted regression results
#'   - ModelSummaries: A list of lm model summaries
#' @export
#'
UnivariateRegressionTable <- function (Data, OutcomeVars, PredictorVars, Covars = NULL, Standardize = FALSE)
{
  Wide_tbl_list <- list()
  Wide_mod_list <- list()
  Wide_tblformatted_list <- list()
  for (yVarIndex in seq_along(OutcomeVars)) {
    YVar <- OutcomeVars[yVarIndex]
    tbl_list <- list()
    mod_list <- list()
    tblformatted_list <- list()
    for (xVarIndex in seq_along(PredictorVars)) {
      xVar <- PredictorVars[xVarIndex]
      tryCatch(expr = {
        f <- as.formula(paste(YVar, paste(c(xVar, Covars),
                                          collapse = "+"), sep = " ~ "))
        if (Standardize) {
          ModelData <- Data %>% select(all_of(c(YVar,
                                                xVar, Covars)))
          numeric_cols <- sapply(ModelData, is.numeric)
          ModelData[, numeric_cols] <- scale(ModelData[,
                                                       numeric_cols])
          mod <- lm(formula = f, data = ModelData)
        }else {
          mod <- lm(formula = f, data = Data)
        }
        mod_list[[xVar]] <- mod
        labels <- get_label(Data, def.value = colnames(Data))
        labels <- labels[c(PredictorVars, OutcomeVars)]
        label_list <- setNames(as.list(labels), c(PredictorVars,
                                                  OutcomeVars))
        modTableP <- tbl_regression(mod, pvalue_fun = ~style_pvalue(.x,
                                                                    digits = 2), label = label_list[xVar]) %>%
          bold_p() %>% bold_labels() %>% italicize_levels()
        modTableP$table_body <- modTableP$table_body %>%
          filter(variable %!in% Covars)
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
  return(list(FormattedTable = FinalFormattedTable, LargeTable = FinalTable,
              ModelSummaries = Wide_mod_list))
}
