#' Univariate Regression Table
#'
#' Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified).
#'
#' @param Data Dataframe containing the variables
#' @param YVars Character vector of outcome variable names
#' @param XVars Character vector of predictor variable names
#' @param Covars Character vector of covariate variable names (default: NULL)
#' @param Standardize Logical indicating whether to standardize numeric variables (default: FALSE)
#' @importFrom sjlabelled get_label set_label
#' @return A list containing:
#'   - FormattedTable: A merged table with formatted regression results
#'   - LargeTable: A merged table with unformatted regression results
#'   - ModelSummaries: A list of lm model summaries
#' @export
#'
UnivariateRegressionTable <- function(Data, YVars, XVars, Covars = NULL, Standardize = FALSE) {
  # Initialize empty lists to store results
  Wide_tbl_list <- list()
  Wide_mod_list <- list()
  Wide_tblformatted_list <- list()

  # Loop through each outcome variable
  for (yVarIndex in seq_along(YVars)) {
    YVar <- YVars[yVarIndex]

    # Initialize empty lists to store results for this outcome variable
    tbl_list <- list()
    mod_list <- list()
    tblformatted_list <- list()

    # Loop through each predictor variable
    for (xVarIndex in seq_along(XVars)) {
      xVar <- XVars[xVarIndex]

      tryCatch(
        expr = {
          # Create formula for regression model
          f <- as.formula(paste(YVar, paste(c(xVar, Covars), collapse = "+"), sep = " ~ "))

          # Standardize numeric variables if specified
          if (Standardize) {
            ModelData <- Data %>% select(all_of(c(YVar, xVar, Covars)))
            numeric_cols <- sapply(ModelData, is.numeric)
            ModelData[, numeric_cols] <- scale(ModelData[, numeric_cols])
            mod <- lm(formula = f, data = ModelData)
          } else {
            mod <- lm(formula = f, data = Data)
          }

          # Store model summary
          mod_list[[xVar]] <- mod

          # Get variable labels
          labels <- get_label(Data, def.value = colnames(Data))
          labels <- labels[c(XVars, YVars)]

          # Create regression table with labels
          modTableP <- tbl_regression(mod,
                                      pvalue_fun = ~style_pvalue(.x, digits = 2),
                                      label = labels) %>%
            bold_p() %>%
            bold_labels() %>%
            italicize_levels()
          modTableP$table_body <- modTableP$table_body %>% filter(variable %notin% Covars)
          modTableP$table_body$var_label <- as.character(modTableP$table_body$var_label)

          # Store regression table
          tbl_list[[xVar]] <- modTableP

          # Create formatted regression table
          modTableCombined <- modTableP %>%
            add_significance_stars() %>%
            modify_table_styling(columns = "estimate", cols_merge_pattern = "{estimate} ({std.error}){stars}")

          # Store formatted regression table
          tblformatted_list[[xVar]] <- modTableCombined
        },
        error = function(e) {
          stop(paste("Error processing", YVar, "and", xVar, ": ", e$message))
        }
      )
    }

    # Stack tables for this outcome variable
    Wide_tbl_list[[YVar]] <- tbl_stack(tbl_list) %>% remove_row_type(type = "reference")
    Wide_mod_list[[YVar]] <- mod_list
    Wide_tblformatted_list[[YVar]] <- tbl_stack(tblformatted_list) %>% remove_row_type(type = "reference")
  }

  # Get official names for outcome variables
  # Get official names for outcome variables
  s <- sjlabelled::get_label(Data[YVars], def.value = YVars)

  # Merge tables
  FinalTable <- tbl_merge(Wide_tbl_list, tab_spanner = unname(s))
  FinalFormattedTable <- tbl_merge(Wide_tblformatted_list, tab_spanner = unname(s))

  # Return results
  return(list(FormattedTable = FinalFormattedTable, LargeTable = FinalTable, ModelSummaries = Wide_mod_list))
}
