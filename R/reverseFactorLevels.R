#' Reverse Levels of Categorical Factors
#'
#' This function reverses the levels of specified categorical variables in a given dataframe.
#'
#' @param df A dataframe containing the categorical variables to be reversed.
#' @param variables A character vector of column names in the dataframe to reverse levels.
#'   These columns must be categorical factors.
#'
#' @return A dataframe with the levels of the specified factors reversed.
#'   Columns not specified in `variables` remain unchanged.
#' @export
reverseFactorLevels <- function(df, variables) {
  # Check if all specified variables are in the dataframe
  missing_vars <- setdiff(variables, colnames(df))
  if (length(missing_vars) > 0) {
    stop("The following variables are not in the dataframe: ", paste(missing_vars, collapse = ", "))
  }

  # Check if all specified variables are factors
  non_factors <- variables[!sapply(df[variables], is.factor)]
  if (length(non_factors) > 0) {
    stop("The following variables are not factors: ", paste(non_factors, collapse = ", "))
  }

  # Reverse levels of specified variables
  df[variables] <- lapply(df[variables], function(var) {
    factor(var, levels = rev(levels(var)))
  })

  return(df)
}
