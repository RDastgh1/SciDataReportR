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

  # Reverse levels of specified variables, retaining labels
  df[variables] <- lapply(variables, function(var_name) {
    var <- df[[var_name]]

    # Get the label if it exists
    var_label <- attr(var, "label")

    # Reverse the factor levels
    reversed_var <- factor(var, levels = rev(levels(var)))

    # Reassign the label
    if (!is.null(var_label)) {
      attr(reversed_var, "label") <- var_label
    }

    return(reversed_var)
  })

  return(df)
}
