#' Create a Mapping Table for Binary Variables
#'
#' This function identifies binary variables within a given dataframe and generates a mapping table.
#' The mapping table includes the variable name, its associated label, and the positive value.
#'
#' @param Data A dataframe containing the data to analyze.
#' @param CatVars A character vector of column names representing categorical variables in the dataframe.
#' @return A dataframe with columns: `VariableName`, `Label`, and `PositiveValue` for binary variables found in `CatVars`.
#' @export
createBinaryMapping <- function(Data, CatVars) {
  # Check for binary variables and create the mapping table
  binaryMappingTable <- do.call(rbind, lapply(CatVars, function(var) {
    if (is.factor(Data[[var]])){
      unique_levels <- levels(na.omit(Data[[var]]))
    } else{
    unique_levels <- unique(na.omit(Data[[var]]))
    }
    if (length(unique_levels) == 2) {
      # Get the label for the variable (default to variable name if no label is set)
      label <- sjlabelled::get_label(Data[[var]], def.value = var)
      # Assign the second unique level as the positive value
      positive_value <- unique_levels[2]
      return(data.frame(
        VariableName = var,
        Label = label,
        PositiveValue = positive_value,
        stringsAsFactors = FALSE
      ))
    } else {
      return(NULL) # Exclude non-binary variables
    }
  }))

  # Handle the case where no binary variables are found
  if (is.null(binaryMappingTable)) {
    stop("No binary variables found in the provided CatVars.")
  }

  return(binaryMappingTable)
}
