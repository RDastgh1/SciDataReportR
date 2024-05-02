#' ConvertOrdinalToNumeric
#'
#' Convert ordinal variables in a dataframe to numeric if they contain numeric values in their character representation.
#'
#' @param Data The dataframe containing the variables.
#' @param Variables A character vector specifying the names of variables to consider. If NULL, all columns of the dataframe will be considered.
#' @importFrom sjlabelled get_label set_label
#' @return The dataframe with ordinal variables potentially converted to numeric.
#'
#' @export
ConvertOrdinalToNumeric <- function(Data, Variables = NULL) {

  # If Variables argument is NULL, consider all columns of the dataframe
  if(is.null(Variables)){
    Variables = colnames(Data)
  }

  # Identify ordered variables
  orderedVars <- names(Data[Variables])[sapply(Data, is.ordered)] %>% na.omit()

  # get original labels to reset them later
  l <- get_label(Data)

  # Iterate through ordered variables
  for (col in orderedVars) {
    # Convert ordered variable to character
    x <- as.character(Data[[col]])

    # Check if character values are numeric
    numeric_values <- grepl("^\\d+\\.?\\d*$", x)

    # If all character values are numeric, convert to numeric
    if(sum(!numeric_values) == 0){
      # preserve label


      Data[[col]] <- as.numeric(x)
    }
  }

  # Readd labels
  set_label(Data)<- l
  return(Data)
}
