#' Convert ordinal variables to numeric
#'
#' Convert ordinal variables in a dataframe to numeric if they contain numeric values in their character representation.
#'
#' @param data The dataframe containing the variables.
#' @param variables A character vector specifying the names of variables to consider. If NULL, all columns of the dataframe will be considered.
#' @importFrom sjlabelled get_label set_label
#' @return The dataframe with ordinal variables potentially converted to numeric.
#'
#' @examples
#' # An ordered factor with numeric levels, and one with non-numeric levels
#' df <- data.frame(
#'   id     = 1:5,
#'   likert = factor(c("1", "2", "3", "2", "1"),
#'                   levels = c("1", "2", "3"), ordered = TRUE),
#'   grade  = factor(c("A", "B", "A", "C", "B"),
#'                   levels = c("A", "B", "C"), ordered = TRUE)
#' )
#'
#' out <- ConvertOrdinalToNumeric(df)
#'
#' # likert becomes numeric; grade stays an ordered factor (levels are not numeric)
#' sapply(out, class)
#'
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @export
ConvertOrdinalToNumeric <- function(data,
    variables = NULL,
    Data = lifecycle::deprecated(),
    Variables = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "ConvertOrdinalToNumeric(Data)", "ConvertOrdinalToNumeric(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "ConvertOrdinalToNumeric(Variables)", "ConvertOrdinalToNumeric(variables)")
    variables <- Variables
  }
  Variables <- variables


  # If Variables argument is NULL, consider all columns of the dataframe
  if(is.null(Variables)){
    Variables = colnames(Data)
  }

  # Identify ordered variables
  orderedVars <- names(Data[Variables])[sapply(Data, is.ordered)] %>% na.omit()

  # get original labels to reset them later
  l <- sjlabelled::get_label(Data)

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
  sjlabelled::set_label(Data)<- l
  return(Data)
}
