#' Get Numeric Variables
#'
#' Extracts numeric variables from a data frame.
#'
#' @param DataFrame The data frame from which to extract numeric variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of numeric variables.
#' @importFrom dplyr select where
#' @export
getNumVars <- function(DataFrame, Ordinal = FALSE) {
  NumVars <- DataFrame %>% select(where(is.numeric)) %>% names()
  if (Ordinal) {
    NumVars <- DataFrame %>% select(where(is.numeric) | where(is.ordered)) %>% names()
  }
  return(NumVars)
}
