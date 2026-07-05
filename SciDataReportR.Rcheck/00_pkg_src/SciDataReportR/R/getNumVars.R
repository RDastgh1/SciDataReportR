#' Get Numeric Variables
#'
#' Extracts numeric variables from a data frame.
#'
#' @param data The data frame from which to extract numeric variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of numeric variables.
#' @importFrom dplyr select where
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
getNumVars <- function(data,
    Ordinal = FALSE,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "getNumVars(DataFrame)", "getNumVars(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  NumVars <- DataFrame %>% select(where(is.numeric)) %>% names()
  if (Ordinal) {
    NumVars <- DataFrame %>% select(where(is.numeric) | where(is.ordered)) %>% names()
  }
  return(NumVars)
}
