#' Get Categorical Variables
#'
#' Extracts categorical variables from a data frame.
#'
#' @param data The data frame from which to extract categorical variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of categorical variables.
#' @importFrom dplyr select where
#' @importFrom magrittr "%>%"
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
getCatVars <- function(data,
    Ordinal = TRUE,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "getCatVars(DataFrame)", "getCatVars(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  if (Ordinal) {
    CatVars <- DataFrame %>%
      select(where(is.factor)) %>%
      names()
  } else {
    CatVars <- DataFrame %>%
      select(where(is.factor) & !where(is.ordered)) %>%
      names()
  }
  return(CatVars)
}
