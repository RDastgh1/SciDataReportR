#' Get Categorical Variables
#'
#' Extracts categorical variables from a data frame.
#'
#' @param DataFrame The data frame from which to extract categorical variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of categorical variables.
#' @importFrom dplyr select where
#' @export
getCatVars <- function(DataFrame, Ordinal = TRUE) {
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
