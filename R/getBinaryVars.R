#' @title Identify Binary Variables
#' @description This function identifies and returns a list of binary variables
#' in a dataframe. Binary variables are defined as having exactly two unique values
#' or levels. The function supports options for handling ordinal factors and
#' revalued data.
#'
#' @param DataFrame A dataframe to analyze for binary variables.
#' @param Ordinal Logical. If TRUE, ordinal factors are included in the search
#' for binary variables. Default is TRUE.
#' @param Revalued Logical. If TRUE, the function checks factors and their levels;
#' otherwise, it checks for variables with two unique values.
#' @return A character vector containing the names of binary variables.
#' @export
getBinaryVars <- function(DataFrame, Ordinal = TRUE, Revalued = TRUE) {
  if (Revalued) {
    # For revalued data, focus on factors with exactly 2 levels
    if (Ordinal) {
      # Include all factors (ordered and unordered)
      CatVars <- DataFrame %>%
        select(where(is.factor)) %>%
        names()
    } else {
      # Exclude ordered factors
      CatVars <- DataFrame %>%
        select(where(is.factor) & !where(is.ordered)) %>%
        names()
    }

    # Select binary variables from categorical variables (2 levels)
    BinaryVars <- DataFrame %>%
      select(all_of(CatVars)) %>%
      select(where(~ length(levels(.)) == 2)) %>%
      names()
  } else {
    # For non-revalued data, check for variables with 2 unique values
    BinaryVars <- DataFrame %>%
      select(where(~ length(unique(.)) == 2)) %>%
      names()
  }

  # Return the names of binary variables
  return(BinaryVars)
}
