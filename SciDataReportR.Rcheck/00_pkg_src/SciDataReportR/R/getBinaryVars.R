#' @title Identify Binary Variables
#' @description This function identifies and returns a list of binary variables
#' in a dataframe. Binary variables are defined as having exactly two unique values
#' or levels. The function supports options for handling ordinal factors and
#' revalued data.
#'
#' @param data A dataframe to analyze for binary variables.
#' @param Ordinal Logical. If TRUE, ordinal factors are included in the search
#' for binary variables. Default is TRUE.
#' @param Revalued Logical. If TRUE, the function checks factors and their levels;
#' otherwise, it checks for variables with two unique values.
#' @return A character vector containing the names of binary variables.
#' @seealso [createBinaryMapping()] to fix which level counts as positive, and
#'   [getCatVars()] / [getNumVars()] for the other partitions.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Two-level factors, which are the ones that can be modelled as 0/1
#' vars_Binary <- getBinaryVars(Labelled)
#' vars_Binary
#'
#' # `Revalued = FALSE` looks for any column with two distinct values instead
#' # of two factor levels, for frames that have not been through RevalueData().
#' getBinaryVars(SampleData, Revalued = FALSE)
#'
#' # Which level each one is scored against
#' createBinaryMapping(Labelled, vars_Binary)
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
getBinaryVars <- function(data,
    Ordinal = TRUE,
    Revalued = TRUE,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "getBinaryVars(DataFrame)", "getBinaryVars(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  if (Revalued) {
    # For revalued data, focus on factors with exactly 2 levels
    if (Ordinal) {
      # Include all factors (ordered and unordered)
      CatVars <- DataFrame %>%
        dplyr::select(dplyr::where(is.factor)) %>%
        names()
    } else {
      # Exclude ordered factors
      CatVars <- DataFrame %>%
        dplyr::select(dplyr::where(is.factor) & !dplyr::where(is.ordered)) %>%
        names()
    }

    # Select binary variables from categorical variables (2 levels)
    BinaryVars <- DataFrame %>%
      dplyr::select(dplyr::all_of(CatVars)) %>%
      dplyr::select(dplyr::where(~ length(levels(.)) == 2)) %>%
      names()
  } else {
    # For non-revalued data, check for variables with 2 unique values
    BinaryVars <- DataFrame %>%
      dplyr::select(dplyr::where(~ length(unique(.)) == 2)) %>%
      names()
  }

  # Return the names of binary variables
  return(BinaryVars)
}
