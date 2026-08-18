#To do: add warnings if the strings were not in the original

#' Remove Strings from a Vector
#'
#' This function removes strings from a vector that are present in another vector.
#'
#' @param Orig The original vector containing strings.
#' @param Remove The vector of strings to be removed from the original vector.
#' @return A vector containing the strings from the original vector that were not present in the removal vector.
#' @details
#' A set-difference helper for variable vectors, equivalent to
#' `setdiff(Orig, Remove)` but keeping duplicates in `Orig`. It exists for the
#' common case of taking an auto-detected variable set and dropping the few
#' columns that should not be analyzed - identifiers, dates, the outcome
#' itself - without retyping the ones that should.
#'
#' Names in `Remove` that do not appear in `Orig` are ignored silently, so a
#' stale exclusion list will not error; it will simply stop having an effect.
#'
#' @seealso [getNumVars()], [getCatVars()], and [getBinaryVars()], which
#'   produce the vectors this usually trims.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Every numeric column, including ones that are not really measurements
#' vars_Numeric <- getNumVars(Labelled)
#' length(vars_Numeric)
#'
#' # Drop the ones that should not be screened as biomarkers
#' vars_Biomarkers <- removeString(vars_Numeric, c("age", "tau", "p_tau"))
#' length(vars_Biomarkers)
#' utils::head(vars_Biomarkers, 4)
#'
#' # Names that are not present are ignored, so an over-broad exclusion list
#' # is harmless
#' removeString(c("age", "sex", "AXL"), c("sex", "not_a_column"))
#'
#' # Unlike setdiff(), duplicates in the original are kept
#' removeString(c("a", "a", "b", "c"), "b")
#' setdiff(c("a", "a", "b", "c"), "b")
#'
#' @export
removeString <- function(Orig, Remove) {
  leftovers <- Orig[Orig %in% Remove == FALSE]
  return(leftovers)
}
