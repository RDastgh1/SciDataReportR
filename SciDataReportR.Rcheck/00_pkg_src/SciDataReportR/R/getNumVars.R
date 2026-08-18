#' Get Numeric Variables
#'
#' Extracts numeric variables from a data frame.
#'
#' @param data The data frame from which to extract numeric variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of numeric variables.
#' @details
#' The variable-selector helpers exist so that a variable set is derived from
#' the data rather than typed out and left to drift. `getNumVars()`,
#' [getCatVars()], and [getBinaryVars()] partition a labelled data frame the way
#' the analysis functions expect, and each takes an `Ordinal` argument because
#' an ordered factor can legitimately count as either continuous or categorical
#' depending on the analysis.
#'
#' Run them on a frame that has already been through [RevalueData()]. Before
#' relabelling, a 0/1 diagnosis is still numeric and would be picked up here as
#' continuous.
#'
#' @seealso [getCatVars()], [getBinaryVars()], and [ConvertOrdinalToNumeric()]
#'   for the ordinal policy these share.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # 128 numeric columns: age plus the biomarker panel
#' vars_Numeric <- getNumVars(Labelled)
#' length(vars_Numeric)
#' utils::head(vars_Numeric)
#'
#' # Ordinal variables are excluded by default; include them when they should
#' # be modelled on their numeric scale.
#' length(getNumVars(Labelled, Ordinal = TRUE))
#'
#' # The point of deriving the set: it feeds straight into an analysis without
#' # a hand-typed vector that can fall out of date.
#' MakeTable1(Labelled, variables = utils::head(vars_Numeric, 5))
#'
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

  NumVars <- DataFrame %>% dplyr::select(dplyr::where(is.numeric)) %>% names()
  if (Ordinal) {
    NumVars <- DataFrame %>% dplyr::select(dplyr::where(is.numeric) | dplyr::where(is.ordered)) %>% names()
  }
  return(NumVars)
}
