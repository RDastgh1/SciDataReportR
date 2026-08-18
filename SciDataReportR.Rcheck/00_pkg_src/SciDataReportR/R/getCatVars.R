#' Get Categorical Variables
#'
#' Extracts categorical variables from a data frame.
#'
#' @param data The data frame from which to extract categorical variables.
#' @param Ordinal Logical, indicating whether to include ordinal variables.
#' @return A character vector containing the names of categorical variables.
#' @seealso [getNumVars()] and [getBinaryVars()] for the other partitions.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Every factor in the frame
#' getCatVars(Labelled)
#'
#' # `Ordinal = FALSE` drops ordered factors, which is what you want when the
#' # ordered variables are going to be analyzed on their numeric scale instead.
#' getCatVars(Labelled, Ordinal = FALSE)
#'
#' # Only meaningful after RevalueData(): in the raw extract `sex` is still a
#' # bare 0/1 numeric column and is not detected as categorical.
#' getCatVars(SampleData)
#'
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
      dplyr::select(dplyr::where(is.factor)) %>%
      names()
  } else {
    CatVars <- DataFrame %>%
      dplyr::select(dplyr::where(is.factor) & !dplyr::where(is.ordered)) %>%
      names()
  }
  return(CatVars)
}
