#' Convert ordinal variables to numeric
#'
#' Convert ordinal variables in a dataframe to numeric scores.
#'
#' @param data The dataframe containing the variables.
#' @param variables A character vector specifying the names of variables to consider. If NULL, all columns of the dataframe will be considered.
#' @importFrom sjlabelled get_label set_label
#' @details Numeric scores recorded by [RevalueData()] from the codebook are
#' used when available. Otherwise the ordered-factor ranks are used.
#' @return The dataframe with ordinal variables converted to numeric scores.
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


  if (is.null(Variables)) Variables <- colnames(Data)
  missing_vars <- setdiff(Variables, names(Data))
  if (length(missing_vars)) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  ordinal_vars <- Variables[vapply(Data[Variables], ScidrIsOrdinal, logical(1))]
  for (var in ordinal_vars) {
    Data[[var]] <- ScidrOrdinalAsNumeric(Data[[var]])
  }
  Data
}
