#' Replace Missing Labels in Dataframe Columns
#'
#' This function iterates through the columns of a dataframe and assigns the column name as the label to any column that does not have a label.
#'
#' @param data A dataframe.
#' @return The input dataframe with missing labels replaced.
#' @importFrom sjlabelled get_label set_label
#' @importFrom labelled var_label
#' @examples
#' # A data frame where only some columns carry a variable label
#' df <- data.frame(
#'   age    = c(52, 61, 77),
#'   bmi    = c(24.1, 29.7, 22.0),
#'   smoker = c(0, 1, 0)
#' )
#' labelled::var_label(df$age) <- "Age (years)"
#'
#' # Before: bmi and smoker have no label
#' sapply(df, function(x) sjlabelled::get_label(x))
#'
#' # After: unlabelled columns are labelled with their column name
#' filled <- ReplaceMissingLabels(df)
#' sapply(filled, function(x) sjlabelled::get_label(x))
#'
#' @param df \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export

ReplaceMissingLabels <- function(data,
    df = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(df)) {
    lifecycle::deprecate_warn("19.15.0", "ReplaceMissingLabels(df)", "ReplaceMissingLabels(data)")
    data <- df
  }
  if (!missing(data)) df <- data

  # Find columns without labels
  cols_without_labels <- names(df)[sapply(df, function(x) is.null(sjlabelled::get_label(x)))]

  # Assign column name as label
  for (col_name in cols_without_labels) {
    labelled::var_label(df[[col_name]]) <- col_name
  }

  return(df)
}
