#' Replace Missing Labels in Dataframe Columns
#'
#' This function iterates through the columns of a dataframe and assigns the column name as the label to any column that does not have a label.
#'
#' @param data A dataframe.
#' @return The input dataframe with missing labels replaced.
#' @importFrom sjlabelled get_label set_label
#' @importFrom labelled var_label
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
