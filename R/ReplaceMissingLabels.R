#' Replace Missing Labels in Dataframe Columns
#'
#' This function iterates through the columns of a dataframe and assigns the column name as the label to any column that does not have a label.
#'
#' @param df A dataframe.
#' @return The input dataframe with missing labels replaced.
#' @importFrom sjlabelled get_label set_label
#' @importFrom labelled var_label
#' @export

ReplaceMissingLabels <- function(df) {
  # Find columns without labels
  cols_without_labels <- names(df)[sapply(df, function(x) is.null(sjlabelled::get_label(x)))]

  # Assign column name as label
  for (col_name in cols_without_labels) {
    labelled::var_label(df[[col_name]]) <- col_name
  }

  return(df)
}
