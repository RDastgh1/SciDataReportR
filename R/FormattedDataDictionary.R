#' Create a formatted data dictionary table
#'
#' This function generates a formatted data dictionary table using the specified data frame.
#' The table includes variable names, labels, types, and additional formatting based on variable types.
#'
#' @param DataFrame The data frame for which the data dictionary is to be created.
#' @param numdecimals Number of decimals to display for numeric variables (default: 2).
#'
#' @return A formatted data dictionary table (gt object).
#'
#' @details This function requires the `gt` package. If not installed, the function will return an error.

#' @export
FormattedDataDictionary <- function(DataFrame, numdecimals = 2) {
  # Ensure `gt` is installed
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("The 'gt' package is required but not installed. Please install it using install.packages('gt').")
  }

  # Generate the initial data dictionary
  CB <- Make_DataDictionary(DataFrame, numdecimals)

  # Create the gt object
  g <- gt::gt(CB) %>%
    # Apply styles based on variable types
    gt::tab_style(style = gt::cell_fill(color = '#d2dfe6'), locations = gt::cells_body(rows = Type == "factor")) %>%
    gt::tab_style(style = gt::cell_fill(color = '#92b4c4'), locations = gt::cells_body(rows = Type == "factor" & `Ordered Factor` == TRUE)) %>%
    gt::tab_style(style = gt::cell_fill(color = '#c4929b'), locations = gt::cells_body(rows = Type == "character")) %>%
    gt::tab_style(style = gt::cell_fill(color = '#c4a292'), locations = gt::cells_body(rows = Type == "Date")) %>%

    # Apply style to column headers
    gt::tab_style(
      locations = gt::cells_column_labels(columns = gt::everything()),
      style = list(
        gt::cell_borders(sides = "bottom", weight = gt::px(3)),
        gt::cell_text(weight = "bold")
      )) %>%

    # Apply border style to cells
    gt::tab_style(
      style = gt::cell_borders(
        side = "left",
        color = "grey",
        weight = gt::px(0.1),
        style = "dashed"
      ),
      locations = list(
        gt::cells_body(gt::everything()),
        gt::cells_column_labels(gt::everything())
      )) %>%

    # Add footnotes with proper location
    gt::tab_footnote(
      footnote = paste("Data Frame has", nrow(DataFrame), "rows and", ncol(DataFrame), "columns"),
      locations = gt::cells_column_labels(columns = 1)
    ) %>%
    gt::tab_footnote(
      footnote = "Color indicates type of variable",
      locations = gt::cells_column_labels(columns = 1)
    )

  return(g)
}
