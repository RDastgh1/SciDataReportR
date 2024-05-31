#' Create a formatted data dictionary table
#'
#' This function generates a formatted data dictionary table using the specified data frame.
#' The table includes variable names, labels, types, and additional formatting based on variable types.
#'
#' @param DataFrame The data frame for which the data dictionary is to be created.
#' @param numdecimals Number of decimals to display for numeric variables (default: 2).
#'
#' @return A formatted data dictionary table.
#'
#' @importFrom gt gt tab_style tab_footnote cells_body cells_column_labels everything cell_fill cell_borders px
#'
#' @export
FormattedDataDictionary <- function(DataFrame, numdecimals = 2) {

  # Generate the initial data dictionary
  CB <- Make_DataDictionary(DataFrame, numdecimals)

  # Create the gt object
  g <- gt::gt(CB)

  # Apply styles based on variable types
  g <- g %>%
    tab_style(style = cell_fill(color = '#d2dfe6'), locations = cells_body(rows = Type == "factor")) %>%
    tab_style(style = cell_fill(color = '#92b4c4'), locations = cells_body(rows = Type == "factor" & `Ordered Factor` == TRUE)) %>%
    tab_style(style = cell_fill(color = '#c4929b'), locations = cells_body(rows = Type == "character")) %>%
    tab_style(style = cell_fill(color = '#c4a292'), locations = cells_body(rows = Type == "Date")) %>%

    # Apply style to column headers
    tab_style(
      locations = cells_column_labels(columns = everything()),
      style = list(
        cell_borders(sides = "bottom", weight = gt::px(3)),
        cell_text(weight = "bold")
      )) %>%

    # Apply border style to cells
    tab_style(
      style = cell_borders(
        side = "left",
        color = "grey",
        weight = gt::px(0.1),
        style = "dashed"
      ),
      locations = list(
        cells_body(
          everything()
        ),
        cells_column_labels(
          everything()
        )
      )) %>%

    # Add footnotes
    tab_footnote(
      footnote = paste("Data Frame had", nrow(DataFrame), "rows and", ncol(DataFrame), "columns")
    ) %>%
    tab_footnote(
      footnote = "Color indicates type of variable"
    )

  return(g)
}
