#' Create a formatted data dictionary table
#'
#' This function generates a formatted data dictionary table using the specified data frame.
#' The table includes variable names, labels, types, and additional formatting based on variable types.
#'
#' @param data The data frame for which the data dictionary is to be created.
#' @param digits Number of decimals to display for numeric variables (default: 2).
#'
#' @return A formatted data dictionary table (gt object).
#'
#' @section Giving the table room:
#' A data dictionary is wide by nature - name, label, type, levels, and summary
#' statistics all belong on one row - so it reads badly squeezed into a narrow
#' text column. In a Quarto report the chunk option `column: screen` handles
#' this. The equivalent on a reference page is to let the table spill into the
#' page margins and scroll horizontally for the rest, which the example's
#' `ColumnScreen()` helper does.
#'
#' The `max()` clamp in that helper is what keeps it safe: the table breaks out
#' as far as the margins allow but never past the edge of the window, so no
#' column is pushed off screen and the page itself never scrolls sideways.
#'
#' @details This function requires the `gt` and `codebook` packages. If either
#' is not installed, the function will return an error.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels first
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Formatted data dictionary for a subset of variables
#' dictionary <- FormattedDataDictionary(
#'   Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin")]
#' )
#'
#' # Break the table out to the full width of the window
#' ColumnScreen <- function(x) {
#'   htmltools::browsable(htmltools::div(
#'     style = paste(
#'       "position: relative;",
#'       "margin-left: max(-5rem, calc(50% - 50vw));",
#'       "margin-right: max(-5rem, calc(50% - 50vw));",
#'       "box-sizing: border-box;",
#'       "overflow-x: auto;"
#'     ),
#'     htmltools::HTML(gt::as_raw_html(x))
#'   ))
#' }
#'
#' ColumnScreen(dictionary)
#' }
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param numdecimals \strong{Deprecated} (since 19.15.0). Use \code{digits} instead.
#' @export
FormattedDataDictionary <- function(data,
    digits = 2,
    DataFrame = lifecycle::deprecated(),
    numdecimals = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "FormattedDataDictionary(DataFrame)", "FormattedDataDictionary(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(numdecimals)) {
    lifecycle::deprecate_warn("19.15.0", "FormattedDataDictionary(numdecimals)", "FormattedDataDictionary(digits)")
    digits <- numdecimals
  }
  numdecimals <- digits

  # Ensure `gt` is installed
  if (!ScidrPackageAvailable("gt")) {
    stop("The 'gt' package is required but not installed. Please install it using install.packages('gt').")
  }
  if (!ScidrPackageAvailable("codebook")) {
    stop("The 'codebook' package is required but not installed. Please install it using install.packages('codebook').")
  }

  # Generate the initial data dictionary
  CB <- MakeDataDictionary(DataFrame, numdecimals)

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

ScidrPackageAvailable <- function(package) {
  requireNamespace(package, quietly = TRUE)
}
