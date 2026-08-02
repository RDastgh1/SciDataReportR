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
#' @details This function requires the `gt` and `codebook` packages. If either
#' is not installed, the function will return an error.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach variable labels and factor levels first so the dictionary shows
#' # human-readable labels and correct variable types.
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Formatted data dictionary for a subset of variables
#' FormattedDataDictionary(
#'   Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin")]
#' )
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
