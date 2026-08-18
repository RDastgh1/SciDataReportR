#' Freeze the header row of a long table when scrolling
#'
#' Render a table with a "sticky" header row that stays visible while the reader
#' scrolls, so the column meanings are never lost in a long table. This is handy
#' for the tall tables produced by [MakeComparisonTable()] and other
#' `gtsummary` helpers when they are knit into HTML Quarto or R Markdown
#' documents.
#'
#' @details
#' The table is rendered through \pkg{kableExtra} with
#' `kable_styling(fixed_thead = ...)`, which pins the header row using CSS
#' `position: sticky`. Two scrolling modes are supported:
#'
#' * With `height = NULL` (the default) the header sticks to the top of the
#'   viewport as the whole page scrolls past the table.
#' * With a `height` (for example `"400px"`) the table is placed in its own
#'   scroll box of that height and the header stays frozen at the top of the box.
#'
#' Sticky headers only take effect in HTML output. In PDF or Word output the
#' table renders normally without a frozen header.
#'
#' @param x A `gtsummary` object (such as the output of [MakeComparisonTable()]),
#'   a data frame or tibble, or an existing \pkg{kableExtra} / `knitr_kable`
#'   object.
#' @param height Optional CSS height for a scroll box, for example `"400px"` or
#'   `"60vh"`. When supplied, the table scrolls within a box of this height with
#'   the header frozen. When `NULL` (default), the header sticks during
#'   page-level scrolling instead.
#' @param width Optional CSS width for the scroll box, for example `"100%"`.
#'   Only used when `height` is supplied.
#' @param header_background Background color for the frozen header row. Defaults
#'   to `"white"` so the header cleanly covers rows scrolling underneath it.
#' @param bootstrap_options Character vector of \pkg{kableExtra} bootstrap
#'   styling options passed to `kableExtra::kable_styling()`. Defaults to
#'   `c("striped", "hover", "condensed")`.
#' @param full_width Logical; passed to `kableExtra::kable_styling()`. Default
#'   `FALSE`.
#' @param font_size Optional numeric font size passed to
#'   `kableExtra::kable_styling()`.
#' @param ... Additional arguments passed to `kableExtra::kable_styling()`.
#'
#' @return A \pkg{kableExtra} HTML table object with a frozen header, suitable
#'   for printing in a Quarto or R Markdown chunk.
#'
#' @seealso [MakeComparisonTable()]
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Repeat the full comparison table to make scrolling and the frozen header
#' # obvious on the reference page.
#' tbl_Base <- MakeComparisonTable(
#'   data = Labelled,
#'   group_var = "Diagnosis",
#'   variables = setdiff(names(Labelled), "Diagnosis")
#' )
#'
#' tbl_All <- dplyr::bind_rows(rep(list(tbl_Base$table_body), 8))
#'
#' # Inside a fixed-height vertical scroll box
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(tbl_All, height = "450px", width = "100%")
#' )))
#' }
#'
#' @export
FreezeTableHeader <- function(x,
    height = NULL,
    width = NULL,
    header_background = "white",
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size = NULL,
    ...) {

  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop(
      "FreezeTableHeader() requires the 'kableExtra' package. ",
      "Install it with install.packages('kableExtra')."
    )
  }

  # Convert the supported input types into a knitr_kable / kableExtra object.
  if (inherits(x, "kableExtra") || inherits(x, "knitr_kable")) {
    kbl <- x
  } else if (inherits(x, "gtsummary")) {
    if (!requireNamespace("gtsummary", quietly = TRUE)) {
      stop("FreezeTableHeader() requires the 'gtsummary' package for gtsummary input.")
    }
    kbl <- gtsummary::as_kable_extra(x)
  } else if (inherits(x, "gt_tbl")) {
    stop(
      "FreezeTableHeader() cannot freeze 'gt' objects directly. ",
      "Pass the gtsummary object (before as_gt()), a data frame, or a ",
      "kableExtra object instead."
    )
  } else if (is.data.frame(x)) {
    kbl <- knitr::kable(x, format = "html")
  } else {
    stop(
      "x must be a gtsummary object, a data frame, or a kableExtra/knitr_kable ",
      "object. Received an object of class: ",
      paste(class(x), collapse = ", ")
    )
  }

  styled <- kableExtra::kable_styling(
    kbl,
    bootstrap_options = bootstrap_options,
    full_width = full_width,
    font_size = font_size,
    fixed_thead = list(enabled = TRUE, background = header_background),
    ...
  )

  if (!is.null(height) || !is.null(width)) {
    styled <- kableExtra::scroll_box(styled, height = height, width = width)
  }

  styled
}
