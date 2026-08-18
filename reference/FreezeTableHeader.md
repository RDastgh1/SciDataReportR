# Freeze the header row of a long table when scrolling

Render a table with a "sticky" header row that stays visible while the
reader scrolls, so the column meanings are never lost in a long table.
This is handy for the tall tables produced by
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
and other `gtsummary` helpers when they are knit into HTML Quarto or R
Markdown documents.

## Usage

``` r
FreezeTableHeader(
  x,
  height = NULL,
  width = NULL,
  header_background = "white",
  bootstrap_options = c("striped", "hover", "condensed"),
  full_width = FALSE,
  font_size = NULL,
  ...
)
```

## Arguments

- x:

  A `gtsummary` object (such as the output of
  [`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)),
  a data frame or tibble, or an existing kableExtra / `knitr_kable`
  object.

- height:

  Optional CSS height for a scroll box, for example `"400px"` or
  `"60vh"`. When supplied, the table scrolls within a box of this height
  with the header frozen. When `NULL` (default), the header sticks
  during page-level scrolling instead.

- width:

  Optional CSS width for the scroll box, for example `"100%"`. Only used
  when `height` is supplied.

- header_background:

  Background color for the frozen header row. Defaults to `"white"` so
  the header cleanly covers rows scrolling underneath it.

- bootstrap_options:

  Character vector of kableExtra bootstrap styling options passed to
  [`kableExtra::kable_styling()`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).
  Defaults to `c("striped", "hover", "condensed")`.

- full_width:

  Logical; passed to
  [`kableExtra::kable_styling()`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).
  Default `FALSE`.

- font_size:

  Optional numeric font size passed to
  [`kableExtra::kable_styling()`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).

- ...:

  Additional arguments passed to
  [`kableExtra::kable_styling()`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).

## Value

A kableExtra HTML table object with a frozen header, suitable for
printing in a Quarto or R Markdown chunk.

## Details

The table is rendered through kableExtra with
`kable_styling(fixed_thead = ...)`, which pins the header row using CSS
`position: sticky`. Two scrolling modes are supported:

- With `height = NULL` (the default) the header sticks to the top of the
  viewport as the whole page scrolls past the table.

- With a `height` (for example `"400px"`) the table is placed in its own
  scroll box of that height and the header stays frozen at the top of
  the box.

Sticky headers only take effect in HTML output. In PDF or Word output
the table renders normally without a frozen header.

## See also

[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Repeat the full comparison table to make scrolling and the frozen header
# obvious on the reference page.
tbl_Base <- MakeComparisonTable(
  data = Labelled,
  group_var = "Diagnosis",
  variables = setdiff(names(Labelled), "Diagnosis")
)

tbl_All <- dplyr::bind_rows(rep(list(tbl_Base$table_body), 8))

# Inside a fixed-height vertical scroll box
htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(tbl_All, height = "450px", width = "100%")
)))


 variable 
```
