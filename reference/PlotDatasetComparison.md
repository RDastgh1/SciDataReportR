# Plot dataset comparison diagnostics

Create diagnostic plots from a
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
result object. This function visualizes dataset-version changes,
including check status, summary metrics, structure changes, and
variable-level value changes.

## Usage

``` r
PlotDatasetComparison(
  CompareObj,
  Plot = c("All", "Checks", "SummaryMetrics", "StructureChanges", "VariableChanges",
    "TopChangedVariables"),
  Interactive = TRUE,
  TopN = 10
)
```

## Arguments

- CompareObj:

  A list returned by
  [`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md).

- Plot:

  Character value specifying which plot to return. Options are `"All"`,
  `"Checks"`, `"SummaryMetrics"`, `"StructureChanges"`,
  `"VariableChanges"`, and `"TopChangedVariables"`. Default is `"All"`.

- Interactive:

  Logical; if `TRUE`, plots are converted to interactive `plotly`
  objects using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
  Default is `TRUE`.

- TopN:

  Integer number of variables or records to preview in plots and hover
  text. Default is `10`.

## Value

If `Plot = "All"`, a named list of plots. Otherwise, a single plot
object. Plot objects are either `ggplot` objects or `plotly`
htmlwidgets, depending on `Interactive`.

## Details

Use this function after running
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
to quickly inspect what changed between two versions of a dataset.
