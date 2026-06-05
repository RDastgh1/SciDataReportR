# Plot merge validation diagnostics

Create diagnostic plots from a
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
result object. This function visualizes key merge-audit outputs,
including validation check status, key coverage, join-variable auditing,
duplicate-variable agreement, and duplicate-variable conflict counts.
Use this after running
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
to quickly inspect whether a merged dataset appears trustworthy.

## Usage

``` r
PlotMergeValidation(
  MergeObj,
  Plot = c("All", "Checks", "Coverage", "JoinAudit", "Agreement", "Conflicts"),
  Interactive = TRUE
)
```

## Arguments

- MergeObj:

  A list returned by
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).

- Plot:

  Character value specifying which plot to return. Options are `"All"`,
  `"Checks"`, `"Coverage"`, `"JoinAudit"`, `"Agreement"`, and
  `"Conflicts"`. Default is `"All"`.

- Interactive:

  Logical; if `TRUE`, plots are converted to interactive `plotly`
  objects using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
  Default is `TRUE`.

## Value

If `Plot = "All"`, a named list of plots. Otherwise, a single plot
object. Plot objects are either `ggplot` objects or `plotly`
htmlwidgets, depending on `Interactive`.

## Details

The function expects the object returned by
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).
It does not re-run any merge validation checks.
