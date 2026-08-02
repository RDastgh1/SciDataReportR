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
  interactive = TRUE,
  Interactive = lifecycle::deprecated()
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

- interactive:

  Logical; if `TRUE`, plots are converted to interactive `plotly`
  objects using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
  Default is `TRUE`.

- Interactive:

  **Deprecated** (since 19.15.0). Use `interactive` instead.

## Value

If `Plot = "All"`, a named list of plots. Otherwise, a single plot
object. Plot objects are either `ggplot` objects or `plotly`
htmlwidgets, depending on `Interactive`.

## Details

The function expects the object returned by
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).
It does not re-run any merge validation checks.

## Examples

``` r
set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

diagnostics <- PlotMergeValidation(
  validation,
  Plot = "All",
  interactive = FALSE
)

# Merge-check status, key coverage, join audit, agreement, and conflicts
diagnostics$Checks

diagnostics$Coverage

diagnostics$JoinAudit

diagnostics$Agreement

diagnostics$Conflicts
```
