# Create a Forest Plot from Univariate Regression Tables

This function generates a forest plot from the results of
[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md).

`plotForestFromTable()` was renamed to `PlotForestFromTable()` in
SciDataReportR 20.5.0 to match the package's `Plot*` naming convention.
It remains available as a backwards-compatible synonym.

## Usage

``` r
PlotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)

plotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)
```

## Arguments

- UnivariateRegressionTables:

  Either the full list returned by
  [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  (its `Results` dataframe is used directly), or a dataframe with the
  `Results` columns. Passing a dataframe lets you filter, reorder, or
  relabel `Results` before plotting; required columns are
  `OutcomeLabel`, `TermLabel`, `Estimate`, `ConfLow`, `ConfHigh`, and
  `PValue` (`Significant` and `ReferenceValue` are recomputed if
  absent). Lists created by older package versions (without a `Results`
  element) are still supported.

- pSize:

  Numeric. Size of the points in the plot. Default is 2.

- Flip:

  Logical. If `FALSE`, outcomes are facets and predictors/terms are
  rows. If `TRUE`, predictors/terms are facets and outcomes are rows.

## Value

A ggplot object representing the forest plot.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Build univariate regression tables to plot
urt <- MakeUnivariateRegressionTable(
  data = Labelled,
  outcome_vars = "AXL",
  predictor_vars = c("age", "Adiponectin", "Alpha_1_Antitrypsin")
)

PlotForestFromTable(urt)


# Or plot a filtered subset of the results
PlotForestFromTable(subset(urt$Results, Predictor != "age"))

# }
```
