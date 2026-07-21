# Create a Z-score plot with statistical significance

Compatibility alias for
[`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md).
Prefer
[`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md)
in new code because this function creates a scientific visualization.

## Usage

``` r
CreateZScorePlot(...)
```

## Arguments

- ...:

  Arguments passed to
  [`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md).

## Value

A ggplot object returned by
[`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md).

## See also

[`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md)
for the canonical function and full examples.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

CreateZScorePlot(
  Labelled,
  TargetVar = "Diagnosis",
  variables = c("age", "AXL", "Adiponectin")
)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```
