# Plot a diagnostic likelihood-ratio heatmap

Visualizes results created by
[`DiagnosticLikelihoodRatioTable()`](https://rdastgh1.github.io/SciDataReportR/reference/DiagnosticLikelihoodRatioTable.md)
without recomputing diagnostic statistics. Tile fill is log2(LR), making
reciprocal evidence (for example, LR 0.25 and LR 4) equally distant from
LR 1.

## Usage

``` r
PlotDiagnosticLRHeatmap(
  x,
  result = c("all", "positive", "negative"),
  predictor_order = c("original", "alphabetical", "strength", "cluster"),
  outcome_order = c("original", "alphabetical", "strength", "cluster"),
  orientation = c("auto", "predictors_rows", "outcomes_rows"),
  show_values = "auto",
  show_ci_marker = TRUE,
  facet_strata = TRUE,
  cap = NULL,
  na_color = "grey90"
)
```

## Arguments

- x:

  An object returned by
  [`DiagnosticLikelihoodRatioTable()`](https://rdastgh1.github.io/SciDataReportR/reference/DiagnosticLikelihoodRatioTable.md),
  or a tidy data frame compatible with its `Results` element.

- result:

  Diagnostic result levels to display: `"all"`, `"positive"`, or
  `"negative"`. Positive and negative selections require the result
  object.

- predictor_order:

  Predictor row ordering: `"original"`, `"alphabetical"`, `"strength"`,
  or `"cluster"`.

- outcome_order:

  Outcome column ordering with the same choices.

- orientation:

  Tile orientation: `"auto"`, `"predictors_rows"`, or `"outcomes_rows"`.

- show_values:

  `"auto"`, `TRUE`, or `FALSE`; auto labels at most 150 tiles.

- show_ci_marker:

  Logical; append `*` when the unadjusted LR CI excludes 1.

- facet_strata:

  Logical; facet separate diagnostic strata when present.

- cap:

  Optional positive maximum absolute log2(LR) for color scaling.

- na_color:

  Fill color for unavailable LR values.

## Value

A static `ggplot`. Its `DiagnosticLRData` attribute and tile `text`
aesthetic contain complete hover-ready information for optional Plotly
use.

## See also

[`DiagnosticLikelihoodRatioTable()`](https://rdastgh1.github.io/SciDataReportR/reference/DiagnosticLikelihoodRatioTable.md)
to calculate the displayed LRs.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)
df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
df_Labelled$DiagnosisBinary <- factor(df_Labelled$Diagnosis,
  levels = c("Control", "Impaired"))
lr <- DiagnosticLikelihoodRatioTable(df_Labelled, "DiagnosisBinary",
  c("sex", "Genotype"))
#> DiagnosisBinary: 'Impaired' treated as outcome-positive.
PlotDiagnosticLRHeatmap(lr, result = "all")
```
