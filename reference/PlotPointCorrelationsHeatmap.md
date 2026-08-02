# Plot Point-Biserial Correlations Between Binary and Continuous Variables

Calculates point-biserial correlations (binary vs continuous) with
explicit 0/1 coding where 1 == PositiveLevel from
[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md),
and renders heatmap-style tiles.

## Usage

``` r
PlotPointCorrelationsHeatmap(
  data,
  CatVars,
  ContVars,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  binary_map = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  Data = lifecycle::deprecated(),
  Covariates = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe.

- CatVars:

  Character vector of binary categorical variables.

- ContVars:

  Character vector of continuous variables.

- covariates:

  Optional covariates (reserved).

- Relabel:

  Logical; use sjlabelled variable labels for axes.

- Ordinal:

  Logical; reserved for future use.

- binary_map:

  Optional mapping as returned by
  [`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md).
  If NULL, a mapping is created internally for `CatVars`.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all p-values at once (historical behavior).
  `"per_outcome"` corrects separately within each continuous variable:
  outcomes are the continuous variables (`ContVars`).

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A list with Unadjusted, FDRCorrected, method ("R_pb"), Relabel,
Covariates, BinaryMapping.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# CatVars must be binary (exactly two unique non-NA values)
result <- PlotPointCorrelationsHeatmap(
  Labelled,
  CatVars = c("Diagnosis", "sex"),
  ContVars = c("age", "AXL", "Adiponectin")
)

# Raw p-value point-biserial heatmap
result$Unadjusted$plot
#> Warning: Removed 5 rows containing missing values or values outside the scale range
#> (`geom_text()`).


# FDR-adjusted point-biserial heatmap
result$FDRCorrected$plot
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_text()`).
```
