# Plot Point-Biserial Correlations Between Binary and Continuous Variables

Calculates point-biserial correlations (binary vs continuous) with
explicit 0/1 coding where 1 == PositiveLevel from
[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md),
and renders heatmap-style tiles.

## Usage

``` r
PlotPointCorrelationsHeatmap(
  Data,
  CatVars,
  ContVars,
  Covariates = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  binary_map = NULL
)
```

## Arguments

- Data:

  A dataframe.

- CatVars:

  Character vector of binary categorical variables.

- ContVars:

  Character vector of continuous variables.

- Covariates:

  Optional covariates (reserved).

- Relabel:

  Logical; use sjlabelled variable labels for axes.

- Ordinal:

  Logical; reserved for future use.

- binary_map:

  Optional mapping as returned by
  [`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md).
  If NULL, a mapping is created internally for `CatVars`.

## Value

A list with Unadjusted, FDRCorrected, method ("R_pb"), Relabel,
Covariates, BinaryMapping.
