# Plot Phi Correlations Between Binary Variables

Computes pairwise phi coefficients (Φ) between binary categorical
variables with explicit 0/1 coding (1 == PositiveLevel from
[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md)),
then renders heatmap-style plots with raw and FDR-adjusted significance.

## Usage

``` r
PlotPhiHeatmap(Data, CatVars, Relabel = TRUE, binary_map = NULL)
```

## Arguments

- Data:

  A dataframe.

- CatVars:

  Character vector of binary categorical variable names.

- Relabel:

  Logical; if TRUE, uses sjlabelled variable labels for axes.

- binary_map:

  Optional mapping as returned by
  [`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md).
  If NULL, a mapping is created internally for `CatVars`.

## Value

A list with:

- `Unadjusted`: list(PvalTable, plot)

- `FDRCorrected`: list(PvalTable, plot)

- `method` = "Phi"

- `Relabel`

- `BinaryMapping` (used)
