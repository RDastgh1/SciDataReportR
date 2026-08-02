# Plot Phi Correlations Between Binary Variables

Computes pairwise phi coefficients between binary categorical variables
with explicit 0/1 coding (1 == PositiveLevel from
[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md)),
then renders heatmap-style plots with raw and FDR-adjusted significance.

## Usage

``` r
PlotPhiHeatmap(
  data,
  CatVars,
  Relabel = TRUE,
  binary_map = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  Data = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe.

- CatVars:

  Character vector of binary categorical variable names.

- Relabel:

  Logical; if TRUE, uses sjlabelled variable labels for axes.

- binary_map:

  Optional mapping as returned by
  [`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md).
  If NULL, a mapping is created internally for `CatVars`.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all p-values at once (historical behavior).
  `"per_outcome"` corrects separately within each y-axis variable
  (`YVar`); the Phi matrix is symmetric, so this treats each variable's
  row of tiles as one family.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A list with:

- `Unadjusted`: list(PvalTable, plot)

- `FDRCorrected`: list(PvalTable, plot)

- `method` = "Phi"

- `Relabel`

- `BinaryMapping` (used)

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# CatVars must be binary (exactly two unique non-NA values)
result <- PlotPhiHeatmap(Labelled, CatVars = c("Diagnosis", "sex"))

# Raw p-value phi heatmap
result$Unadjusted$plot


# FDR-adjusted phi heatmap
result$FDRCorrected$plot
```
