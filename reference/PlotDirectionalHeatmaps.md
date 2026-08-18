# Create directional heatmaps across continuous & binary variables

Combines:

- Continuous~Continuous (Pearson/Spearman)

- Binary~Binary (Phi; 1 == PositiveLevel)

- Binary~Continuous (r_pb; 1 == PositiveLevel) into a single square
  heatmap with raw and FDR-star overlays.

## Usage

``` r
PlotDirectionalHeatmaps(
  data,
  variables = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  Data = lifecycle::deprecated(),
  xVars = lifecycle::deprecated(),
  yVars = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe.

- variables:

  Character vector of variables to include (subset of `data` columns).
  The analysis is symmetric: every variable is related to every other,
  so a single variable set defines both axes. If NULL, uses all detected
  continuous + binary vars.

- Relabel:

  Logical; use sjlabelled variable labels if present.

- Ordinal:

  Logical; passed to
  [`PlotPointCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPointCorrelationsHeatmap.md)
  for the binary~continuous block, where it controls whether ordinal
  variables are treated as continuous. Defaults to `TRUE`.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, threaded through to
  the three sub-analyses
  ([`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md),
  [`PlotPhiHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPhiHeatmap.md),
  [`PlotPointCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPointCorrelationsHeatmap.md)).
  Correction is applied within each sub-analysis block
  (continuous~continuous, binary~binary, binary~continuous), matching
  historical behavior; each sub-function's documented outcome
  orientation applies within its block.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- xVars:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- yVars:

  **Deprecated** (since 19.15.0). Use `variables` instead. If supplied,
  the old rectangular x-by-y display is still honored.

## Value

list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)

## Details

Constant variables (no variation in the current data) are automatically
excluded before computing any tiles.

## Note

The analysis covers continuous and binary variables. Multi-level
categorical variables (more than two levels) are not placed on the
heatmap. `Ordinal` affects only the binary~continuous block, which is
the one sub-analysis that accepts it.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels so binary variables are detected and
# axis labels are readable
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A mix of binary categorical (Diagnosis, sex) and continuous variables
result <- PlotDirectionalHeatmaps(
  Labelled,
  variables = c("Diagnosis", "sex", "age", "AXL", "Adiponectin",
                "Alpha_1_Antitrypsin", "C_Reactive_Protein", "Cortisol",
                "Insulin", "Leptin")
)

# Raw p-value directional heatmap
result$Unadjusted$plot
#> Warning: Removed 24 rows containing missing values or values outside the scale range
#> (`geom_text()`).


# FDR-adjusted directional heatmap
result$FDRCorrected$plot
#> Warning: Removed 26 rows containing missing values or values outside the scale range
#> (`geom_text()`).


# How binary variables were coded (which level counts as the positive one)
result$BinaryMapping
#>    Variable     Label PositiveLevel NegativeLevel
#> 1 Diagnosis Diagnosis      Impaired       Control
#> 2       sex       Sex          Male        Female
```
