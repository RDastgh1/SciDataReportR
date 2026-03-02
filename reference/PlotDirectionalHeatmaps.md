# Create directional heatmaps across continuous & binary variables

Combines:

- Continuous~Continuous (Pearson/Spearman)

- Binary~Binary (Phi; 1 == PositiveLevel)

- Binary~Continuous (r_pb; 1 == PositiveLevel) into a single square
  heatmap with raw and FDR-star overlays.

## Usage

``` r
PlotDirectionalHeatmaps(
  Data,
  xVars = NULL,
  yVars = NULL,
  Relabel = TRUE,
  Ordinal = TRUE
)
```

## Arguments

- Data:

  A dataframe.

- xVars:

  Character vector of variables for x-axis (subset of Data cols). If
  NULL, uses all detected continuous + binary vars.

- yVars:

  Character vector for y-axis (defaults to xVars if NULL).

- Relabel:

  Logical; use sjlabelled variable labels if present.

- Ordinal:

  Logical; reserved for future use.

## Value

list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)

## Details

Constant variables (no variation in the current data) are automatically
excluded before computing any tiles.
