# Apply multiple-comparison correction across a p-value matrix

Central helper used by the SciDataReportR heatmap/matrix family to apply
multiple-comparison correction (FDR by default) with a selectable scope.
All non-finite p-values (`NA`, `NaN`) are left untouched and excluded
from the correction, matching the long-standing behavior of the plotting
functions that now delegate to this helper.

## Usage

``` r
ApplyFDRCorrection(
  pmat,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  outcome_margin = 2,
  method = "fdr",
  outcome_ids = NULL,
  predictor_ids = NULL
)
```

## Arguments

- pmat:

  A numeric matrix or data frame of p-values, or a plain numeric vector.
  Matrices and data frames keep their dimensions and dimnames.

- fdr_scope:

  Either `"matrix"` (default), `"per_outcome"`, or `"per_predictor"`.
  `"matrix"` corrects across all p-values at once (one family).
  `"per_outcome"` corrects separately within each outcome's p-values:
  for matrix input, groups run along `outcome_margin`; for vector input,
  groups are defined by `outcome_ids`.

- outcome_margin:

  For matrix/data-frame input with `fdr_scope = "per_outcome"`: `2`
  (default) if outcomes are columns, `1` if outcomes are rows. Ignored
  for `"matrix"` scope and for vector input.

- method:

  Correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"fdr"` (Benjamini-Hochberg).

- outcome_ids:

  Optional vector (same length as `pmat`) identifying the outcome each
  p-value belongs to. Only used - and then required - when `pmat` is a
  vector and `fdr_scope = "per_outcome"`. This is how the long-format
  table functions (for example
  [`PlotPhiHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPhiHeatmap.md)
  or
  [`PlotChiSqCovar()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotChiSqCovar.md))
  group their p-values by outcome.

- predictor_ids:

  Optional vector (same length as `pmat`) identifying the predictor each
  p-value belongs to. Only used - and then required - when `pmat` is a
  vector and `fdr_scope = "per_predictor"`.

## Value

An object of the same shape as `pmat` (matrix, data frame, or vector)
containing adjusted p-values. Non-finite entries remain `NA`.

## Examples

``` r
pm <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
             nrow = 2,
             dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3")))

# One family across the whole matrix (classic behavior)
ApplyFDRCorrection(pm)
#>       out1 out2 out3
#> pred1 0.06 0.06 0.06
#> pred2 0.06 0.06 0.06

# Correct each outcome (column) separately
ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 2)
#>       out1 out2 out3
#> pred1 0.02 0.04 0.06
#> pred2 0.02 0.04 0.06

# Vector input with explicit outcome grouping
ApplyFDRCorrection(c(0.01, 0.04, 0.02, 0.03),
                   fdr_scope = "per_outcome",
                   outcome_ids = c("y1", "y1", "y2", "y2"))
#> [1] 0.02 0.04 0.03 0.03
```
