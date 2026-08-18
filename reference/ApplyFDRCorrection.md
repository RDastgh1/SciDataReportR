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
  predictor_ids = NULL,
  symmetric = "auto",
  include_diagonal = FALSE
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

- symmetric:

  How to treat a square matrix whose two triangles hold the same
  p-values. `"auto"` (default) detects symmetry and corrects each pair
  once; `TRUE` requires it (and errors if `pmat` is not symmetric);
  `FALSE` restores the whole-matrix behavior that counts each pair
  twice. Ignored for vector input. See the Symmetric matrices section.

- include_diagonal:

  Logical; for symmetric input, whether the diagonal (self-comparisons)
  joins the family being corrected. Default `FALSE`, which excludes it
  and returns `NA` on the diagonal.

## Value

An object of the same shape as `pmat` (matrix, data frame, or vector)
containing adjusted p-values. Non-finite entries remain `NA`.

## Symmetric matrices

When a variable set is correlated against itself, the resulting p-value
matrix is symmetric: every pair appears twice, above and below the
diagonal, and the diagonal holds self-comparisons that were never
tested. Correcting across every cell therefore describes a family of
`n * (n - 1)` tests when only `n * (n - 1) / 2` were actually run - 90
instead of 45 for ten variables.

How much that matters depends on the method, and on what sits on the
diagonal:

- Benjamini-Hochberg (`method = "fdr"`, the default) is invariant to
  exact duplication, because doubling both a p-value's rank and the
  family size cancels out. Adjusted values are unchanged.

- `"bonferroni"`, `"holm"`, and `"hochberg"` are not invariant. Counting
  each pair twice made every adjusted p-value exactly twice as large as
  it should be.

- A diagonal carrying real values is the damaging case in every method.
  A self-correlation p-value of 0 enters the family as the most
  significant test there is, which drags the whole ranking down and
  makes the off-diagonal results look *stronger* than they are.

`symmetric = "auto"` (the default) detects symmetry and corrects the
unique pairs once, then mirrors the adjusted values back into the lower
triangle so the matrix stays symmetric. The diagonal is excluded from
the family and returned as `NA` unless `include_diagonal = TRUE`. Set
`symmetric = FALSE` to force the old whole-matrix behavior, or
`symmetric = TRUE` to require symmetric handling and error out if `pmat`
is not symmetric.

Detection ignores the diagonal, since the heatmap functions routinely
blank it before correcting, and requires the two triangles to agree on
both their values and their missing cells. An asymmetric matrix is never
affected.

Under `"per_outcome"` or `"per_predictor"` scope each row or column is
already a family of unique tests, so symmetric input only has its
diagonal excluded. The result of a per-group correction on a symmetric
matrix is not itself symmetric.

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

# A symmetric matrix: six filled cells, but only three pairs tested
vars <- c("mmse", "trails", "digit_span")
pm_sym <- matrix(NA_real_, 3, 3, dimnames = list(vars, vars))
pm_sym[lower.tri(pm_sym)] <- c(0.001, 0.020, 0.040)
pm_sym[upper.tri(pm_sym)] <- t(pm_sym)[upper.tri(pm_sym)]
pm_sym
#>             mmse trails digit_span
#> mmse          NA  0.001       0.02
#> trails     0.001     NA       0.04
#> digit_span 0.020  0.040         NA

# Detected automatically: three tests in the family, diagonal excluded
ApplyFDRCorrection(pm_sym)
#>             mmse trails digit_span
#> mmse          NA  0.003       0.03
#> trails     0.003     NA       0.04
#> digit_span 0.030  0.040         NA

# Bonferroni, corrected once per pair and then against all six cells
ApplyFDRCorrection(pm_sym, method = "bonferroni")
#>             mmse trails digit_span
#> mmse          NA  0.003       0.06
#> trails     0.003     NA       0.12
#> digit_span 0.060  0.120         NA
ApplyFDRCorrection(pm_sym, method = "bonferroni", symmetric = FALSE)
#>             mmse trails digit_span
#> mmse          NA  0.006       0.12
#> trails     0.006     NA       0.24
#> digit_span 0.120  0.240         NA
```
