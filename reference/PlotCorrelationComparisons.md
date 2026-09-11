# Compare correlations between two independent groups

Computes correlations or partial correlations separately within two
independent groups, compares corresponding correlations, and visualizes
the between-group difference as a heatmap. Tile color represents
`DeltaR = r_comparison - r_reference`, significance stars represent the
statistical test comparing the two correlations, and striped tiles
indicate correlations with opposite signs between groups.

## Usage

``` r
PlotCorrelationComparisons(
  data,
  predictor_vars = NULL,
  outcome_vars = NULL,
  group_var,
  comparison_group = NULL,
  reference_group = NULL,
  covariates = NULL,
  method = "pearson",
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  min_n = 4,
  eps = 1e-12,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  reversal_style = c("outline", "stripe", "none"),
  interactive = c("none", "plotly", "girafe", "both"),
  low_color = "#B2182B",
  mid_color = "white",
  high_color = "#2166AC",
  color_limits = c(-2, 2),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  triangle = c("full", "upper", "lower")
)
```

## Arguments

- data:

  A data frame.

- predictor_vars:

  Character vector of predictor variables. If `NULL`, variable selection
  is inherited from
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).

- outcome_vars:

  Optional character vector of outcome variables. If `NULL`, the same
  variables are used on both axes.

- group_var:

  Character string naming the grouping variable.

- comparison_group:

  Optional level of `group_var` used as the comparison group. Positive
  DeltaR values indicate a more positive correlation in this group
  relative to the reference group.

- reference_group:

  Optional level of `group_var` used as the reference group.

- covariates:

  Optional character vector of covariates used to calculate partial
  correlations within each group.

- method:

  Correlation method. Either `"pearson"` or `"spearman"`.

- Relabel:

  Logical indicating whether variable labels should be used when
  available.

- TreatOrdinalAs:

  Passed to
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).

- min_n:

  Minimum number of complete observations required for an individual
  correlation.

- eps:

  Variance tolerance passed to
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).

- fdr_scope:

  Scope for FDR correction of correlation-comparison tests. One of
  `"matrix"`, `"per_outcome"`, or `"per_predictor"`.

- reversal_style:

  How correlations with opposite signs should be shown. One of
  `"outline"` (default), `"stripe"`, or `"none"`. Stripes are a
  static-only display option and require `ggpattern`.

- interactive:

  Optional interactive output. One of `"none"` (default), `"plotly"`,
  `"girafe"`, or `"both"`. Static ggplots are always retained in
  `Unadjusted$plot` and `FDRCorrected$plot`; requested widgets are added
  under `Interactive`.

- low_color:

  Color representing negative DeltaR values.

- mid_color:

  Color representing DeltaR = 0.

- high_color:

  Color representing positive DeltaR values.

- color_limits:

  Limits for the DeltaR color scale. The theoretical range is -2 to 2.

- triangle:

  Display "full" (default), "upper", or "lower" half of a symmetric
  comparison matrix. This affects plots only; returned matrices and
  Results remain complete.

## Value

A list containing:

- Correlations:

  The original
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
  objects for the comparison and reference groups.

- Unadjusted:

  Matrices and heatmap using raw comparison p-values.

- FDRCorrected:

  Matrices and heatmap using FDR-adjusted comparison p-values.

- Results:

  A tibble with one row per correlation pair.

- DirectionReversal:

  Logical matrix indicating opposite correlation signs between groups.

- Metadata:

  Comparison settings, group information, and the inferential
  approximation used.

- Interactive:

  Optional Plotly and/or ggiraph widgets.

## Details

Correlations are calculated using
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
so variable handling, covariate adjustment, ordinal handling, labels,
and missing-data behavior remain consistent with the SciDataReportR
correlation workflow.

Pearson correlations without covariates use the usual
independent-samples Fisher-z comparison. Spearman and residualized
partial-correlation comparisons use an approximate Fisher-z calculation;
inspect `Results$InferenceStatus` or `Metadata$InferenceStatus` before
interpreting those p-values.

When `comparison_group` and `reference_group` are both omitted for a
two-level factor, the first factor level is used as the reference group
and the second factor level as the comparison group. For non-factor
grouping variables, observed order is used. When more than two groups
are present, both groups must be specified explicitly.

## Examples

``` r
data(SampleData)

# If Sex is a factor with levels c("Male", "Female"),
# Male is automatically the reference and Female the comparison.
#
# res <- PlotCorrelationComparisons(
#   data = SampleData,
#   predictor_vars = c("age", "AXL", "Ferritin", "IL_6"),
#   outcome_vars = c("Cortisol", "Insulin"),
#   group_var = "Sex",
#   covariates = "education",
#   method = "spearman",
#   triangle = "upper"
# )
#
# res$FDRCorrected$plot
```
