# Plot a mixed-type pairwise mining matrix

Builds referent-centred phenotype contrasts for continuous and
categorical measures. It follows
[`MakePairwiseHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/MakePairwiseHeatmap.md)
conventions while retaining categorical category-level contrasts instead
of silently discarding them.

## Usage

``` r
PlotPairwiseMiningMatrix(
  data,
  group_var,
  variables,
  Referent,
  covariates = NULL,
  adjust_scope = c("per_group", "per_variable", "matrix", "none"),
  p_adjust_method = c("fdr", "bonferroni", "holm", "none"),
  star_p = c("raw", "adjusted", "none"),
  variable_metadata = NULL,
  max_levels = 30,
  continuous_fill_limits = NULL,
  categorical_fill_limits = NULL,
  x_axis_text_angle = 0,
  row_label_width = 58
)
```

## Arguments

- data:

  A data frame with labelled variables.

- group_var:

  Character scalar naming the grouping variable.

- variables:

  Character vector of continuous or categorical variables.

- Referent:

  Character scalar naming the referent level of `group_var`.

- covariates:

  Optional character vector of covariates.

- adjust_scope:

  Multiple-comparison correction scope: `"per_group"`, `"per_variable"`,
  `"matrix"`, or `"none"`.

- p_adjust_method:

  Method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- star_p:

  Which p-values drive cell stars: `"raw"`, `"adjusted"`, or `"none"`.

- variable_metadata:

  Optional data frame with `Variable` and optional `AnchorN`,
  `AlignedN`, `FilledPrior`, `FilledFuture`, and `TimeExtended` columns.
  These fields are joined to the returned audit table.

- max_levels:

  Maximum categorical levels allowed for one variable.

- continuous_fill_limits, categorical_fill_limits:

  Optional symmetric plotting limits for the continuous and categorical
  matrices.

- x_axis_text_angle:

  Numeric angle for phenotype labels.

- row_label_width:

  Approximate character width used to wrap displayed variable and
  category labels.

## Value

An object of class `"SciDataReportRPairwiseMiningMatrix"` with `Plots`,
`Results`, `OmnibusResults`, `Settings`, `Models`, and `Warnings`.
`Plots$Continuous` is a reference-SD mean-difference matrix;
`Plots$Categorical` is a percentage-point prevalence-difference matrix;
and `Plots$Omnibus` summarizes the overall phenotype association for
every source variable.

## Details

Continuous cells are `Group - Referent` contrasts after scaling each
outcome to the referent mean and standard deviation. Categorical
variables expand to one row per observed level, with cells equal to
`Group - Referent` prevalence in percentage points. The two plot types
intentionally have separate scales. Do not compare their colour
magnitudes as if they were the same effect size.
