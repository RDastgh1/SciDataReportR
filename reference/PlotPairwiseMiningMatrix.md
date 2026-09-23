# Plot non-directional mixed-type pairwise mining matrices

Builds referent-centred phenotype screening matrices that show the
magnitude of pairwise separation only. Signed continuous contrasts and
category-level prevalence contrasts are retained in returned audit
objects for follow-up, but are deliberately not displayed in the mining
matrices.

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
  adjusted_outline = TRUE,
  adjusted_significance_threshold = 0.05,
  adjusted_outline_color = "black",
  adjusted_outline_linewidth = 1,
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

  Optional character vector of covariates. Continuous contrasts use
  these covariates; Cramer's V screening remains unadjusted.

- adjust_scope:

  Multiple-comparison correction scope: `"per_group"`, `"per_variable"`,
  `"matrix"`, or `"none"`.

- p_adjust_method:

  Method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- star_p:

  Which p-values drive cell stars: `"raw"`, `"adjusted"`, or `"none"`.

- adjusted_outline:

  Logical; outline cells significant after adjustment.

- adjusted_significance_threshold:

  Threshold for adjusted-significant outlines.

- adjusted_outline_color, adjusted_outline_linewidth:

  Appearance of the adjusted-significant outline.

- variable_metadata:

  Optional data frame with `Variable` and optional `AnchorN`,
  `AlignedN`, `FilledPrior`, `FilledFuture`, and `TimeExtended`.

- max_levels:

  Maximum categorical levels allowed for one variable.

- continuous_fill_limits, categorical_fill_limits:

  Optional non-negative plotting limits for the continuous and
  categorical magnitude matrices.

- x_axis_text_angle:

  Numeric angle for phenotype labels.

- row_label_width:

  Approximate character width used to wrap row labels.

## Value

An object of class `"SciDataReportRPairwiseMiningMatrix"` with `Plots`,
`Results`, `ContinuousAudit`, `CategoryLevelResults`, `Settings`,
`Models`, and `Warnings`. `Plots$Continuous` displays absolute
referent-SD mean differences and `Plots$Categorical` displays pairwise
Cramer's V. Both use an independent pale-gray-to-navy magnitude scale.

## Details

The required Referent determines the phenotype comparison columns.
Continuous cells are absolute `Group - Referent` contrasts in
reference-SD units. Categorical cells are pairwise Cramer's V values
from contingency tables, one row per variable. These metrics
intentionally have separate legends and should not be compared as
interchangeable effect sizes.
