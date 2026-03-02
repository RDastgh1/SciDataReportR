# Split violin with aligned half-boxplots and a centered significance label

Draws a split (left/right) violin for up to two groups at a single
x-position, overlays per-group boxplots aligned with each half, and
optionally annotates the plot with a p-value significance label (e.g.,
"*", "**", "***", "ns"). The y-axis title is taken from the variable's
label attribute when available (via labelled or Hmisc); otherwise the
column name is used.

## Usage

``` r
PlotSplitViolin(
  data,
  Var,
  Group,
  covars = NULL,
  nonparametric = FALSE,
  annotation_text = NULL,
  show_ns = FALSE,
  left_group = NULL,
  color_palette = NULL,
  box_offset = 0.11,
  box_width = 0.15,
  star_from = c("quantile", "data_max", "whisker"),
  star_quantile = 0.995,
  star_pad = 0.03,
  headroom = 0.08,
  star_size = 6
)
```

## Arguments

- data:

  A data frame containing `Var`, `Group`, and any covariates.

- Var:

  Column of the numeric outcome to plot. Tidy-eval friendly (unquoted or
  string).

- Group:

  Column of the grouping variable. Must have 1 or 2 unique values. No
  factor relabeling is performed inside; set levels externally if
  needed.

- covars:

  Character vector of covariate column names (default `"Age"`). Use
  `character(0)` for no covariates.

- nonparametric:

  Logical. If `FALSE` (default), compute the p-value using a linear
  model and an emmeans contrast on `Group`. If `TRUE`: with no
  covariates, use Wilcoxon rank-sum on `Var ~ Group`; with covariates,
  first residualize `Var` on the covariates (no `Group`), then apply
  Wilcoxon rank-sum to the residuals by `Group`. With only one group, no
  test is performed.

- annotation_text:

  Optional label to draw (e.g., `"*"`, `"**"`, `"***"`, `"ns"`, `"★"`).
  If `NULL`, the label is derived from the computed p-value.

- show_ns:

  Logical; if `TRUE`, show `"ns"` for non-significant results when
  `annotation_text` is `NULL`.

- left_group:

  Optional. The `Group` value to place on the **left** half. If omitted,
  the first factor level (or first observed value) is used.

- color_palette:

  Optional color specification. Supply a named vector keyed by `Group`
  values; if unnamed, values are interpreted as `(left, right)`.

- box_offset:

  Horizontal offset for the left/right boxplots from the center (single
  x = 1). Small values (e.g., `0.10-0.14`) sit the boxes inside each
  half.

- box_width:

  Boxplot width (default `0.15`). For clean separation, keep
  `box_width <= 2 * box_offset`.

- star_from:

  Where to anchor the vertical position of the significance label: one
  of `"quantile"` (default), `"data_max"`, or `"whisker"`.

- star_quantile:

  Quantile used when `star_from = "quantile"` (default `0.995`).

- star_pad:

  Fraction of the y-range added above the anchor for the label (default
  `0.03`).

- headroom:

  Fraction of the y-range added to the top limit to avoid clipping
  (default `0.08`).

- star_size:

  Text size for the significance label (default `6`).

## Value

A ggplot2 object.
