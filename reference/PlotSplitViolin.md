# Split violin with aligned half-boxplots, significance label, sample sizes, and label-aware title

Draws a split (left/right) violin for up to two groups at a single
x-position, overlays per-group boxplots aligned with each half, and
optionally annotates the plot with a p-value significance label.
Supports displaying sample sizes (n) and automatically using variable
labels for axis and title when available.

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
  plot_title = NULL,
  use_var_label_as_title = FALSE,
  show_n = TRUE,
  n_position = c("legend", "top"),
  n_size = 3.5,
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

  A data frame containing `Var`, `Group`, and optional covariates.

- Var:

  Numeric outcome variable (tidy-eval).

- Group:

  Grouping variable (\<= 2 unique values).

- covars:

  Character vector of covariates (default `NULL`).

- nonparametric:

  Logical. If `FALSE`, uses linear model + emmeans contrast. If `TRUE`,
  uses Wilcoxon (with residualization if covariates are present).

- annotation_text:

  Optional manual annotation (e.g., "\*", "ns").

- show_ns:

  Logical; if `TRUE`, display "ns" for non-significant results.

- plot_title:

  Optional custom plot title.

- use_var_label_as_title:

  Logical; if `TRUE`, uses variable label as title.

- show_n:

  Logical; if `TRUE`, display sample size per group.

- n_position:

  Where to display n: "legend" or "top".

- n_size:

  Text size for n labels when shown on top.

- left_group:

  Optional group to force on left side.

- color_palette:

  Optional named vector of colors.

- box_offset:

  Horizontal offset for boxplots.

- box_width:

  Width of boxplots.

- star_from:

  Position method ("quantile","data_max","whisker").

- star_quantile:

  Quantile used for placement.

- star_pad:

  Padding above anchor.

- headroom:

  Extra space above data.

- star_size:

  Text size for annotation.

## Value

A ggplot2 object.
