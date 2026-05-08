# Plot Continuous Distributions

Creates rain-cloud plots (half-violin + box/median + scatter) for one or
more continuous variables, with optional group-wise colouring.

## Usage

``` r
PlotContinuousDistributions(
  DataFrame,
  Variables = NULL,
  Fill = NULL,
  Relabel = TRUE,
  FacetLabelStyle = c("both", "label_only", "variable_only", "auto"),
  ncol = 3,
  Ordinal = TRUE
)
```

## Arguments

- DataFrame:

  A data frame containing the variables to be plotted.

- Variables:

  Character vector of column names to plot.

- Fill:

  Optional column name for grouping.

- Relabel:

  Logical; use variable labels when available.

- FacetLabelStyle:

  One of "both", "label_only", "variable_only", "auto".

- ncol:

  Number of columns in the facet grid.

- Ordinal:

  Logical; include labelled-ordinal variables as numeric.

## Value

A ggplot object.
