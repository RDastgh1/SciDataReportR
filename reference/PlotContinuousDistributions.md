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
  ncol = 3,
  Ordinal = TRUE
)
```

## Arguments

- DataFrame:

  A data frame containing the variables to be plotted.

- Variables:

  Character vector of column names to plot. If NULL, all numeric
  variables (or numeric + ordinal if Ordinal = TRUE) are used.

- Fill:

  *Optional.* Single column name to map to the `fill`/`colour` aesthetic
  (e.g., a factor such as sex, treatment group, etc.).

- Relabel:

  Logical; if TRUE, facet labels use variable labels when available
  (sjlabelled).

- ncol:

  Number of columns in the facet grid.

- Ordinal:

  Logical; include labelled-ordinal variables as numeric?

## Value

         A `ggplot` object.
