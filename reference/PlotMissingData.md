# Plot Missing Data

A custom function to plot missing data with variables as rows and
observations as columns.

## Usage

``` r
PlotMissingData(
  DataFrame,
  Variables = NULL,
  Relabel = TRUE,
  show_perc = TRUE,
  show_perc_var = TRUE,
  cluster = FALSE
)
```

## Arguments

- DataFrame:

  The dataframe containing missing data.

- Variables:

  A character vector specifying the variables to be included in the
  plot. If NULL (default), all variables in the dataframe will be used.

- Relabel:

  Logical indicating whether to relabel the columns based on their
  labels. Default is TRUE.

- show_perc:

  Logical indicating whether to show percentage labels in legend.
  Default is TRUE.

- show_perc_var:

  Logical indicating whether to show percentage of missing for each
  variable. Default is TRUE.

- cluster:

  Logical indicating whether to cluster rows by missingness pattern.
  Default is FALSE.

## Value

A ggplot object visualizing missing data with variables as rows.
