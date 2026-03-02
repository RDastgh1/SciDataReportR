# Plot Categorical Distributions

This function creates plots to visualize the distributions of
categorical variables in a dataframe.

## Usage

``` r
PlotCategoricalDistributions(
  DataFrame,
  Variables = NULL,
  Relabel = TRUE,
  Ordinal = TRUE
)
```

## Arguments

- DataFrame:

  The dataframe containing the variables to be plotted.

- Variables:

  Optional. A character vector specifying the names of the categorical
  variables to be plotted. If NULL, categorical variables are
  automatically detected.

- Relabel:

  Logical. If TRUE, missing labels in the dataframe are replaced with
  column names as labels for plotting.

- Ordinal:

  Logical, indicating whether ordinal variables should be included.

## Value

A plot visualizing the distributions of categorical variables.
