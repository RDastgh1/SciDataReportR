# Plot categorical distributions

This function creates plots to visualize the distributions of
categorical variables in a dataframe.

## Usage

``` r
PlotCategoricalDistributions(
  DataFrame,
  Variables = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  LabelType = "percent",
  MissingLabel = "Missing"
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

- LabelType:

  Character. Either "percent" or "count", indicating what should be
  shown on the x-axis and inside the bars.

- MissingLabel:

  Character label to use for missing values.

## Value

A ggplot object visualizing the distributions of categorical variables.
