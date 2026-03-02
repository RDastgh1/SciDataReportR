# Plot Associations

This function generates scatter plots or box plots to visualize the
relationship between two variables.

## Usage

``` r
PlotAssociations(DataFrame, Var1, Var2, Ordinal = FALSE)
```

## Arguments

- DataFrame:

  The data frame containing the variables of interest.

- Var1:

  The name of the first variable.

- Var2:

  The name of the second variable.

- Ordinal:

  Logical, indicating whether ordinal variables should be included.

## Value

A ggplot object representing the relationship between the variables.
