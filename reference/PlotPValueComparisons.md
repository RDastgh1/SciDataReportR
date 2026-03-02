# Plot P-Value Comparisons

This function generates a plot comparing p-values for different
variables across two or more groups.

## Usage

``` r
PlotPValueComparisons(
  Data,
  GroupVariable,
  Variables = NULL,
  VariableCategories = NULL,
  Relabel = TRUE
)
```

## Arguments

- Data:

  Data frame containing the variables to compare.

- GroupVariable:

  Character string specifying the name of the column in `Data` that
  contains the group labels.

- Variables:

  Character vector specifying the names of the columns in `Data` to
  include in the comparison. If `NULL`, all columns except
  `GroupVariable` are included.

- VariableCategories:

  Character vector specifying the categories for each variable. If
  `NULL`, no categories are used.

- Relabel:

  Logical indicating whether to replace missing labels with the column
  names.

## Value

A ggplot object displaying the p-value comparisons.
