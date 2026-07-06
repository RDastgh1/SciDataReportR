# Plot P-Value Comparisons

This function generates a plot comparing p-values for different
variables across two or more groups.

## Usage

``` r
PlotPValueComparisons(
  data,
  group_var,
  variables = NULL,
  VariableCategories = NULL,
  Relabel = TRUE,
  Data = lifecycle::deprecated(),
  GroupVariable = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  Data frame containing the variables to compare.

- group_var:

  Character string specifying the name of the column in `Data` that
  contains the group labels.

- variables:

  Character vector specifying the names of the columns in `Data` to
  include in the comparison. If `NULL`, all columns except
  `GroupVariable` are included.

- VariableCategories:

  Character vector specifying the categories for each variable. If
  `NULL`, no categories are used.

- Relabel:

  Logical indicating whether to replace missing labels with the column
  names.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- GroupVariable:

  **Deprecated** (since 19.15.0). Use `group_var` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A ggplot object displaying the p-value comparisons.

## Examples

``` r
data(SampleData)

PlotPValueComparisons(
  SampleData,
  group_var = "Diagnosis",
  variables = c("age", "AXL", "Adiponectin")
)
```
