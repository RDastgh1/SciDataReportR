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
  alternative = c("two.sided", "greater", "less"),
  IncludeTrendSymbol = FALSE,
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

- alternative:

  Hypothesis alternative for two-group tests. `"greater"` tests whether
  the second group factor level is greater than the first.

- IncludeTrendSymbol:

  Logical; add a plus sign for `0.05 < p <= 0.10`.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- GroupVariable:

  **Deprecated** (since 19.15.0). Use `group_var` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A ggplot object displaying the p-value comparisons.

## Details

Categorical input safeguards apply to `group_var`. Numeric columns are
converted locally to factors, with a warning, only when all observed
values are finite whole-number codes and there are at most 10 distinct
non-missing values. Other numeric inputs error: explicitly convert the
column to a factor first if categorical treatment is intended (for
example, `df_All$cohort <- factor(df_All$cohort)`). Few unique values
alone do not establish categorical meaning. Factors, character columns,
and logical columns are accepted; more than 10 observed categories warn
and continue. Entirely missing or unsupported columns error. Caller data
are unchanged, and existing restrictions on the required number of
groups still apply.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

PlotPValueComparisons(
  Labelled,
  group_var = "Diagnosis",
  variables = c(
    "age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Cortisol",
    "Ferritin", "GRO_alpha", "MMP10", "MMP7", "NT_proBNP", "PAI_1",
    "TRAIL_R3", "VEGF", "Ab_42", "p_tau", "tau"
  )
)
```
