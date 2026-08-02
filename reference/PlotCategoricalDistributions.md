# Plot categorical distributions

This function creates plots to visualize the distributions of
categorical variables in a dataframe.

## Usage

``` r
PlotCategoricalDistributions(
  data,
  variables = NULL,
  Relabel = TRUE,
  Ordinal = lifecycle::deprecated(),
  TreatOrdinalAs = "Categorical",
  LabelType = "percent",
  MissingLabel = "Missing",
  DataFrame = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataframe containing the variables to be plotted.

- variables:

  Optional. A character vector specifying the names of the categorical
  variables to be plotted. If NULL, categorical variables are
  automatically detected.

- Relabel:

  Logical. If TRUE, missing labels in the dataframe are replaced with
  column names as labels for plotting.

- Ordinal:

  Deprecated logical compatibility option; use `TreatOrdinalAs` instead.

- TreatOrdinalAs:

  How ordinal variables are handled. This categorical plot accepts
  `"Categorical"` or `"Exclude"`.

- LabelType:

  Character. Either "percent" or "count", indicating what should be
  shown on the x-axis and inside the bars.

- MissingLabel:

  Character label to use for missing values.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A ggplot object visualizing the distributions of categorical variables.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

PlotCategoricalDistributions(
  Labelled,
  variables = c("Diagnosis", "Genotype")
)
```
