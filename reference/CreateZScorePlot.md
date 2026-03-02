# Create a Z-Score Plot with Statistical Significance

This function generates a Z-score plot to compare multiple variables
across different groups. It offers options for parametric or
non-parametric tests, ordinal variable conversion, and custom labeling.
Significant p-values and FDR-adjusted p-values are highlighted on the
plot.

## Usage

``` r
CreateZScorePlot(
  Data,
  TargetVar,
  Variables,
  VariableCategories = NULL,
  Relabel = TRUE,
  sort = TRUE,
  RemoveXAxisLabels = TRUE,
  Ordinal = TRUE,
  Parametric = TRUE,
  SigP_YCoord = 1.5,
  SigFDR_YCoord = 1.6
)
```

## Arguments

- Data:

  A dataframe containing the data to be analyzed.

- TargetVar:

  A string specifying the column name of the grouping variable.

- Variables:

  A vector of strings specifying the column names of the variables to be
  analyzed.

- VariableCategories:

  An optional vector categorizing the variables.

- Relabel:

  Logical; if TRUE, variables will be relabeled using their labels from
  the dataframe.

- sort:

  Logical; if TRUE, variables will be sorted by category and p-value.

- RemoveXAxisLabels:

  Logical; if TRUE, X-axis labels will be removed.

- Ordinal:

  Logical; if TRUE, ordinal variables will be converted to numeric.

- Parametric:

  Logical; if TRUE, parametric tests (t-test/ANOVA) will be used;
  otherwise, non-parametric tests (Wilcoxon/Kruskal-Wallis) will be
  used.

- SigP_YCoord:

  Numeric; the y-coordinate for marking significant p-values.

- SigFDR_YCoord:

  Numeric; the y-coordinate for marking significant FDR-adjusted
  p-values.

## Value

A ggplot object representing the Z-score plot.
