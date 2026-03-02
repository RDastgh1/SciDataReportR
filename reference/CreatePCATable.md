# Create PCA table and visualization

Perform principal component analysis (PCA) on specified variables and
create visualizations.

## Usage

``` r
CreatePCATable(
  Data,
  VarsToReduce,
  VariableCategories = NULL,
  minThresh = 0.85,
  scale = TRUE,
  center = TRUE,
  Ordinal = FALSE,
  numComponents = NULL
)
```

## Arguments

- Data:

  The dataset containing the variables for PCA.

- VarsToReduce:

  A character vector specifying the variables to include in the PCA.

- VariableCategories:

  Optional categorical vector, used to color the lollipop plot. Must be
  the same length and order as VarsToReduce if provided.

- minThresh:

  The minimum threshold for cumulative proportion of variance (default
  0.85).

- scale:

  Logical, indicating whether to scale the data (divide by SD). Default
  TRUE.

- center:

  Logical, indicating whether to center the data (subtract mean).
  Default TRUE.

- Ordinal:

  Logical, indicating whether ordinal variables should be handled.
  Currently not used inside this function.

- numComponents:

  Number of principal components to compute. If NULL, chosen by
  minThresh.

## Value

A list containing PCA results and visualizations:

- p_scree:

  Scree and variance plot (ggplot object).

- pcaresults:

  psych::principal result (after psych::fa.sort).

- LoadingTable:

  Data frame of loadings with variable names and labels.

- Scores:

  Matrix or data frame of component scores.

- CombinedData:

  Original Data with component scores appended.

- Lollipop:

  Lollipop style loading plot (ggplot object).

- ScaleParams:

  List with elements means, sds, center, scale.

- VarsUsed:

  Character vector of variables actually used in PCA.

- Center:

  Logical flag indicating if centering was used.

- Scale:

  Logical flag indicating if scaling was used.
