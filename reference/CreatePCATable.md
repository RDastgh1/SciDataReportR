# Create PCA table and visualization

Perform principal component analysis (PCA) on specified variables and
create visualizations. The default `"classic"` mode preserves the
original
[`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
workflow for compatibility with existing SciDataReportR analyses. The
optional `"omics"` mode is designed for high-dimensional data such as
proteomics, metabolomics, and other settings where the number of
variables is much larger than the number of participants.

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
  numComponents = NULL,
  Mode = c("classic", "omics"),
  backend = c("psych", "prcomp", "irlba"),
  rotate = c("varimax", "none"),
  maxComponents = 20,
  maxScreeComponents = 20,
  VarianceFilter = NULL,
  VarianceFilterMethod = c("top_n", "variance_quantile"),
  MissingnessWarningThreshold = 0.2,
  ParticipantMissingnessWarningThreshold = 0.2,
  imputeMethod = c("missRanger", "median"),
  SuppressWarnings = FALSE
)
```

## Arguments

- Data:

  The dataset containing the variables for PCA.

- VarsToReduce:

  A character vector specifying the variables to include in the PCA.

- VariableCategories:

  Optional categorical vector, used to color the lollipop plot. Must be
  the same length and order as `VarsToReduce` if provided.

- minThresh:

  The minimum threshold for cumulative proportion of variance. Default
  is `0.85`. In `"omics"` mode, this threshold is evaluated only across
  the capped scree components.

- scale:

  Logical, indicating whether to scale the data by standard deviation.
  Default is `TRUE`.

- center:

  Logical, indicating whether to center the data by subtracting the
  mean. Default is `TRUE`.

- Ordinal:

  Logical, indicating whether ordinal variables should be handled.
  Currently not used inside this function.

- numComponents:

  Number of principal components to compute. If `NULL`, chosen by
  `minThresh` in `"classic"` mode and by capped scree results in
  `"omics"` mode.

- Mode:

  Character. Either `"classic"` or `"omics"`. `"classic"` preserves the
  original
  [`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
  behavior. `"omics"` uses capped PCA logic for high-dimensional data.

- backend:

  Character. PCA backend to use. Options are `"psych"`, `"prcomp"`, and
  `"irlba"`. Default is `"psych"` to preserve backward compatibility. In
  `"omics"` mode, `"irlba"` is recommended for speed.

- rotate:

  Character. Rotation method. Options are `"varimax"` and `"none"`.
  Default is `"varimax"` to preserve interpretability and existing
  loading-table behavior.

- maxComponents:

  Maximum number of final components to retain in `"omics"` mode when
  `numComponents` is `NULL`. Default is `20`.

- maxScreeComponents:

  Maximum number of components used to estimate and plot the scree curve
  in `"omics"` mode. Default is `20`.

- VarianceFilter:

  Optional numeric value for variance filtering before PCA in `"omics"`
  mode. If `VarianceFilterMethod = "top_n"`, this keeps the top
  `VarianceFilter` most variable variables. If
  `VarianceFilterMethod = "variance_quantile"`, this keeps variables
  with variance at or above the specified quantile.

- VarianceFilterMethod:

  Character. Either `"top_n"` or `"variance_quantile"`.

- MissingnessWarningThreshold:

  Numeric threshold for warning about high variable-level missingness.
  Default is `0.20`.

- ParticipantMissingnessWarningThreshold:

  Numeric threshold for warning about high participant-level
  missingness. Default is `0.20`.

- imputeMethod:

  Character. Missing-data imputation method. Options are `"missRanger"`
  and `"median"`. Default is `"missRanger"` to preserve the original
  behavior.

- SuppressWarnings:

  Logical. If `TRUE`, suppresses PCA-specific warning messages from this
  function. Default is `FALSE`.

## Value

A list containing PCA results and visualizations:

- p_scree:

  Scree and variance plot as a ggplot object.

- pcaresults:

  PCA result object. In `"classic"` mode this is the
  [`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
  result after
  [`psych::fa.sort()`](https://rdrr.io/pkg/psych/man/fa.sort.html). In
  `"omics"` mode this is a harmonized list with rotated and sorted
  loadings, scores, and variance information.

- LoadingTable:

  Data frame of rotated loadings with variable names and labels.

- Scores:

  Matrix or data frame of component scores.

- CombinedData:

  Original data with component scores appended.

- Lollipop:

  Lollipop-style loading plot as a ggplot object.

- ScaleParams:

  List with means, standard deviations, center, and scale.

- VarsUsed:

  Character vector of variables actually used in PCA.

- VarianceTable:

  Table of variance explained by component.

- Mode:

  PCA mode used.

- Backend:

  PCA backend used.

- Center:

  Logical flag indicating if centering was used.

- Scale:

  Logical flag indicating if scaling was used.

## Details

In `"omics"` mode, the function avoids computing one component per input
variable. Instead, it computes a capped number of components using a
faster PCA backend, applies optional variance filtering, and then
applies rotation and sorting so loadings remain interpretable. This is
useful when PCA is being used for exploratory structure, visualization,
latent biological signal, or clustering rather than exhaustive factor
decomposition.

## Examples

``` r
PCA <- CreatePCATable(
  Data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

PCA_Omics <- CreatePCATable(
  Data = mtcars,
  VarsToReduce = names(mtcars),
  Mode = "omics",
  backend = "prcomp",
  maxComponents = 5,
  maxScreeComponents = 5,
  imputeMethod = "median"
)
#> Warning: STATS is longer than the extent of 'dim(x)[MARGIN]'
```
