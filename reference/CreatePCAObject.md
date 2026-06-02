# Create a reusable PCA object and visualizations

Perform principal component analysis (PCA) on specified variables and
return reusable PCA results, scores, loading tables, combined data, and
plots.

## Usage

``` r
CreatePCAObject(
  Data,
  VarsToReduce,
  VariableCategories = NULL,
  Relabel = TRUE,
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

  A data frame containing the variables for PCA.

- VarsToReduce:

  Character vector of raw variable names to include in PCA.

- VariableCategories:

  Optional vector used to color variables in the lollipop loading plot.
  If supplied, it should be the same length and order as `VarsToReduce`.

- Relabel:

  Logical. If `TRUE`, variable labels are used in output tables and
  plots when available. If `FALSE`, raw variable names are used. Default
  is `TRUE`.

- minThresh:

  Numeric threshold for cumulative variance used to select the number of
  components when `numComponents = NULL`. Default is `0.85`.

- scale:

  Logical. If `TRUE`, variables are scaled by their standard deviation
  before PCA. Default is `TRUE`.

- center:

  Logical. If `TRUE`, variables are centered before PCA. Default is
  `TRUE`.

- Ordinal:

  Logical. Placeholder retained for backward compatibility. Currently
  not used inside this function.

- numComponents:

  Optional integer number of components to retain. If `NULL`, the number
  is selected using `minThresh`.

- Mode:

  Character. Either `"classic"` or `"omics"`. `"classic"` preserves the
  original
  [`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
  behavior. `"omics"` uses capped PCA logic for higher-dimensional data.

- backend:

  Character. PCA backend to use. Options are `"psych"`, `"prcomp"`, and
  `"irlba"`. Default is `"psych"` for compatibility.

- rotate:

  Character. Rotation method. Options are `"varimax"` and `"none"`.
  Default is `"varimax"`.

- maxComponents:

  Integer. Maximum number of final components to retain in `"omics"`
  mode when `numComponents = NULL`. Default is `20`.

- maxScreeComponents:

  Integer. Maximum number of components used to estimate and plot the
  scree curve in `"omics"` mode. Default is `20`.

- VarianceFilter:

  Optional numeric value for variance filtering before PCA in `"omics"`
  mode. If `VarianceFilterMethod = "top_n"`, keeps the top
  `VarianceFilter` most variable variables. If
  `VarianceFilterMethod = "variance_quantile"`, keeps variables with
  variance at or above the specified quantile.

- VarianceFilterMethod:

  Character. Either `"top_n"` or `"variance_quantile"`. Default is
  `"top_n"`.

- MissingnessWarningThreshold:

  Numeric threshold for warning about variable-level missingness.
  Default is `0.20`.

- ParticipantMissingnessWarningThreshold:

  Numeric threshold for warning about participant-level missingness.
  Default is `0.20`.

- imputeMethod:

  Character. Missing-data imputation method. Options are `"missRanger"`
  and `"median"`. Default is `"missRanger"`.

- SuppressWarnings:

  Logical. If `TRUE`, suppresses PCA-specific warning messages from this
  function. Default is `FALSE`.

## Value

A list with the following elements:

- p_scree:

  A ggplot scree and cumulative variance plot.

- pcaresults:

  The PCA result object. In `"classic"` mode this is a
  [`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
  result after
  [`psych::fa.sort()`](https://rdrr.io/pkg/psych/man/fa.sort.html). In
  `"omics"` mode with `"prcomp"` or `"irlba"`, this is a harmonized list
  containing loadings, scores, variance information, backend, mode, and
  rotation.

- LoadingTable:

  A data frame of component loadings with raw variable names, labels,
  and duplicate-safe plot labels.

- Scores:

  A data frame of component scores.

- CombinedData:

  The original input data with component scores appended.

- Lollipop:

  A ggplot lollipop loading plot.

- ScaleParams:

  A list containing centering and scaling parameters.

- VarsUsed:

  Character vector of variables actually used in PCA after
  preprocessing.

- VarianceTable:

  A variance table used for optional omics variance filtering, or
  `NULL`.

- Preprocessing:

  A list documenting variables dropped during preprocessing and
  missingness summaries.

- Mode:

  The PCA mode used.

- Backend:

  The PCA backend used.

- Center:

  Logical flag indicating whether centering was used.

- Scale:

  Logical flag indicating whether scaling was used.

## Details

The default `"classic"` mode preserves the original
[`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
workflow for compatibility with existing SciDataReportR analyses. The
optional `"omics"` mode is designed for high-dimensional data such as
proteomics, metabolomics, flow cytometry, FACS, transcriptomics, and
other settings where the number of variables may be large relative to
the number of participants.

The function uses raw variable names internally, but by default uses
variable labels in human-facing outputs when labels are available.
Variables with zero or undefined standard deviation after preprocessing
are dropped automatically with a warning so PCA can continue without
manual preprocessing.

## Examples

``` r
PCA <- CreatePCAObject(
  Data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3,
  imputeMethod = "median"
)

PCA_raw_names <- CreatePCAObject(
  Data = mtcars,
  VarsToReduce = names(mtcars),
  Relabel = FALSE,
  numComponents = 3,
  imputeMethod = "median"
)

PCA_omics <- CreatePCAObject(
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
