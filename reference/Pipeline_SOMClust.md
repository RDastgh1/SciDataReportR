# SOM + latent profile clustering pipeline (with AHP and distance baselines)

End-to-end pipeline to:

- Standardize variables using SciDataReportR::CalcZScore() or a supplied
  Z-score object.

- Fit a Self-Organizing Map (SOM; kohonen) on complete cases.

- Generate aweSOM visualizations (Circular, Line, Cloud) with optional
  relabeling using variable labels from the original data frame.

- Cluster SOM codebook vectors using latent profile analysis (tidyLPA /
  mclust backend).

- In `method = "exploratory"`, fit a grid of models and select a
  recommended solution using an Analytic Hierarchy Process (AHP)-style
  index combining AIC, BIC, and Entropy.

- In `method = "finalize"`, fit a user-specified model and number of
  profiles.

- Map node-level clusters and posterior probabilities back to
  individuals.

Missing data:

- SOM and clustering are fit only on rows with complete Z-scores.

- The returned `df_with_clusters` has the full original data with
  `.scidr_rowid` and a single cluster column appended; rows not used in
  SOM/LPA get NA.

Stable row id:

- `.scidr_rowid` is added to the input data and carried into
  `df_with_clusters` and `ProbFit$individual`.

- `ProbFit$individual$RowID` is set equal to `.scidr_rowid` so merges do
  not rely on row order.

Z-score behavior:

- `ZScoreType = "Center and Scale"/"Center Only"/"Scale Only"` computes
  Z-scores from `df` via
  [`CalcZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/CalcZScore.md).

- `ZScoreType = "ZScoreObj"` projects Z-scores using an external
  `ZScoreObj` via
  [`Project_ZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/Project_ZScore.md).

- `ZScoreType = "PreZScored"` uses existing Z-score columns in `df`
  as-is and does not re-zscore.

## Usage

``` r
Pipeline_SOMClust(
  df,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:15,
  models = c(1, 2, 3),
  final_k = NULL,
  final_model = NULL,
  ClusterName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  seed_som = 934521L,
  seed_lpa = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_col = NULL
)
```

## Arguments

- df:

  Data frame containing the variables to be used in SOM and clustering.

- variables:

  Optional character vector of variable names. If NULL, numeric
  variables are auto-detected using
  `SciDataReportR::getNumVars(df, Ordinal = FALSE)`. In
  `ZScoreType = "PreZScored"`, this can also be NULL if you supply
  `ZScoreVars` or if Z-score columns can be auto-detected by prefix.

- method:

  One of `"exploratory"` (default) or `"finalize"`. In `"exploratory"`,
  a grid of models is fit and AHP chooses the recommended solution. In
  `"finalize"`, the user must specify `final_k` and `final_model`.

- k_range:

  Integer vector of numbers of clusters/profiles to consider in
  exploratory mode. Default `2:15`.

- models:

  Integer vector of model specifications for tidyLPA (mclust backend).
  Default `c(1, 2, 3)`.

- final_k:

  Integer; number of profiles for `method = "finalize"`.

- final_model:

  Integer; model specification for `method = "finalize"` (should be one
  of `models`).

- ClusterName:

  Name of the cluster column in the output. Defaults to `"Cluster"`. If
  this column already exists in `df`, it is overwritten (with a
  message).

- ZScoreType:

  One of:

  - `"Center and Scale"` (default)

  - `"Center Only"`

  - `"Scale Only"`

  - `"ZScoreObj"` (use an existing ZScore object)

  - `"PreZScored"` (use existing Z-score columns in df as-is)

- ZScoreObject:

  Optional ZScoreObj (from
  [`CalcZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/CalcZScore.md)
  or
  [`Project_ZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/Project_ZScore.md))
  to use when `ZScoreType = "ZScoreObj"`.

- som_xdim, som_ydim:

  Optional integers for SOM grid dimensions. If NULL, a square grid with
  side length `ceiling(n_complete^(1/3))` is used.

- som_topo:

  SOM topology for `kohonen::somgrid()`, default `"hexagonal"`.

- som_neigh:

  SOM neighbourhood function, default `"gaussian"`.

- seed_som, seed_lpa:

  Integer seeds for SOM and LPA steps (defaults 934521 and 93421).

- Relabel:

  Logical; if TRUE (default), aweSOM plots are relabeled using variable
  labels from the *original* `df` (via Hmisc or sjlabelled when
  available) by stripping the Z-score prefix.

- ZScorePrefix:

  Character prefix used for Z-score columns when
  `ZScoreType = "PreZScored"`. Default `"Z_"`.

- ZScoreVars:

  Optional character vector of Z-score column names to use when
  `ZScoreType = "PreZScored"`. If NULL, the function attempts to infer
  them from `variables` or by detecting columns starting with
  `ZScorePrefix`.

- id_col:

  Optional character scalar. If provided and present in `df`, this
  column is carried into `ProbFit$individual` for convenience.

## Value

A list of class `"Pipeline_SOMClust"` with components:

- `method`, `vars_used`, `ZScoreType`, `ZScoreObject`, `ZScoreVars`,
  `ClusterName`

- `complete_rows`: logical vector (rows used for SOM/LPA)

- `df_with_clusters`: original `df` with `.scidr_rowid` and only the
  cluster column appended

- `fit_plot`: ggplot of AIC/BIC/Entropy vs k and model

- `ModelInfo_SOM`: list with `som_model`, `som_codes`, `som_grid`,
  `SOMFit` (distance diagnostics, baselines, and per-cluster flags),
  `plots` (aweSOM plots)

- `ModelInfo_MClust`: list with `lpa_models`, `fit_table`, and `AHP`
  information

- `ProbFit`: list with `node` (node-level posterior probabilities),
  `individual` (per-person mapping and probabilities, including
  `.scidr_rowid`), and probability plots

## Details

The AHP-style index is computed by:

1.  Scaling AIC, BIC, and Entropy across candidate solutions (AIC/BIC
    are negated so that lower values correspond to better fit; higher
    scaled scores are preferred).

2.  Taking the mean of the three scaled indices. The model with the
    highest AHP index is recommended.

## References

Saaty TL. *The Analytic Hierarchy Process*. McGraw–Hill, 1980.
