# Project new data onto an existing SOM clinical phenotype space

Train once, project many: a reusable unsupervised clinical phenotyping
framework that learns phenotype structure in a training cohort and
projects new participants into the fixed phenotype space while
quantifying membership uncertainty and projection fit. Given a fitted
[`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md)
object, project a new data frame onto:

- The same Z-score scaling (via SciDataReportR::ProjectZScore()) when
  the training object used computed or projected Z-scores.

- The same SOM (kohonen), using
  [`kohonen::map()`](https://rdrr.io/pkg/kohonen/man/map.kohonen.html).

- The same node-level latent profiles and posterior probabilities.

If the training object used `ZScoreType = "PreZScored"`, this function
expects the new data already contains the same Z-score columns (by name)
stored in `object$ZScoreVars`. No re-zscoring is performed.

For each projected case, the function:

- Assigns a SOM node and distance to BMU.

- Maps node-level cluster/posterior probabilities to the individual.

- Computes an overall SOM-distance z-score and flags:

  - `Flag_SOMDist_overallHigh` - distance \> overall training p95.

  - `Flag_SOMDist_clusterHigh` - distance \> cluster-specific training
    p95 for that cluster.

Stable row id:

- `.scidr_rowid` is added to `new_df` (if not already present) and
  carried into `df_with_clusters` and `ProbFit$individual`.

- `ProbFit$individual$RowID` is set equal to `.scidr_rowid`.

Missing data:

- Only rows with complete Z-scores are mapped to the SOM.

- Rows with missing Z-scores receive NA for SOM_Node, SOM_Distance, and
  the cluster label.

`Project_SOMClust()` has been superseded by `ProjectSOMCluster()`. It
remains available as a backwards-compatible alias and returns the same
projection object.

## Usage

``` r
ProjectSOMCluster(
  object,
  new_df,
  ClusterName = NULL,
  high_dist_quantile = 0.95,
  low_prob_threshold = 0.7
)

Project_SOMClust(...)
```

## Arguments

- object:

  A SOM cluster model object from
  [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md).

- new_df:

  Data frame of new cases to project.

- ClusterName:

  Optional name for the cluster column; defaults to
  `object$ClusterName`. If that column already exists in `new_df` it is
  overwritten (with a message).

- high_dist_quantile:

  Numeric value between 0 and 1 used to define high SOM-distance flags
  from the training distance distribution. Default is `0.95`.

- low_prob_threshold:

  Numeric posterior probability threshold used to flag uncertain
  phenotype membership. Default is `0.70`.

## Value

A list of class `"Project_SOMClust"` with components:

- `vars_used`, `ClusterName`, `complete_rows`

- `df_with_clusters`: `new_df` with `.scidr_rowid` and only the cluster
  column appended.

- `SOMProj`: list with training and projected distance summaries,
  cluster-level flag summaries, comparison, and plots.

- `ProbFit`: list with `node` (training node-level posterior info),
  `individual` (projection-level info including distance flags and
  z-scores), and probability plots.

- `ModelInfo_SOM`, `ModelInfo_MClust`: references to the original model
  objects for convenience.

## Examples

``` r
if (FALSE) { # \dontrun{
# NOTE: Not run - projection requires a trained SOM model, and
# CreateSOMClusterModel() currently errors on the tracked get_data() bug
# (see CreateSOMClusterModel() for details).
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

model <- CreateSOMClusterModel(
  data = Labelled,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  method = "finalize",
  final_k = 3,
  final_model = 1
)

projected <- ProjectSOMCluster(object = model, new_df = Labelled)
} # }
```
