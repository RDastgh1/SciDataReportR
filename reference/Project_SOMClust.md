# Project new data onto an existing SOM + cluster solution

Given a fitted `Pipeline_SOMClust` object, project a new data frame
onto:

- The same Z-score scaling (via SciDataReportR::Project_ZScore()) when
  the training object used computed or projected Z-scores.

- The same SOM (kohonen), using `kohonen::map()`.

- The same node-level latent profiles and posterior probabilities.

If the training object used `ZScoreType = "PreZScored"`, this function
expects the new data already contains the same Z-score columns (by name)
stored in `object$ZScoreVars`. No re-zscoring is performed.

For each projected case, the function:

- Assigns a SOM node and distance to BMU.

- Maps node-level cluster/posterior probabilities to the individual.

- Computes an overall SOM-distance z-score and flags:

  - `Flag_SOMDist_overallHigh` – distance \> overall training p95.

  - `Flag_SOMDist_clusterHigh` – distance \> cluster-specific training
    p95 for that cluster.

Stable row id:

- `.scidr_rowid` is added to `new_df` (if not already present) and
  carried into `df_with_clusters` and `ProbFit$individual`.

- `ProbFit$individual$RowID` is set equal to `.scidr_rowid`.

Missing data:

- Only rows with complete Z-scores are mapped to the SOM.

- Rows with missing Z-scores receive NA for SOM_Node, SOM_Distance, and
  the cluster label.

## Usage

``` r
Project_SOMClust(object, new_df, ClusterName = NULL)
```

## Arguments

- object:

  A `Pipeline_SOMClust` object from a previous run.

- new_df:

  Data frame of new cases to project.

- ClusterName:

  Optional name for the cluster column; defaults to
  `object$ClusterName`. If that column already exists in `new_df` it is
  overwritten (with a message).

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
