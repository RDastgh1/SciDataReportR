# Fit a projectable HDBSCAN model

Best for irregularly shaped numeric clusters, variable density, and data
where a meaningful noise/outlier group is expected.

## Usage

``` r
CreateClusterModel_HDBSCAN(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  minPts_range = 2:10,
  cluster_selection_epsilon_range = c(0, 0.05, 0.1),
  final_minPts = NULL,
  final_cluster_selection_epsilon = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

## Arguments

- data:

  Data frame containing numeric clustering variables.

- variables:

  Variables used for clustering.

- method:

  Either `"exploratory"` or `"finalize"`.

- minPts_range:

  Candidate minimum-points settings in exploratory mode.

- cluster_selection_epsilon_range:

  Candidate extraction epsilon values.

- final_minPts, final_cluster_selection_epsilon:

  Finalized density settings.

- ZScoreType:

  Frozen numeric preprocessing. `Scaling` is a compatibility alias.

- Scaling:

  Compatibility alias for `ZScoreType`.

- ClusterVariableName:

  Output cluster column name.

- seed:

  Random seed retained for reproducibility.

- stability_resamples:

  Number of 90% participant subsample refits used to estimate candidate
  reproducibility. Subsamples are drawn without replacement. Use `0` to
  disable stability analysis.

- stability_seed:

  Seed controlling participant subsampling.

- stability_progress:

  Whether to print subsample progress messages.

- stability_cores:

  Number of workers for stability refits. `NULL` uses all detected
  physical cores minus one, capped at the requested resamples and any
  scheduler limit. Each worker holds a refit in memory, so lower this
  setting for large SOM or high-dimensional analyses.

## Value

A fitted HDBSCAN model with cluster/noise assignments, membership
probabilities, outlier scores, frozen nearest-core support thresholds,
and subsample ARI/Jaccard and noise-recovery metrics. Persistence is
higher-is-better, noise proportion is lower-is-better, and the extracted
class count is data-derived; membership probability and outlier score
are assignment diagnostics rather than candidate fit metrics. Figures
sit beside what they describe: `fit_plot` reviews the density grid;
`ModelInfo$plots` holds `density_review`, `persistence`,
`cluster_persistence`, and `profiles`; `ModelInfo$FitDiagnostics$plots`
holds the nearest-core-distance histogram and
`outlier_score_by_cluster`, both measured against the same frozen
reference a projected case is triaged on; `ProbFit$plots` holds
membership-probability figures; and `Stability$plots` holds subsample
agreement, per-cluster recovery, and noise recovery.

## Stability output

Stability assesses internal reproducibility by full-pipeline 90%
participant subsampling without replacement. For each replicate, 90% of
complete participants are selected once, all preprocessing, any
reduction (PCA, MCA, or SOM), and the selected clustering method are
refit. Assignments from that refit are compared with the full-data
reference assignments for the same sampled participants, joined by
`.row_id`; projection is not used. It is an internal sensitivity
analysis, not independent-cohort validation.

`Stability$settings` records the analysis provenance:

- `resamples`: requested number of 90% subsample refits.

- `seed`: seed used to derive reproducible sampling and model seeds.

- `refit_scope`: always `"full_pipeline_in_sample"`, meaning
  preprocessing, reduction where applicable, and clustering were all
  refit.

- `comparison_scope`: always `"sampled_participants"`; only the
  participants used in a refit enter that replicate's agreement metrics.

- `resample_type`: `"subsample_without_replacement"` for the primary
  stability analysis.

- `resample_fraction`: the retained participant fraction, `0.90`.

- `coassignment_limit`: maximum number of complete training participants
  (2,000) for which the full pairwise co-assignment matrix is
  calculated.

- `noise_policy`: whether the method has noise labels. For ordinary
  methods it is `"all clusters included"`; HDBSCAN variants retain noise
  in ARI but exclude it from phenotype-level Jaccard recovery.

`Stability$replicates` has one row per requested refit. `Model` and
`Classes` identify the selected candidate (for HDBSCAN, `Classes` is the
data-derived extracted count); `Replicate` is its sequence number;
`Status` is `"success"` or a failure status; `SamplingSeed` and
`ModelSeed` identify the independently varied participant draw and
stochastic refit; and `Error` contains the error message for an
unsuccessful refit. `ARI` is the primary successful-refit metric: it is
invariant to cluster numbering, higher is better, and it can be negative
when agreement is worse than chance.

`Stability$cluster_recovery` has one row for each reference `Cluster` in
each successful `Replicate`. `Jaccard` is the recovery of that reference
cluster after matching it to the subsample cluster with the largest
Jaccard overlap; higher is better and one is exact recovery. `Model` and
`Classes` again identify the fitted candidate.

`Stability$summary` combines successful subsample refits:
`StabilitySuccessRate` is successful refits divided by requested refits,
an operational reliability measure that does not enter the
reproducibility score. `StabilityARI_Mean` and `StabilityARI_P05` are
the mean and fifth percentile of ARI. `StabilityJaccard_Mean` and
`StabilityJaccard_Min` are respectively the mean cluster-level recovery
and the minimum of the cluster-specific mean recoveries across
replicates. `ReproducibilityScore` is the mean of the finite
`StabilityARI_Mean` and `StabilityJaccard_Mean` values only; it does not
include success rate.

`Stability$failures` repeats the replicate columns for unsuccessful
refits, making fit failures auditable without mixing them with
successful metrics.

`Stability$participant_inclusion` is one row per complete reference
participant. `.row_id` identifies the merge-safe input row, `Cluster`
its reference assignment, `SuccessfulRefits` the number of usable
refits, and `InclusionProbability` the proportion of those refits in
which the participant returned to that cluster's label-matched subsample
cluster. `Model` and `Classes` identify the candidate. Higher inclusion
is better.

`Stability$cluster_inclusion` summarizes inclusion within each reference
`Cluster`: `MeanInclusion`, `P05Inclusion`, and `MinInclusion` are the
mean, fifth percentile, and minimum participant inclusion probabilities;
`Model` and `Classes` identify the candidate. Higher values indicate
that all, not only the average, of a cluster is recovered consistently.

`Stability$coassignment` is available only when the complete training
cohort has at most `coassignment_limit` participants. Each candidate
entry has a `status` of `"available"`, `"skipped"`, or
`"not_available"`; `reason` explains a non-available result; `matrix` is
the pairwise probability that two complete reference participants are
assigned together across successful subsample refits; and `.row_id` maps
matrix rows and columns to the same identifiers returned in
`DataWithClusters` and `ProbFit$individual`. Higher matrix values mean
more consistent pairwise co-membership. Where more than one candidate is
summarized, entries are named by its `Model_Classes` key. The matrix is
diagnostic only and is never used for selection.

`Stability$plots` contains `cluster_recovery` (per-cluster Jaccard),
`partition_metrics` (the ARI distribution), and `cluster_inclusion`; it
also contains a co-assignment heatmap when the matrix is available.
These diagnostics complement rather than replace ARI and Jaccard: none
can turn a poorly reproducible cluster into a stable phenotype.

Metric sources: Hubert and Arabie (1985) define ARI; Jaccard (1901)
defines the overlap coefficient; and Monti et al. (2003) describe
resampling-based consensus co-assignment.

## References

McInnes L, Healy J, Astels S. hdbscan: Hierarchical density based
clustering. *J Open Source Softw.* 2017;2(11):205.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
review <- CreateClusterModel_HDBSCAN(
  df_Training, c("DensityX", "DensityY"), minPts_range = c(6, 10),
  cluster_selection_epsilon_range = c(0, 0.05),
  stability_resamples = 2
)
review$ModelInfo$fit_table
#> # A tibble: 4 × 26
#>   Classes MinPts Epsilon Persistence NoiseProportion MeanMembershipProbability
#>     <int>  <dbl>   <dbl>       <dbl>           <dbl>                     <dbl>
#> 1       3      6    0           774.          0.0688                     0.710
#> 2       2     10    0          1112.          0.134                      0.654
#> 3       3      6    0.05        774.          0.0688                     0.710
#> 4       2     10    0.05       1112.          0.134                      0.654
#> # ℹ 20 more variables: MinClusterN <int>, SizeBalance <dbl>,
#> #   StabilitySuccessRate <dbl>, StabilityARI_Mean <dbl>,
#> #   StabilityARI_P05 <dbl>, StabilityJaccard_Mean <dbl>,
#> #   StabilityJaccard_Min <dbl>, NoiseSensitivity <dbl>, NoiseSpecificity <dbl>,
#> #   ReproducibilityScore <dbl>, Persistence_scaled <dbl>,
#> #   MeanMembershipProbability_scaled <dbl>, MinClusterN_scaled <dbl>,
#> #   SizeBalance_scaled <dbl>, ReproducibilityScore_scaled <dbl>, …
review$ModelInfo$AHP$recommendation
#> [1] "AHP-style review recommends HDBSCAN density setting (Classes = 2, MinPts = 10, Epsilon = 0.05). Review this advisory choice alongside the candidate plots."
review$fit_plot

review$ModelInfo$plots$density_review

model <- CreateClusterModel_HDBSCAN(
  df_Training, c("DensityX", "DensityY"), method = "finalize",
  final_minPts = 10, final_cluster_selection_epsilon = 0,
  stability_resamples = 2
)
model$ModelInfo$plots$density_review
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?

model$ModelInfo$plots$cluster_persistence

model$ModelInfo$FitDiagnostics$plots$outlier_score_by_cluster

projected <- ProjectCluster(model, df_Projection)
projected$ProjectionFit$plots$nearest_core_support

# }
```
