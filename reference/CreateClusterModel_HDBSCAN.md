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
  stability_progress = FALSE
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

  Number of 80% participant subsample refits used to estimate candidate
  reproducibility. Subsamples are drawn without replacement. Use `0` to
  disable stability analysis.

- stability_seed:

  Seed controlling participant subsampling.

- stability_progress:

  Whether to print subsample progress messages.

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

Stability assesses internal reproducibility by full-pipeline 80%
participant subsampling without replacement. For each replicate, 80% of
complete participants are selected once, all preprocessing and any
reduction (PCA, MCA, or SOM) are refit, the selected clustering method
is refit, and the original complete training participants are projected
into that subsample fit for comparison with the original fitted
partition. It is an internal sensitivity analysis, not
independent-cohort validation.

`Stability$settings` records the analysis provenance:

- `resamples`: requested number of 80% subsample refits.

- `seed`: seed used to select subsamples.

- `refit_scope`: always `"full_pipeline"`, meaning preprocessing,
  reduction where applicable, and clustering were all refit.

- `resample_type`: `"subsample_without_replacement"` for the primary
  stability analysis.

- `resample_fraction`: the retained participant fraction, `0.80`.

- `coassignment_limit`: maximum number of complete training participants
  (2,000) for which the full pairwise co-assignment matrix is
  calculated.

- `noise_policy`: whether the method has noise labels. For ordinary
  methods it is `"all clusters included"`; HDBSCAN variants retain noise
  in global partition metrics but exclude it from per-cluster inclusion
  and co-assignment summaries.

`Stability$replicates` has one row per requested refit. `Model` and
`Classes` identify the selected candidate (for HDBSCAN, `Classes` is the
data-derived extracted count); `Replicate` is its sequence number;
`Status` is `"success"` or a failure status; and `Error` contains the
error message for an unsuccessful refit. Successful rows contain these
partition metrics:

- `ARI`: adjusted Rand index, agreement corrected for chance; higher is
  better and can be negative when agreement is worse than chance.

- `VI`: variation of information, the information lost or gained when
  changing partitions; lower is better and zero is identical.

- `NMI`: normalized mutual information; higher is better and one is
  identical.

- `FowlkesMallows`: pairwise clustering agreement; higher is better and
  one is identical.

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
`StabilityJaccard_Min` are respectively the mean and minimum
label-matched Jaccard recovery. `ReproducibilityScore` is the mean of
the finite `StabilityARI_Mean` and `StabilityJaccard_Mean` values only;
it does not include success rate, VI, NMI, or Fowlkes–Mallows.

`Stability$failures` repeats the replicate columns for unsuccessful
refits, making fit failures auditable without mixing them with
successful metrics.

`Stability$participant_inclusion` is one row per complete reference
participant. `RowIndex` identifies its original row position, `Cluster`
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
assigned together across successful subsample refits; and `row_ids` maps
matrix rows and columns to original training-row positions. Higher
matrix values mean more consistent pairwise co-membership. Where more
than one candidate is summarized, entries are named by its
`Model_Classes` key. The matrix is diagnostic only and is never used for
selection.

`Stability$plots` contains `cluster_recovery` (per-cluster Jaccard),
`partition_metrics` (ARI, VI, NMI, and Fowlkes–Mallows distributions),
and `cluster_inclusion`; it also contains a co-assignment heatmap when
the matrix is available. These diagnostics complement rather than
replace ARI and Jaccard: none can turn a poorly reproducible cluster
into a stable phenotype.

Metric sources: Hubert and Arabie (1985) define ARI; Jaccard (1901)
defines the overlap coefficient; Meila (2005) defines VI; Strehl and
Ghosh (2002) describe NMI for partition comparison; Fowlkes and Mallows
(1983) define their pairwise index; and Monti et al. (2003) describe
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
