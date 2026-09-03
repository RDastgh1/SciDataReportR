# Fit a projectable K-means clustering model

Best for approximately spherical, similarly sized groups in a numeric
feature space; use scaling unless all variables are commensurate.

## Usage

``` r
CreateClusterModel_KMeans(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  nstart = 50L,
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

- k_range:

  Candidate cluster counts in exploratory mode.

- final_k:

  Number of clusters for a finalized K-means solution.

- ZScoreType:

  Frozen numeric preprocessing. `Scaling` is a compatibility alias.

- Scaling:

  Compatibility alias for `ZScoreType`.

- ClusterVariableName:

  Output cluster column name.

- seed:

  Random seed retained for reproducibility.

- nstart:

  Number of random K-means starts.

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

A fitted model with `ModelInfo$fit_table` (WSS, between-cluster sum of
squares, silhouette, Calinski-Harabasz index, and minimum cluster size)
and subsample stability. WSS is lower-is-better; between-cluster sum of
squares, silhouette, and Calinski-Harabasz are higher-is-better. Figures
sit beside what they describe: `fit_plot` reviews candidates;
`ModelInfo$plots` holds `elbow`, `silhouette` (the per-participant
silhouette profile of the selected solution), `silhouette_by_k`,
`calinski_harabasz`, `centre_heatmap`, `centre_profile`, and `profiles`;
`ModelInfo$FitDiagnostics$plots` holds the distance-to-centroid
histogram; `ProbFit$plots` holds assignment-margin figures; and
`Stability$plots` holds cluster-recovery and complementary stability
diagnostics. `ModelInfo$ReviewSpace` freezes the two-dimensional review
space shared by the training and projection maps; it is diagnostic only
and does not affect clustering.

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

MacQueen JB. Some methods for classification and analysis of
multivariate observations. In: *Proceedings of the Fifth Berkeley
Symposium*; 1967.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
vars_Numeric <- paste0("Var", 1:12)
review <- CreateClusterModel_KMeans(
  df_Training, vars_Numeric, k_range = 2:5, nstart = 10,
  stability_resamples = 2
)
review$ModelInfo$fit_table
#> # A tibble: 4 × 23
#>   Classes   WSS BetweenSS CalinskiHarabasz Silhouette MinClusterN SizeBalance
#>     <int> <dbl>     <dbl>            <dbl>      <dbl>       <int>       <dbl>
#> 1       2 2914.      914.             99.8      0.272          80       0.333
#> 2       3 2073.     1755.            134.       0.373          80       0.5  
#> 3       4 1264.     2564.            214.       0.490          80       1    
#> 4       5 1118.     2710.            191.       0.455          23       0.288
#> # ℹ 16 more variables: StabilitySuccessRate <dbl>, StabilityARI_Mean <dbl>,
#> #   StabilityARI_P05 <dbl>, StabilityJaccard_Mean <dbl>,
#> #   StabilityJaccard_Min <dbl>, NoiseSensitivity <dbl>, NoiseSpecificity <dbl>,
#> #   ReproducibilityScore <dbl>, Silhouette_scaled <dbl>,
#> #   CalinskiHarabasz_scaled <dbl>, MinClusterN_scaled <dbl>,
#> #   SizeBalance_scaled <dbl>, ReproducibilityScore_scaled <dbl>,
#> #   WSS_scaled <dbl>, ahp_index <dbl>, Recommended <lgl>
review$ModelInfo$AHP$recommendation
#> [1] "AHP-style review recommends K-means k (Classes = 4). Review this advisory choice alongside the candidate plots."
review$fit_plot

review$ModelInfo$plots$elbow

review$ModelInfo$plots$silhouette_by_k

model <- CreateClusterModel_KMeans(
  df_Training, vars_Numeric, method = "finalize", final_k = 4,
  nstart = 10, stability_resamples = 2
)
model$ModelInfo$plots$silhouette

model$ModelInfo$plots$centre_heatmap

projected <- ProjectCluster(model, df_Projection)
projected$ProjectionFit$plots$projection_fit_class_bar

# }
```
