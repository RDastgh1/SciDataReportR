# Fit covariance-aware K-means clusters on a frozen SOM

Trains the established SOM, then fits K-means to equally weighted,
whitened node profiles. The frozen regularized node covariance defines a
shared Mahalanobis metric for node clustering. SOM training and
best-matching-node assignment retain their original geometry.

## Usage

``` r
CreateClusterModel_SOM_KMeans(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  nstart = 50L,
  covariance_shrinkage = 0.05,
  ClusterVariableName = "Cluster",
  seed_som = 934521L,
  seed_kmeans = 93421L,
  stability_resamples = 0L,
  seed_stability = seed_som + 1L,
  stability_progress = FALSE,
  stability_cores = NULL,
  ...
)
```

## Arguments

- data:

  Data frame containing numeric clustering variables.

- variables:

  Character vector of clustering variables.

- method:

  Either `"exploratory"` or `"finalize"`.

- k_range:

  Candidate cluster counts in exploratory mode.

- final_k:

  Explicit cluster count in finalize mode.

- nstart:

  Number of K-means random starts; defaults to 50.

- covariance_shrinkage:

  Number from zero to one, default 0.05. Covariance is shrunk toward a
  spherical matrix with the same average variance. Zero requires
  positive-definite node covariance.

- ClusterVariableName:

  Name of the appended cluster column.

- seed_som:

  Seed for SOM training.

- seed_kmeans:

  Seed for node K-means fits; candidate k uses seed_kmeans + k.

- stability_resamples:

  Number of full-pipeline 90% participant subsample refits.

- seed_stability:

  Seed for participant subsampling.

- stability_progress:

  Whether to print stability progress.

- stability_cores:

  Worker count, or `NULL` for the shared memory-aware policy.

- ...:

  SOM preprocessing and map controls passed to
  [`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md),
  including `ZScoreType`, `som_xdim`, `som_ydim`, `som_epochs`, and
  `Relabel`. Assignment is hard and node weighting is uniform; mixture
  and stability controls in `...` are rejected.

## Value

A `Pipeline_SOM_KMeans` object with `ModelInfo`, `ModelInfo_KMeans`,
`ModelInfo_SOM`, `fit_plot`, `Specification`, and `Stability`. Finalized
models include `DataWithClusters` and `ProbFit`. `ModelInfo$Whitening`
freezes the node mean, raw and regularized covariance, and transforms.
`ModelInfo$fit_table` reports whitened-node WSS, silhouette,
Calinski-Harabasz, minimum node and participant counts, and size
balance. `ModelInfo$centers` and profile figures use original SOM
feature units; `kmeans_model$centers` uses whitened coordinates.
`cluster_sizes` reports node and participant counts, including empty
nodes.

## Details

Participants inherit their frozen best-matching node's cluster.
`DistanceToCentroid` measures participant distance to that cluster's
centroid in the frozen Mahalanobis metric. `AssignmentMargin` is
distance to the closest alternative centroid minus assigned-centroid
distance; it can be negative. `Flag_ProjectionDisagreement` compares
inherited assignment to nearest-centroid assignment without changing the
label. `Node*` fields describe the corresponding node profile. These
measures are not membership probabilities. SOM distance diagnostics
remain in `ModelInfo_SOM$SOMFit`; centroid diagnostics are in
`ModelInfo$FitDiagnostics`. Stability retrains scaling, the SOM,
whitening, and K-means on participant subsamples with reference grid
dimensions. It accompanies the initial recommendation without changing
node labels after map assembly. External preprocessing supplied by the
caller remains externally fixed during refits.

## Stability output

Stability assesses internal reproducibility by full-pipeline 90%
participant subsampling without replacement. For each replicate, 90% of
complete participants are selected once, all preprocessing, any
reduction (PCA, MCA, or SOM), and the selected clustering method are
refit. Most pipelines compare assignments for sampled participants.
SOM + Mclust instead maps the complete reference cohort onto every refit
and compares the resulting full participant partition.
`Stability$settings$comparison_scope` records the policy used. It is an
internal sensitivity analysis, not independent-cohort validation.

`Stability$settings` records the analysis provenance:

- `resamples`: requested number of 90% subsample refits.

- `seed`: seed used to derive reproducible sampling and model seeds.

- `refit_scope`: whether agreement is calculated directly in the
  refitted sample or after refit projection onto the reference cohort.

- `comparison_scope`: the participants entering replicate agreement
  metrics; most pipelines use `"sampled_participants"`, while SOM +
  Mclust uses `"all_reference_participants"`.

- `resample_type`: `"subsample_without_replacement"` for the primary
  stability analysis.

- `resample_fraction`: the retained participant fraction, `0.90`.

- `resample_fractions`: every requested retention fraction. The default
  is `0.90`; additional values create an opt-in deletion-stability curve
  without replacing the 90% primary screen.

- `coassignment_limit`: maximum number of complete training participants
  (2,000) for which the full pairwise co-assignment matrix is
  calculated.

- `noise_policy`: whether the method has noise labels. For ordinary
  methods it is `"all clusters included"`; HDBSCAN variants retain noise
  in ARI but exclude it from phenotype-level Jaccard recovery.

- `requested_cores`, `resolved_cores`, and `worker_policy`: the
  requested and resolved refit workers and the reason for the
  resolution. Automatic resolution is capped by CPU, refit count, and
  the smaller of 50% physical RAM and 75% currently available RAM.
  Pairwise-memory methods are run serially automatically.

- `estimated_bytes_per_worker`, `estimated_total_bytes`,
  `physical_memory_bytes`, `available_memory_bytes`, and
  `safe_memory_budget_bytes`: the conservative memory audit used to
  resolve automatic workers. An explicit worker request is preserved
  (within the CPU/scheduler cap) but warns if it exceeds the estimated
  safe budget.

`Stability$replicates` has one row per requested refit. `Model` and
`Classes` identify the selected candidate (for HDBSCAN, `Classes` is the
data-derived extracted count); `Replicate` is its sequence number;
`Status` is `"success"` or a failure status; `SamplingSeed` and
`ModelSeed` identify the participant draw and refit initialization
policy; and `Error` contains the error message for an unsuccessful
refit. `ARI` is the primary successful-refit metric: it is invariant to
cluster numbering, higher is better, and it can be negative when
agreement is worse than chance.

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
replicates. `ReproducibilityScore` does not include success rate. For
most pipelines it is the mean of finite `StabilityARI_Mean` and
`StabilityJaccard_Mean` values; SOM + Mclust uses `StabilityARI_Mean`
directly because it selects on agreement of the final participant
partition.

SOM + Mclust additionally returns `StabilityARI_NonBoundary`,
`StabilityARI_NonBoundary_P05`, and non-boundary Jaccard summaries.
These are secondary diagnostics calculated only on the fixed
reference-model mask with `SoftMargin` below its threshold or
`SoftEntropy` above its threshold excluded. The primary all-participant
ARI always retains these boundary participants. `boundary_summary`
reports their count and proportion overall and within each reference
cluster. Earlier posterior-based `Core` fields are retained for
compatibility as descriptive diagnostics; they are not the
boundary-aware sensitivity measure.

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

SOM + Mclust additionally returns memory-safe final-partition
diagnostics for every cohort size. `co_membership_by_cluster` reports
the proportion of within-reference-cluster participant pairs retained
together (`Cohesion`) and the probability of co-assignment with
outsiders (`CrossMixing`). `co_membership_between_clusters` reports
cross-reference-cluster co-assignment probabilities.
`split_merge_reference` and `split_merge_refit` identify dominant
sources, shares, and effective numbers of groups, while
`profile_recovery` reports matched-cluster correlation and RMSE for mean
clustering-variable profiles on the frozen reference scale. These
explain the nature of membership changes; they do not supersede ARI.

`deletion_curve` summarizes ARI, VI, NMI, and Fowlkes–Mallows agreement
at each requested participant-retention fraction. It is useful for
describing sensitivity to increasing cohort perturbation, not for
selecting the least demanding fraction.

`Stability$plots` contains `cluster_recovery` (per-cluster Jaccard),
`partition_metrics` (the ARI distribution), and `cluster_inclusion`; it
also contains a co-assignment heatmap when the matrix is available.
These diagnostics complement rather than replace ARI and Jaccard: none
can turn a poorly reproducible cluster into a stable phenotype.

Metric sources: Hubert and Arabie (1985) define ARI; Jaccard (1901)
defines the overlap coefficient; and Monti et al. (2003) describe
resampling-based consensus co-assignment.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
model <- CreateClusterModel_SOM_KMeans(
  SimulatedPhenotypeData, paste0("Var", 1:6), method = "finalize",
  final_k = 4, som_xdim = 6, som_ydim = 6)
#> Warning: Explicit SOM dimensions (6 x 6 = 36 nodes) are below the floor of 10 nodes per clustering variable (60 nodes for 6 variables); diagnostics are retained for review.
projected <- ProjectCluster(model, SimulatedPhenotypeData)
#> far-from-map people, no group found; no far-from-map excess
projected$ProjectionFit$summary
#> # A tibble: 18 × 2
#>    metric                                      value
#>    <chr>                                       <dbl>
#>  1 n_total                                 480      
#>  2 n_complete                              480      
#>  3 n_missing                                 0      
#>  4 training_mean_som_distance                0.612  
#>  5 projected_mean_som_distance               0.567  
#>  6 training_median_som_distance              0.590  
#>  7 projected_median_som_distance             0.545  
#>  8 mean_distance_ratio                       0.926  
#>  9 high_distance_burden                      0.0271 
#> 10 cluster_high_distance_burden              0.0292 
#> 11 phenotype_drift_index                     0.0767 
#> 12 excess_high_distance                     -0.0229 
#> 13 cluster_occupancy_js_divergence           0      
#> 14 node_occupancy_js_divergence              0      
#> 15 one_person_node_share                     0.00417
#> 16 rough_expected_rate_unfamiliar_profiles   0.00417
#> 17 chance_flag_rate                          0.05   
#> 18 far_from_map_share                        0.0312 
projected$ProjectionFit$centroid_fit$summary
#> # A tibble: 1 × 9
#>       n mean_distance median_distance high_distance_cutoff prop_high_distance
#>   <int>         <dbl>           <dbl>                <dbl>              <dbl>
#> 1   480          2.14            2.09                 3.44               0.05
#> # ℹ 4 more variables: prop_good_fit <dbl>, prop_uncertain <dbl>,
#> #   prop_poor_fit <dbl>, prop_potential_novel <dbl>
# }
```
