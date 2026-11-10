# Fit hierarchical clusters on a frozen self-organizing map

Trains the shared SOM workflow and clusters all node codebook profiles
with equal weight. Participants inherit their best-matching node's
dendrogram label. New participants map to the frozen SOM without
refitting.

## Usage

``` r
CreateClusterModel_SOM_Hierarchical(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  linkage = "ward.D2",
  distance = "euclidean",
  ClusterVariableName = "Cluster",
  seed_som = 934521L,
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

- linkage:

  Agglomeration linkage: `"ward.D2"` (default), `"average"`,
  `"complete"`, or `"single"`.

- distance:

  `"euclidean"` (default) or `"manhattan"`. Ward linkage requires
  Euclidean distance.

- ClusterVariableName:

  Name of the appended cluster column.

- seed_som:

  Seed for SOM training.

- stability_resamples:

  Number of full-pipeline 90% participant subsample refits.

- seed_stability:

  Seed for participant subsampling.

- stability_progress:

  Whether to print stability progress.

- stability_cores:

  Worker count, or `NULL` for the shared memory-aware policy.

- ...:

  Additional SOM controls passed to
  [`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md),
  including `som_xdim`, `som_ydim`, `som_epochs`, `ZScoreType`, and
  `Relabel`. Assignment is always hard and node weighting is always
  uniform.

## Value

A `Pipeline_SOM_Hierarchical` object with the shared SOM outputs,
`ModelInfo_Hierarchical`, and `ModelInfo` containing the node
dendrogram, candidate fit table, medoids, profiles, silhouette, and
cluster sizes. Candidate `MinClusterN` counts nodes; `MinParticipantN`
counts people. `cluster_sizes` reports both for the selected cut,
including empty nodes. Finalized objects retain `DataWithClusters`,
`ProbFit`, and SOM distance diagnostics. Probability fields inherited
from the SOM infrastructure are missing because hierarchical clustering
has no posterior probabilities.

## Details

Node occupancy is reported, not used to weight the dendrogram. Stability
refits scaling, the SOM, and the dendrogram on 90% participant
subsamples using the reference grid dimensions. It compares sampled
participant assignments. Stability accompanies the initial candidate
recommendation without changing the selected node labels after map
assembly. Projection retains frozen-map distance percentiles,
high-distance flags, range checks, residuals, and cohort summaries. Node
medoid margins describe node separation, not participant probabilities.

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
model <- CreateClusterModel_SOM_Hierarchical(
  SimulatedPhenotypeData, paste0("Var", 1:6),
  method = "finalize", final_k = 3, som_xdim = 4, som_ydim = 4)
#> Warning: Explicit SOM dimensions (4 x 4 = 16 nodes) are below the floor of 10 nodes per clustering variable (60 nodes for 6 variables); diagnostics are retained for review.
projected <- ProjectCluster(model, SimulatedPhenotypeData)
#> far-from-map people, no group found; no far-from-map excess
projected$ProjectionFit$summary
#> # A tibble: 18 × 2
#>    metric                                     value
#>    <chr>                                      <dbl>
#>  1 n_total                                 480     
#>  2 n_complete                              480     
#>  3 n_missing                                 0     
#>  4 training_mean_som_distance                0.677 
#>  5 projected_mean_som_distance               0.655 
#>  6 training_median_som_distance              0.670 
#>  7 projected_median_som_distance             0.649 
#>  8 mean_distance_ratio                       0.967 
#>  9 high_distance_burden                      0.0312
#> 10 cluster_high_distance_burden              0.0312
#> 11 phenotype_drift_index                     0.0334
#> 12 excess_high_distance                     -0.0188
#> 13 cluster_occupancy_js_divergence           0     
#> 14 node_occupancy_js_divergence              0     
#> 15 one_person_node_share                     0     
#> 16 rough_expected_rate_unfamiliar_profiles   0     
#> 17 chance_flag_rate                          0.05  
#> 18 far_from_map_share                        0.0333
# }
```
