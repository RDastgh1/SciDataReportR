# Fit PCA followed by K-means

Best for high-dimensional correlated continuous measures with compact
clusters in PCA score space.

## Usage

``` r
CreateClusterModel_PCA_KMeans(
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
  pca_variance_threshold = 0.85,
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

- k_range:

  Candidate cluster counts in exploratory mode.

- final_k:

  Number of clusters for a finalized PCA + K-means solution.

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

- pca_variance_threshold:

  Cumulative variance retained by the existing PCA workflow.

- stability_resamples:

  Number of 80% participant subsample refits used to estimate candidate
  reproducibility. Subsamples are drawn without replacement. Use `0` to
  disable stability analysis.

- stability_seed:

  Seed controlling participant subsampling.

- stability_progress:

  Whether to print subsample progress messages.

## Value

A frozen PCA and K-means pipeline with full-pipeline subsample
stability. `ModelInfo$plots` leads with the reduction layer's `scree`
and `loadings`, followed by the K-means structure figures in score space
(including the per-participant `silhouette` profile) and `profiles`
restated in the original measurement scale. WSS/BSS, silhouette, and
Calinski-Harabasz use the same interpretation as K-means.

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

Jolliffe IT, Cadima J. *Phil Trans R Soc A.* 2016;374:20150202.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
vars_Numeric <- paste0("Var", 1:12)
review <- CreateClusterModel_PCA_KMeans(
  df_Training, vars_Numeric, k_range = 2:5, nstart = 10,
  stability_resamples = 2
)
review$ModelInfo$fit_table
#> # A tibble: 4 × 23
#>   Classes   WSS BetweenSS CalinskiHarabasz Silhouette MinClusterN SizeBalance
#>     <int> <dbl>     <dbl>            <dbl>      <dbl>       <int>       <dbl>
#> 1       2  977.      299.             97.5      0.282          80       0.333
#> 2       3  679.      597.            139.       0.405          80       0.5  
#> 3       4  383.      893.            246.       0.543          80       1    
#> 4       5  328.      948.            227.       0.513          24       0.3  
#> # ℹ 16 more variables: Silhouette_scaled <dbl>, CalinskiHarabasz_scaled <dbl>,
#> #   MinClusterN_scaled <dbl>, SizeBalance_scaled <dbl>, WSS_scaled <dbl>,
#> #   ahp_index <dbl>, Recommended <lgl>, StabilitySuccessRate <dbl>,
#> #   StabilityARI_Mean <dbl>, StabilityARI_P05 <dbl>,
#> #   StabilityJaccard_Mean <dbl>, StabilityJaccard_Min <dbl>,
#> #   NoiseSensitivity <dbl>, NoiseSpecificity <dbl>, ReproducibilityScore <dbl>,
#> #   ReproducibilityScore_scaled <dbl>
review$ModelInfo$AHP$recommendation
#> [1] "AHP-style review recommends PCA plus K-means k (Classes = 4). Review this advisory choice alongside the candidate plots."
review$fit_plot


# The reduction layer comes first: how many components were retained, and
# what each one is actually made of.
review$ModelInfo$plots$scree

review$ModelInfo$plots$loadings


# Then the clustering in that score space.
review$ModelInfo$plots$silhouette

# }
```
