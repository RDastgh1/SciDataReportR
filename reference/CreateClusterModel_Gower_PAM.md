# Fit a projectable Gower-distance PAM model for mixed clinical data

Best for mixed continuous, binary, ordinal, and nominal clinical
measures where medoid exemplars are more interpretable than centroids.

## Usage

``` r
CreateClusterModel_Gower_PAM(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
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

  Data frame containing numeric, logical, factor, or ordered variables.

- variables:

  Variables used for clustering.

- method:

  Either `"exploratory"`/`"explore"` or `"finalize"`.

- k_range:

  Candidate medoid counts in exploratory mode.

- final_k:

  Final medoid count in finalized mode.

- ClusterVariableName:

  Output cluster column name.

- seed:

  Random seed.

- stability_resamples:

  Number of 90% participant subsample refits.

- stability_seed:

  Seed controlling participant subsampling.

- stability_progress:

  Whether to print subsample progress messages.

- stability_cores:

  Number of workers for stability refits. `NULL` uses all detected
  physical cores minus one, capped at the requested resamples and any
  scheduler limit. Each worker holds a refit in memory, so lower this
  setting for large or high-dimensional analyses.

## Value

A fitted PAM model with frozen medoids, numeric ranges, categorical
levels, silhouette and medoid-distance metrics in `ModelInfo$fit_table`,
and subsample stability. Mean silhouette is higher-is-better and mean
assigned-medoid Gower distance is lower-is-better (a dissimilarity, not
a probability). Figures sit beside what they describe: `fit_plot`
reviews candidates; `ModelInfo$plots` holds `silhouette` (the
per-participant profile from the selected PAM solution),
`silhouette_by_k`, `gower_map`, `categorical_composition`,
`categorical_composition_by_cluster`, `categorical_enrichment`, and
`profiles`; `ModelInfo$FitDiagnostics$plots` holds the medoid-distance
histogram; `ProbFit$plots` holds assignment-margin figures; and
`Stability$plots` holds cluster-recovery and complementary stability
diagnostics.

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

Gower JC. A general coefficient of similarity and some of its
properties. *Biometrics.* 1971;27:857-871. Kaufman L, Rousseeuw PJ.
*Finding Groups in Data*. Wiley; 1990.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
vars_Mixed <- c("Var1", "Var2", "Var3", "CatVar1", "CatVar2", "CatVar3")
review <- CreateClusterModel_Gower_PAM(
  df_Training, vars_Mixed, k_range = 2:5, stability_resamples = 2
)
review$ModelInfo$fit_table
#> # A tibble: 4 × 20
#>   Classes Silhouette MinClusterN SizeBalance MeanMedoidDistance
#>     <int>      <dbl>       <dbl>       <dbl>              <dbl>
#> 1       2      0.581         160      1                  0.291 
#> 2       3      0.716          80      0.503              0.169 
#> 3       4      0.855          80      1                  0.0508
#> 4       5      0.810           3      0.0375             0.0491
#> # ℹ 15 more variables: StabilitySuccessRate <dbl>, StabilityARI_Mean <dbl>,
#> #   StabilityARI_P05 <dbl>, StabilityJaccard_Mean <dbl>,
#> #   StabilityJaccard_Min <dbl>, NoiseSensitivity <dbl>, NoiseSpecificity <dbl>,
#> #   ReproducibilityScore <dbl>, Silhouette_scaled <dbl>,
#> #   MinClusterN_scaled <dbl>, SizeBalance_scaled <dbl>,
#> #   ReproducibilityScore_scaled <dbl>, MeanMedoidDistance_scaled <dbl>,
#> #   ahp_index <dbl>, Recommended <lgl>
review$ModelInfo$AHP$recommendation
#> [1] "AHP-style review recommends Gower/PAM k (Classes = 4). Review this advisory choice alongside the candidate plots."
review$fit_plot

review$ModelInfo$plots$silhouette_by_k

model <- CreateClusterModel_Gower_PAM(
  df_Training, vars_Mixed, method = "finalize", final_k = 4,
  stability_resamples = 2
)
model$ModelInfo$plots$silhouette

model$ModelInfo$plots$gower_map
#> NULL
model$ModelInfo$plots$categorical_composition

model$ModelInfo$medoids
#>          Var1      Var2      Var3   CatVar1   CatVar2   CatVar3
#> 21   2.287020  2.401218  2.555956 Item 1: A Item 2: A Item 3: A
#> 87   2.421441  2.768372  2.604292 Item 1: B Item 2: B Item 3: B
#> 189 -2.285361 -2.838339 -2.431604 Item 1: C Item 2: C Item 3: C
#> 260 -2.304594 -2.759347 -2.201481 Item 1: D Item 2: D Item 3: D
projected <- ProjectCluster(model, df_Projection)
projected$ProjectionFit$plots$projection_fit_class_bar

# }
```
