# Fit a projectable latent class model for categorical measures

Use latent class analysis when the clustering variables are nominal or
ordinal questionnaire, symptom, diagnosis, or assay-call items. It
estimates class-specific response probabilities and posterior
membership.

## Usage

``` r
CreateClusterModel_LatentClass(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ClusterVariableName = "Cluster",
  nrep = 20L,
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

## Arguments

- data:

  Data frame containing categorical variables.

- variables:

  Variables included in the latent class model.

- method:

  Either `"exploratory"`/`"explore"` or `"finalize"`.

- k_range:

  Candidate class counts.

- final_k:

  Optional class count; when supplied, fit only this solution.

- ClusterVariableName:

  Output class column name.

- nrep:

  Number of random starts per candidate.

- seed:

  Random seed controlling latent-class random starts.

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

A fitted latent-class model with BIC, AIC, log likelihood, entropy,
class-size and subsample stability metrics in `ModelInfo$fit_table`. Log
likelihood and entropy are higher-is-better; AIC and BIC are
lower-is-better. Figures sit beside what they describe: `fit_plot`
reviews candidates; `ModelInfo$plots` holds `response_probabilities`,
`item_profiles`, `bic`, `entropy`, `posterior_map`, and
`categorical_composition`; `ProbFit$plots` holds posterior-confidence
figures; and `Stability$plots` holds cluster-recovery and complementary
recovery.

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

Lazarsfeld PF, Henry NW. *Latent Structure Analysis*. Houghton Mifflin;
1968. Linzer DA, Lewis JB. poLCA: An R package for polytomous variable
latent class analysis. *J Stat Softw.* 2011;42(10):1-29.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
vars_Categorical <- paste0("CatVar", 1:3)
review <- CreateClusterModel_LatentClass(
  df_Training, vars_Categorical, k_range = 2:4, nrep = 5,
  stability_resamples = 2
)
review$ModelInfo$fit_table
#> # A tibble: 3 × 23
#>   Classes   BIC   AIC LogLik Entropy MinClassN SizeBalance StabilitySuccessRate
#>     <int> <dbl> <dbl>  <dbl>   <dbl>     <int>       <dbl>                <dbl>
#> 1       4 1328. 1181.  -552.   0.999        80         1                      1
#> 2       3 1643. 1534.  -738.   0.999        80         0.5                    1
#> 3       2 1983. 1911.  -937.   0.999       160         1                      1
#> # ℹ 15 more variables: StabilityARI_Mean <dbl>, StabilityARI_P05 <dbl>,
#> #   StabilityJaccard_Mean <dbl>, StabilityJaccard_Min <dbl>,
#> #   NoiseSensitivity <dbl>, NoiseSpecificity <dbl>, ReproducibilityScore <dbl>,
#> #   Entropy_scaled <dbl>, MinClassN_scaled <dbl>, SizeBalance_scaled <dbl>,
#> #   ReproducibilityScore_scaled <dbl>, BIC_scaled <dbl>, AIC_scaled <dbl>,
#> #   ahp_index <dbl>, Recommended <lgl>
review$ModelInfo$AHP$recommendation
#> [1] "AHP-style review recommends latent class count (Classes = 4). Review this advisory choice alongside the candidate plots."
review$fit_plot

review$ModelInfo$plots$bic

model <- CreateClusterModel_LatentClass(
  df_Training, vars_Categorical, method = "finalize", final_k = 4,
  nrep = 5, stability_resamples = 2
)
model$ModelInfo$plots$response_probabilities

model$ModelInfo$plots$item_profiles

model$ProbFit$plots$confidence_density

projected <- ProjectCluster(model, df_Projection)
projected$ProjectionFit$plots$projection_fit_class_bar

# }
```
