# Interpret stability output from a finalized clustering model

This shared section is included in the reference documentation for every
finalized clustering constructor. `Stability` is `NULL` when
`stability_resamples = 0`; the fields below are returned when subsample
stability analysis is enabled.

## Details

Shared clustering stability-output reference

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

Hubert L, Arabie P. Comparing partitions. *Journal of Classification*.
1985;2:193-218. Jaccard P. Nouvelles recherches sur la distribution
florale. *Bulletin de la Societe Vaudoise des Sciences Naturelles*.
1901;37:547-579. Meila M. Comparing clusterings: an axiomatic view.
*Proceedings of ICML*. 2005:577-584. Strehl A, Ghosh J. Cluster
ensembles—a knowledge reuse framework for combining multiple partitions.
*Journal of Machine Learning Research*. 2002;3:583-617. Fowlkes EB,
Mallows CL. A method for comparing two hierarchical clusterings.
*Journal of the American Statistical Association*. 1983;78:553-569.
Monti S, Tamayo P, Mesirov J, Golub T. Consensus clustering: a
resampling-based method for class discovery and visualization of gene
expression microarray data. *Machine Learning*. 2003;52:91-118.
