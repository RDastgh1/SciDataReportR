# Compare projectable clustering methods with common participant subsamples

Fits each requested candidate on the same full reference cohort, then on
identical 90% participant subsamples. Every refit is projected to the
full complete reference cohort before participant-level agreement is
calculated. This makes ARI comparable across clustering methods rather
than comparing only participants retained by a particular refit.

## Usage

``` r
CompareClusterMethodStability(
  data,
  variables,
  candidates = NULL,
  stability_resamples = 10L,
  stability_seed = 934522L,
  stability_cores = NULL,
  subsample_fraction = 0.9
)
```

## Arguments

- data:

  A data frame containing the clustering variables.

- variables:

  Numeric clustering variables.

- candidates:

  A candidate specification data frame. `NULL` constructs the planned
  SOM/Mclust, Mclust, K-means, HDBSCAN, and Ward-plus-random-forest
  candidate set. See
  [`DefaultClusterBenchmarkCandidates()`](https://rdastgh1.github.io/SciDataReportR/reference/DefaultClusterBenchmarkCandidates.md).

- stability_resamples:

  Number of common 90% subsample refits. Defaults to 10 for screening.

- stability_seed:

  Seed used to generate the shared participant draws.

- stability_cores:

  Optional number of parallel workers.

- subsample_fraction:

  Proportion of complete reference participants retained in each refit.

## Value

A `ClusterMethodStabilityBenchmark` object with `decision_table`,
`replicates`, shared `subsamples`, `participant_inclusion`, and
full-data `models`.

## Details

Candidate rows use `Method` values `"som_mclust"`, `"mclust"`,
`"kmeans"`, `"hdbscan"`, or `"ward_rf"`. K-based candidates require
`Classes`; mixture candidates also require `Model`. SOM candidates use
`NodeWeighting`, `SOM_XDim`, `SOM_YDim`, and `SOM_Topology`; HDBSCAN
uses `MinPts` and `Epsilon`; and Ward-plus-RF uses `WardN` and
`RF_Trees`.
