# Screen fixed SOM map designs before selecting cluster count

In its default `"map"` scope, this function refits only SOM
preprocessing and map training. It selects a representation for later,
fixed-map Mclust comparisons without fitting a mixture or varying `k`.
The legacy joint SOM + Mclust screen remains available with
`screen_scope = "joint"`.

## Usage

``` r
ScreenSOMGridStability(
  data,
  variables = NULL,
  k_range = 2:10,
  models = c(1, 2, 3),
  grid_candidates = NULL,
  topologies = c("hexagonal", "rectangular"),
  node_weighting = c("uniform", "occupancy"),
  min_participants_per_node = 10,
  screen_scope = c("map", "joint"),
  screen_resamples = 3L,
  seed_som = 934521L,
  stability_seed = 934522L,
  stability_cores = NULL,
  pair_sample = 10000L,
  neighborhood_sample = 500L,
  neighborhood_k = 10L,
  ...
)
```

## Arguments

- data:

  A data frame supplied to
  [`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md).

- variables:

  Numeric clustering variables.

- k_range:

  Requested profile counts. In `"map"` scope, the largest requested
  value defines the automatic candidate-resolution ladder but Mclust is
  not fitted and no k is selected.

- models:

  Mclust covariance-model families. In `"map"` scope they are not
  fitted, but their largest requested parameter count supplies a
  candidate-map support floor.

- grid_candidates:

  Optional square side lengths or a data frame with `som_xdim` and
  `som_ydim`. The map-scope default includes approximately 10, 20, and
  30 nodes per maximum requested profile, plus 3x3 and 10x10 anchors and
  a parameter-support floor for the largest requested Mclust family.
  User-supplied candidates are retained unchanged.

- topologies:

  SOM topologies to screen.

- node_weighting:

  Node-mixture definitions. Ignored in `"map"` scope.

- min_participants_per_node:

  Target complete participants per SOM node for automatic candidates.
  Default `10`; anchors and user grids remain in the review table even
  when they do not meet this map-health target.

- screen_scope:

  `"map"` (default) for k-informed, Mclust-free SOM screening or
  `"joint"` for the earlier grid-by-k Mclust workflow.

- screen_resamples:

  Map refits per candidate. The map default is 3; the explicit joint
  workflow retains its prior default of 10.

- seed_som:

  Fixed SOM initialization seed used for every refit of a candidate
  design.

- stability_seed:

  Seed for participant subsampling and metric sampling.

- stability_cores:

  Optional parallel workers for joint-screen refits; map-screen refits
  are intentionally lightweight and serial.

- pair_sample:

  Number of complete-cohort participant pairs used for
  orientation-invariant grid-distance preservation.

- neighborhood_sample:

  Number of participants used for local-neighbor overlap; sampling is
  deterministic from `stability_seed`.

- neighborhood_k:

  Number of local grid neighbors per sampled participant.

- ...:

  Additional arguments passed to
  [`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md).

## Value

A `SOMGridStabilityScreen` object. Map scope returns a decision table
and per-refit `map_replicates`; joint scope returns the prior grid-by-k
Mclust decision table for
[`ConfirmSOMGridStability()`](https://rdastgh1.github.io/SciDataReportR/reference/ConfirmSOMGridStability.md).
