# Review SOM map support before fitting a SOM + Mclust model

Calculates the competing codebook-node support and participant occupancy
constraints for a planned SOM + Mclust search. Mclust is fit to SOM
codebook vectors, so its effective sample size is the number of nodes;
participant count instead limits how sparse a map should become.

## Usage

``` r
SOMMapSupport(
  n,
  p,
  k_range,
  models,
  min_nodes_per_parameter = 3,
  min_participants_per_node = 10
)
```

## Arguments

- n:

  Number of complete training participants.

- p:

  Number of post-preprocessing clustering variables entering Mclust.

- k_range:

  Requested profile counts.

- models:

  TidyLPA Mclust model IDs (1, 2, 3, or 6).

- min_nodes_per_parameter:

  Target codebook nodes per estimated mixture parameter. Default `3`.

- min_participants_per_node:

  Target complete training participants per SOM node. Default `10`; this
  is an engineering map-health target, not a universal statistical
  threshold.

## Value

A tibble with one row per model/k candidate and common map-sizing
diagnostics. `SupportStatus` is `"supported"`, `"under-supported"`, or
`"over-parameterized"` at the recommended automatic map size.

## Examples

``` r
SOMMapSupport(n = 4330, p = 11, k_range = 2:8, models = 1:3)
#> # A tibble: 21 × 18
#>    Model Classes MixtureParameters ParameterFloorSide MixtureNodes
#>    <int>   <int>             <int>              <dbl>        <int>
#>  1     1       2                34                 11          400
#>  2     1       3                46                 12          400
#>  3     1       4                58                 14          400
#>  4     1       5                70                 15          400
#>  5     1       6                82                 16          400
#>  6     1       7                94                 17          400
#>  7     1       8               106                 18          400
#>  8     2       2                45                 12          400
#>  9     2       3                68                 15          400
#> 10     2       4                91                 17          400
#> # ℹ 11 more rows
#> # ℹ 13 more variables: ParticipantsPerNode <dbl>, OccupancyCompliant <lgl>,
#> #   NodesPerParameter <dbl>, Overparameterized <lgl>, SupportStatus <chr>,
#> #   SOMSide <int>, BaseSide <int>, FloorSide <int>, OccupancyCeilingSide <int>,
#> #   RecommendedSide <int>, BindingConstraint <chr>, MinNodesPerParameter <dbl>,
#> #   MinParticipantsPerNode <dbl>
```
