# Confirm method-specific stability leaders

Confirm method-specific stability leaders

## Usage

``` r
ConfirmClusterMethodStability(
  benchmark,
  candidates = NULL,
  confirm_resamples = 30L
)
```

## Arguments

- benchmark:

  A `ClusterMethodStabilityBenchmark` object.

- candidates:

  Optional candidate IDs or decision-table rows. If omitted, confirms
  the ARI-leading eligible candidate within each method variant.

- confirm_resamples:

  Number of common refits, default 30.

## Value

A `ClusterMethodStabilityConfirmation` object.
