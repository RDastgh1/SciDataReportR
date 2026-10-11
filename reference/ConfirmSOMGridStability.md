# Confirm selected SOM grid stability candidates

Refit exact user-shortlisted rows from
`ScreenSOMGridStability(screen_scope = "joint")` with more 90% subsample
refits. For the default map-only screen, choose and freeze a map first,
then compare k with
[`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md).
This function does not select candidates: every requested joint-screen
row is returned as a separate, projectable SOM + Mclust model.

## Usage

``` r
ConfirmSOMGridStability(
  screen,
  candidates,
  confirm_resamples = 30L,
  distance_percentile = 0.95,
  inclusion_probability = 0.6
)
```

## Arguments

- screen:

  A `SOMGridStabilityScreen` object.

- candidates:

  Values of `.ScreenRow`, or rows from the screen decision table
  containing that column.

- confirm_resamples:

  Number of stability refits, default 30.

- distance_percentile:

  Distance percentile threshold for participant influence flags.

- inclusion_probability:

  Inclusion-probability threshold for influence flags.

## Value

A `SOMGridStabilityConfirmation` object with a decision table and named,
projectable candidate models.
