# Refit an explicit participant-exclusion SOM sensitivity model

Creates a separately labelled sensitivity model after a
reviewer-specified exclusion. It never changes the source model and
records the excluded row IDs and stated reason in the returned object.

## Usage

``` r
RefitSOMSensitivity(
  object,
  data,
  exclude_row_ids,
  exclusion_reason,
  stability_resamples = 30L,
  ...
)
```

## Arguments

- object:

  A confirmed or ordinary `Pipeline_SOM_MClust` object.

- data:

  The original analysis data.

- exclude_row_ids:

  Explicit `.row_id` values to exclude.

- exclusion_reason:

  Human-readable reason recorded in provenance.

- stability_resamples:

  Sensitivity-model refits, default 30.

- ...:

  Additional arguments passed to
  [`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md).

## Value

A separate `Pipeline_SOM_MClust` sensitivity model.
