# Validate a merged dataset against its source datasets

Audit a merged dataset against the two source datasets used to create
it. This function is designed to catch common merge problems before
analysis, including incompatible key types, duplicate key combinations,
incomplete key coverage, unexpected overlapping variables, unresolved
`.x` / `.y` variables, and value conflicts between duplicated variables.

## Usage

``` r
ValidateMerge(
  LeftData,
  RightData,
  MergedData,
  keys,
  Keys = lifecycle::deprecated()
)
```

## Arguments

- LeftData:

  A data frame used as one source for the merge.

- RightData:

  A data frame used as the other source for the merge.

- MergedData:

  The merged data frame to audit.

- keys:

  Character vector of key variables intended to define the merge.
  Multiple keys are supported, such as `c("study_id", "visit")`.

- Keys:

  Deprecated. Use `keys` instead.

## Value

A list with merge validation results.

## Details

This function does not perform a merge. It assumes the user has already
created a merged dataset and wants to check whether the merge can be
trusted.
