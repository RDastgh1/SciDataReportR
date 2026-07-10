# Validate a merge between two source data frames and a merged result

`ValidateMerge()` audits a completed merge by comparing the left/source
data, right/add-on data, and merged result. It checks key coverage, key
uniqueness, expected merge relationship, row inflation, overlapping
non-key variables, unresolved `.x`/`.y` or `_x`/`_y` columns, and
duplicated-variable conflicts.

## Usage

``` r
ValidateMerge(
  LeftData,
  RightData,
  MergedData,
  keys,
  Keys = lifecycle::deprecated(),
  expected_relationship = c("one-to-one", "many-to-one", "one-to-many", "many-to-many",
    "auto")
)
```

## Arguments

- LeftData:

  A data frame used as the left side of the merge.

- RightData:

  A data frame used as the right side of the merge.

- MergedData:

  A data frame produced by merging `LeftData` and `RightData`.

- keys:

  Character vector of merge key column names.

- Keys:

  Deprecated. Use `keys`.

- expected_relationship:

  Character. Expected relationship between `LeftData` and `RightData`
  under `keys`. Defaults to `"one-to-one"` to preserve strict historical
  behavior. One of:

  - `"one-to-one"`: left keys and right keys should both be unique.

  - `"many-to-one"`: left keys may repeat, right keys should be unique.
    This is common when merging participant-level data into a
    longitudinal participant-visit master table by `record_id`.

  - `"one-to-many"`: left keys should be unique, right keys may repeat.

  - `"many-to-many"`: both sides may repeat. This is allowed only when
    explicitly expected.

  - `"auto"`: infer and report relationship, but do not fail solely
    because of duplicate keys or relationship type.

## Value

A list containing merge validation summaries, checks, relationship
audits, duplicate-key audits, coverage tables, overlap audits,
duplicated variable audits, conflict tables, and suggested actions.

## Details

This version supports expected merge relationships. The default is
`"one-to-one"` to preserve strict historical behavior. For longitudinal
databases, use `"many-to-one"` when the left/master data frame is a
participant-visit table and the right/add-on data frame is
participant-level.

## Examples

``` r
left <- data.frame(
  record_id = c(1, 1, 2, 2),
  visit_type = c(1, 2, 1, 2),
  age = c(40, 40, 55, 55)
)

right <- data.frame(
  record_id = c(1, 2),
  imaging_score = c(0.4, 0.8)
)

merged <- dplyr::left_join(left, right, by = "record_id")

validation <- ValidateMerge(
  LeftData = left,
  RightData = right,
  MergedData = merged,
  keys = "record_id",
  expected_relationship = "many-to-one"
)

validation$Summary
#> # A tibble: 1 × 33
#>   LeftRows RightRows MergedRows RowInflationFactor LeftColumns RightColumns
#>      <int>     <int>      <int>              <dbl>       <int>        <int>
#> 1        4         2          4                  1           3            2
#> # ℹ 27 more variables: MergedColumns <int>, LeftUniqueKeys <int>,
#> #   RightUniqueKeys <int>, MergedUniqueKeys <int>, MatchingKeys <int>,
#> #   LeftOnlyKeys <int>, RightOnlyKeys <int>, LeftMatchRate <dbl>,
#> #   RightMatchRate <dbl>, MissingKeyRows_Left <int>,
#> #   MissingKeyRows_Right <int>, MissingKeyRows_Merged <int>,
#> #   DuplicateKeyGroups_Left <int>, DuplicateKeyGroups_Right <int>,
#> #   DuplicateKeyGroups_Merged <int>, DuplicateKeyBlockers <int>, …
```
