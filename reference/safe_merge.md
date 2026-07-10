# Safely merge two data frames with relationship-aware validation

`safe_merge()` performs a merge, validates its structure, logs merge
metrics, and returns the merged data plus validation results. It is
designed for reproducible database construction pipelines where row
count, column count, key coverage, duplicate keys, expected merge
relationships, and unresolved duplicate variables need to be audited
every time the merge is run.

## Usage

``` r
safe_merge(
  df_before,
  df_add,
  by,
  name,
  method = c("exact", "closest_time"),
  time_var_before = NULL,
  time_var_add = NULL,
  min_match_rate = 0.95,
  harmonize_keys = TRUE,
  key_parser = NULL,
  stop_on_failed_numeric = TRUE,
  expected_relationship = c("one-to-one", "many-to-one", "one-to-many", "many-to-many",
    "auto"),
  fail_on_new_duplicate_variables = TRUE,
  fail_on_inherited_duplicate_variables = FALSE,
  ...
)
```

## Arguments

- df_before:

  Data frame on the left side of the merge.

- df_add:

  Data frame to add to `df_before`.

- by:

  Character vector of merge key column names.

- name:

  Single character label used in the merge log and summary table.

- method:

  Merge method. `"exact"` uses
  [`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html).
  `"closest_time"` uses
  [`Merge_ByClosestTime()`](https://rdastgh1.github.io/SciDataReportR/reference/Merge_ByClosestTime.md).

- time_var_before:

  Required when `method = "closest_time"`. Time variable in `df_before`.

- time_var_add:

  Required when `method = "closest_time"`. Time variable in `df_add`.

- min_match_rate:

  Numeric between 0 and 1. Merges below this left-key match rate are
  marked as `"WARNING"` unless another structural blocker makes them
  `"FAIL"`.

- harmonize_keys:

  Logical. If `TRUE`, keys are harmonized using
  [`HarmonizeMergeKeys()`](https://rdastgh1.github.io/SciDataReportR/reference/HarmonizeMergeKeys.md)
  before merging.

- key_parser:

  Optional parser passed to
  [`HarmonizeMergeKeys()`](https://rdastgh1.github.io/SciDataReportR/reference/HarmonizeMergeKeys.md).

- stop_on_failed_numeric:

  Logical passed to
  [`HarmonizeMergeKeys()`](https://rdastgh1.github.io/SciDataReportR/reference/HarmonizeMergeKeys.md).

- expected_relationship:

  Character. Expected relationship between `df_before` and `df_add`
  under `by`. Defaults to `"one-to-one"` to preserve strict historical
  behavior. One of:

  - `"one-to-one"`: both sides should be unique by `by`.

  - `"many-to-one"`: left side may repeat keys, right side should be
    unique.

  - `"one-to-many"`: left side should be unique, right side may repeat.

  - `"many-to-many"`: both sides may repeat.

  - `"auto"`: infer and report relationship, but do not enforce it.

- fail_on_new_duplicate_variables:

  Logical. If `TRUE`, unresolved duplicate variable pairs introduced by
  the current merge cause the merge status to be `"FAIL"`. Default is
  `TRUE`.

- fail_on_inherited_duplicate_variables:

  Logical. If `TRUE`, unresolved duplicate variable pairs already
  present in `df_before` cause the merge status to be `"FAIL"`. Default
  is `FALSE`.

- ...:

  Additional arguments passed to
  [`Merge_ByClosestTime()`](https://rdastgh1.github.io/SciDataReportR/reference/Merge_ByClosestTime.md)
  when `method = "closest_time"`.

## Value

A list with:

- `data`: Merged data frame.

- `validation`: Full validation object from
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md),
  with additional current-merge-aware duplicate-variable diagnostics.

- `log`: One-row tibble containing merge log metrics.

- `summary`: A
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) summary
  table.

## Details

The default `expected_relationship` is `"one-to-one"` to preserve strict
historical behavior. For longitudinal databases, use `"many-to-one"`
when the left/master data frame has repeated participant IDs because of
multiple visits and the right/add-on data frame has one row per
participant.

Duplicate-variable handling is current-merge aware. If unresolved
duplicate variable pairs such as `sy_x` / `sy_y` or `sy.x` / `sy.y`
already exist in `df_before`, they are classified as inherited duplicate
variables. By default, inherited duplicate variables are reported but do
not cause the current merge to fail. New unresolved duplicate variables
introduced by the current merge are treated as failures by default.

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

m <- safe_merge(
  df_before = left,
  df_add = right,
  by = "record_id",
  name = "Example imaging merge",
  expected_relationship = "many-to-one"
)

m$log
#> # A tibble: 1 × 27
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 Example ima… PASS   TRUE                      4         4          3         4
#> # ℹ 20 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   ExpectedRelationship <chr>, DetectedRelationship <chr>,
#> #   RelationshipMatchesExpected <lgl>, MatchedKeys <int>, LeftUniqueKeys <int>,
#> #   MatchRate <dbl>, DuplicateKeyBlockers <int>, DuplicateKeyGroups_Left <int>,
#> #   DuplicateKeyGroups_Right <int>, DuplicateKeyGroups_Merged <int>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>,
#> #   NewUnresolvedDupVars <int>, InheritedUnresolvedDupVars <int>, …
```
