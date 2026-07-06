# Merge two data frames and validate the result in one step

Perform a left join or closest-time merge and immediately audit the
result with
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).
`safe_merge()` can also harmonize merge key types before joining, so
common numeric-vs-character key mismatches do not break an otherwise
valid merge.

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
  ...
)
```

## Arguments

- df_before:

  A data frame. The left dataset that rows are added to.

- df_add:

  A data frame. The right dataset being merged in.

- by:

  Character vector of merge keys. Multiple keys are supported.

- name:

  Required. A short label for this merge.

- method:

  Either `"exact"` or `"closest_time"`.

- time_var_before:

  For `method = "closest_time"`, the time variable in `df_before`.

- time_var_add:

  For `method = "closest_time"`, the time variable in `df_add`.

- min_match_rate:

  Minimum acceptable key match rate before a merge is marked as
  `"WARNING"`.

- harmonize_keys:

  Logical. If `TRUE`, harmonize merge key types before joining.

- key_parser:

  Optional named list of parser functions for specific merge keys. Names
  should match values in `by`.

- stop_on_failed_numeric:

  Logical. If `TRUE`, stop when one key is numeric-like and the other
  contains values that cannot be safely converted to numeric.

- ...:

  Additional arguments passed to
  [`Merge_ByClosestTime()`](https://rdastgh1.github.io/SciDataReportR/reference/Merge_ByClosestTime.md)
  when `method = "closest_time"`.

## Value

A list with four elements:

- data:

  The merged data frame.

- validation:

  The full
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
  result, including `KeyHarmonization` and `SummaryTable`.

- log:

  A one-row tibble with merge-level status and metrics.

- summary:

  A compact formatted summary table.
