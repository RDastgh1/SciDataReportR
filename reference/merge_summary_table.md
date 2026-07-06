# Combine safe_merge logs into a single summary table

Bind the one-row `$log` tibbles produced by
[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
into a single table, one row per merge, optionally filtered to merges
that did not pass cleanly. Useful as an end-of-pipeline rollup after a
sequence of merges.

## Usage

``` r
merge_summary_table(merge_log, flagged_only = FALSE)
```

## Arguments

- merge_log:

  A list of `$log` tibbles from
  [`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
  results. Full
  [`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
  result objects are also accepted (the `$log` element is extracted
  automatically), as is a single log tibble.

- flagged_only:

  Logical. If `TRUE`, only rows with `Status != "PASS"` are returned.
  Default is `FALSE`.

## Value

A tibble with one row per merge and the same columns as a
[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
`$log` tibble.

## See also

[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)

## Examples

``` r
baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
vitals <- data.frame(id = 1:4, sbp = c(120, 135, 118, 141))
m1 <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
m2 <- safe_merge(m1$data, vitals, by = "id", name = "+ vitals")
merge_summary_table(list(m1$log, m2$log))
#> # A tibble: 2 × 16
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 Baseline + … WARNI… TRUE                      4         4          2         3
#> 2 + vitals     PASS   TRUE                      4         4          3         4
#> # ℹ 9 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   MatchedKeys <int>, LeftUniqueKeys <int>, MatchRate <dbl>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>, KeyHarmonization <chr>,
#> #   Note <chr>
```
