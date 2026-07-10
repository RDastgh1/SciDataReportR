# Print a plain-text detail report for one safe_merge result

Print static
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) tables (no
plots, no htmlwidgets) for the key diagnostic sections of a
[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
result: the validation checks, unmatched key combinations on each side,
overlapping non-key variables, and suspicious duplicate-variable
conflicts. Sections with nothing to show are skipped entirely.

## Usage

``` r
merge_detail(m, TopN = 10)
```

## Arguments

- m:

  A list returned by
  [`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
  (must contain a `validation` element from
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
  and a `log` tibble).

- TopN:

  Integer. Maximum number of rows shown for the left-only and right-only
  unmatched-key previews. Default is `10`.

## Value

`m`, invisibly. Called for its printed output.

## Details

In R Markdown or Quarto, use a chunk with `results = "asis"` so the
kable markdown renders as tables.

## See also

[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md),
[`merge_summary_table()`](https://rdastgh1.github.io/SciDataReportR/reference/merge_summary_table.md),
[`ExploreMergeValidation()`](https://rdastgh1.github.io/SciDataReportR/reference/ExploreMergeValidation.md)

## Examples

``` r
baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
merge_detail(m)
#> 
#> 
#> Table: Baseline + labs: validation checks
#> 
#> |Check                          | Count|Status  |Details                                                             |
#> |:------------------------------|-----:|:-------|:-------------------------------------------------------------------|
#> |Key Types                      |     0|PASS    |Key storage classes match across datasets.                          |
#> |Missing Keys                   |     0|PASS    |No missing key rows detected.                                       |
#> |Expected Relationship          |     0|PASS    |Detected relationship matches expected_relationship = 'one-to-one'. |
#> |Duplicate Keys                 |     0|PASS    |No duplicate key blockers detected for the expected relationship.   |
#> |Coverage                       |     1|WARNING |Some complete key combinations appear only in one source dataset.   |
#> |Row Count                      |     0|PASS    |MergedData row count matches LeftData row count.                    |
#> |Row Inflation                  |     1|PASS    |No meaningful row inflation detected.                               |
#> |Overlapping Variables          |     0|PASS    |No non-key variables overlap across source datasets.                |
#> |Unresolved Duplicate Variables |     0|PASS    |No unresolved duplicate variable pairs detected.                    |
#> |Variable Conflicts             |     0|PASS    |No duplicated-variable value conflicts detected.                    |
#> |Suspicious Conflicts           |     0|PASS    |No low-agreement or class-mismatched duplicated variables detected. |
#> |Merge Readiness                |     0|PASS    |No major merge-integrity blockers detected.                         |
#> 
#> 
#> 
#> Table: Baseline + labs: keys only in left data (showing up to 10 of 1)
#> 
#> |id |
#> |:--|
#> |3  |
#> 
```
