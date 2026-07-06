# Merge QC workflow

Merging datasets is where many analysis errors are born: silently
duplicated rows, IDs that never match, and `.x` / `.y` column pairs that
nobody resolves. SciDataReportR wraps the merge-and-audit cycle into a
small workflow:

- [`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
  performs a left join (exact-key or closest-time) and immediately
  audits the result with
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).
- [`merge_detail()`](https://rdastgh1.github.io/SciDataReportR/reference/merge_detail.md)
  prints plain-text diagnostic tables for one merge.
- [`ExploreMergeValidation()`](https://rdastgh1.github.io/SciDataReportR/reference/ExploreMergeValidation.md)
  renders an interactive dashboard for one merge.
- [`merge_summary_table()`](https://rdastgh1.github.io/SciDataReportR/reference/merge_summary_table.md)
  stacks the one-row logs from many merges into a pipeline-level rollup.

``` r

library(SciDataReportR)
library(survival)
```

## Example 1: a real longitudinal merge (survival::pbc + pbcseq)

The `survival` package ships two related datasets from the Mayo Clinic
primary biliary cholangitis trial:

- `pbc`: baseline data, one row per patient (418 patients, including 106
  who did not participate in the randomized trial),
- `pbcseq`: longitudinal follow-up labs, multiple rows per patient (only
  the 312 trial participants).

Merging them by patient `id` is a completely ordinary task — and it
surfaces real coverage gaps and duplicate-key situations that were not
staged for this vignette.

``` r

pbc_baseline <- pbc[, c("id", "age", "sex", "trt", "stage")]
pbc_labs <- pbcseq[, c("id", "day", "bili", "albumin", "platelet")]

m_pbc <- safe_merge(
  pbc_baseline,
  pbc_labs,
  by = "id",
  name = "pbc baseline + follow-up labs"
)
```

[`safe_merge()`](https://rdastgh1.github.io/SciDataReportR/reference/safe_merge.md)
returns four elements and prints nothing on its own. The one-row `log`
is the pipeline-friendly view:

``` r

m_pbc$log
#> # A tibble: 1 × 16
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 pbc baselin… FAIL   FALSE                   418      2051          5         9
#> # ℹ 9 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   MatchedKeys <int>, LeftUniqueKeys <int>, MatchRate <dbl>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>, KeyHarmonization <chr>,
#> #   Note <chr>
```

The `summary` element is a ready-made
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html):

``` r

m_pbc$summary
```

| Metric | Value |
|:---|:---|
| Status | \<span style=” font-weight: bold; color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(198, 40, 40, 255) !important;” \>FAIL\</span\> |
| Rows (before -\> after) | 418 -\> 2051 |
| Columns (before -\> after) | 5 -\> 9 |
| Columns added (expected vs actual) | expected +4; actual +4 |
| Keys matched | 312 / 418 |
| Match rate | 74.6% |
| Duplicate key groups | 570 |
| Key harmonization | Key types already compatible. |
| Note | Rows changed after merge. Review for row multiplication or filtering. |

pbc baseline + follow-up labs {.table .table
style="width: auto !important; margin-left: auto; margin-right: auto;"}

The audit correctly reports this merge as `FAIL` — not because the code
is wrong, but because the result needs a decision from you before
analysis:

- **Coverage**: 106 patients exist only in the baseline data (the
  non-trial patients who have no follow-up labs).
- **Duplicate keys**: `pbcseq` has many rows per patient, so `id` alone
  does not uniquely identify a row after the merge. If you intended a
  one-row-per- patient dataset, this is a genuine error; if you intended
  a longitudinal dataset, you now have written confirmation of what
  happened.

### Drilling in with merge_detail()

[`merge_detail()`](https://rdastgh1.github.io/SciDataReportR/reference/merge_detail.md)
prints plain kable tables for the checks, the unmatched keys on each
side, overlapping variables, and suspicious conflicts, skipping any
empty section. In R Markdown, use a chunk with `results = "asis"`:

``` r

merge_detail(m_pbc, TopN = 5)
```

| Check | Count | Status | Details |
|:---|---:|:---|:---|
| Key Types | 0.000 | PASS | Key storage classes match across datasets. |
| Missing Keys | 0.000 | PASS | No missing key rows detected. |
| Duplicate Keys | 570.000 | FAIL | Duplicate complete key combinations were detected. |
| Coverage | 106.000 | WARNING | Some complete key combinations appear only in one source dataset. |
| Row Inflation | 1.054 | WARNING | MergedData has more rows than expected. Review whether row multiplication was intentional. |
| Overlapping Variables | 0.000 | PASS | No non-key variables overlap across source datasets. |
| Unresolved Duplicate Variables | 0.000 | PASS | No unresolved duplicate variable pairs detected. |
| Variable Conflicts | 0.000 | PASS | No duplicated-variable value conflicts detected. |
| Suspicious Conflicts | 0.000 | PASS | No low-agreement or class-mismatched duplicated variables detected. |
| Merge Readiness | 1.000 | FAIL | Major merge-integrity blockers detected. Review duplicate keys and unresolved duplicate variables. |

pbc baseline + follow-up labs: validation checks {.table}

| id  |
|:----|
| 313 |
| 314 |
| 315 |
| 316 |
| 317 |

pbc baseline + follow-up labs: keys only in left data (showing up to 5
of 106) {.table}

### A closest-time variant

Suppose instead you want a one-row-per-patient dataset containing the
lab draw closest to enrollment. `safe_merge(method = "closest_time")`
calls
[`Merge_ByClosestTime()`](https://rdastgh1.github.io/SciDataReportR/reference/Merge_ByClosestTime.md)
under the hood, matching on `id` exactly and picking the `pbc_labs` row
whose `day` is nearest to each baseline row’s time (enrollment is day 0
here).

``` r

pbc_enroll <- pbc_baseline
pbc_enroll$enroll_day <- 0

m_closest <- safe_merge(
  pbc_enroll,
  pbc_labs,
  by = "id",
  name = "closest lab to enrollment",
  method = "closest_time",
  time_var_before = "enroll_day",
  time_var_add = "day"
)

m_closest$log
#> # A tibble: 1 × 16
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 closest lab… FAIL   FALSE                   418       418          6        10
#> # ℹ 9 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   MatchedKeys <int>, LeftUniqueKeys <int>, MatchRate <dbl>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>, KeyHarmonization <chr>,
#> #   Note <chr>
```

Note the caveat:
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
audits duplicate keys using `by` alone (here `id`), not the time
variables. In longitudinal settings, repeated keys in the lab data are
legitimate repeated visits, so a duplicate-key `FAIL` from a
closest-time merge may reflect expected repetition rather than a merge
error. The validation output is deliberately left as-is — read it with
that caveat in mind (here, the merged data itself has one row per
patient, which you can confirm from `RowsBefore` and `RowsAfter`).

## Example 2: a synthetic merge that fails every check

To demonstrate the remaining check types — key-type mismatches,
unresolved `.x` / `.y` pairs, and value conflicts — here is a small
synthetic pair of data frames constructed to go wrong in all the classic
ways: a duplicated key, an ID stored as integer on one side and double
on the other, IDs missing from each side, and a `sex` column present in
both sources with disagreeing values.

``` r

demographics <- data.frame(
  id = 1:6, # integer
  sex = c("F", "M", "F", "F", "M", "F"),
  age = c(54, 61, 47, 66, 58, 50)
)

device_data <- data.frame(
  id = c(1, 2, 2, 5, 7), # double, with id 2 duplicated
  sex = c("F", "F", "M", "M", "F"), # disagrees with demographics for id 2/5
  score = c(0.82, 0.75, 0.71, 0.64, 0.90)
)

m_synth <- safe_merge(
  demographics,
  device_data,
  by = "id",
  name = "demographics + device data"
)

m_synth$log
#> # A tibble: 1 × 16
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 demographic… FAIL   FALSE                     6         7          3         5
#> # ℹ 9 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   MatchedKeys <int>, LeftUniqueKeys <int>, MatchRate <dbl>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>, KeyHarmonization <chr>,
#> #   Note <chr>
```

``` r

m_synth$summary
```

| Metric | Value |
|:---|:---|
| Status | \<span style=” font-weight: bold; color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: rgba(198, 40, 40, 255) !important;” \>FAIL\</span\> |
| Rows (before -\> after) | 6 -\> 7 |
| Columns (before -\> after) | 3 -\> 5 |
| Columns added (expected vs actual) | expected +2; actual +2 |
| Keys matched | 3 / 6 |
| Match rate | 50% |
| Duplicate key groups | 2 |
| Key harmonization | id: integer / numeric -\> numeric |
| Note | Rows changed after merge. Review for row multiplication or filtering. |

demographics + device data {.table .table
style="width: auto !important; margin-left: auto; margin-right: auto;"}

Every check type in
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
now has something to say: the duplicate key and the unresolved `sex.x` /
`sex.y` pair are integrity blockers (`FAIL`), while key-type coercion,
coverage gaps, row inflation, overlapping variables, and value conflicts
appear as warnings.

``` r

merge_detail(m_synth)
```

| Check | Count | Status | Details |
|:---|---:|:---|:---|
| Key Types | 0.000 | PASS | Key storage classes match across datasets. |
| Missing Keys | 0.000 | PASS | No missing key rows detected. |
| Duplicate Keys | 2.000 | FAIL | Duplicate complete key combinations were detected. |
| Coverage | 4.000 | WARNING | Some complete key combinations appear only in one source dataset. |
| Row Inflation | 1.167 | WARNING | MergedData has more rows than expected. Review whether row multiplication was intentional. |
| Overlapping Variables | 1.000 | WARNING | Variables appear in both source datasets but were not specified as keys. |
| Unresolved Duplicate Variables | 1.000 | FAIL | MergedData still contains unresolved .x/.y or \_x/\_y variable pairs. |
| Variable Conflicts | 4.000 | WARNING | At least one duplicated variable pair contains conflicting values. |
| Suspicious Conflicts | 1.000 | WARNING | At least one duplicated variable has low agreement or mismatched classes. |
| Merge Readiness | 1.000 | FAIL | Major merge-integrity blockers detected. Review duplicate keys and unresolved duplicate variables. |

demographics + device data: validation checks {.table}

| id  |
|:----|
| 3   |
| 4   |
| 6   |

demographics + device data: keys only in left data (showing up to 10 of
3) {.table}

| id  |
|:----|
| 7   |

demographics + device data: keys only in right data (showing up to 10 of
1) {.table}

| Variable |
|:---------|
| sex      |

demographics + device data: overlapping non-key variables {.table}

| Variable | XVariable | YVariable | LeftClass | RightClass | Agreement | Conflicts | MissingnessConflicts | BothMissing | TotalRows |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|
| sex | sex.x | sex.y | character | character | 42.86 | 4 | 3 | 0 | 7 |

demographics + device data: suspicious duplicate-variable conflicts
{.table style="width:100%;"}

### Interactive review with ExploreMergeValidation()

For interactive QC sessions, pass the `validation` element to
[`ExploreMergeValidation()`](https://rdastgh1.github.io/SciDataReportR/reference/ExploreMergeValidation.md).
The default `Detail = "Compact"` keeps the checks table front and
center, with the coverage and conflict explorers as click-to-expand
accordion sections labeled with their item counts. Use `Detail = "Full"`
to render them expanded.

``` r

ExploreMergeValidation(
  m_synth$validation,
  Title = "Demographics + device data",
  Detail = "Compact"
)
```

Demographics + device data

Interactive review of merge integrity from ValidateMerge().

Validation checks

Search, filter, sort, and expand rows to inspect merge-integrity
examples. INFO rows summarize rows, columns, and unique keys across the
source and merged datasets.

Coverage explorer (4 unmatched)

Review matching, left-only, and right-only key combinations.

Duplicate-variable conflicts (1 variable)

Expand a variable to review conflicting .x and .y values side by side.

Suggested actions

Recommended next steps generated from the merge audit.

## Pipeline rollup with merge_summary_table()

After a sequence of merges, stack the logs into a single table. With
`flagged_only = TRUE`, only merges whose worst check status is not
`PASS` remain — an end-of-pipeline QC gate:

``` r

merge_summary_table(
  list(m_pbc$log, m_closest$log, m_synth$log),
  flagged_only = TRUE
)
#> # A tibble: 3 × 16
#>   Merge        Status ReadyForAnalysis RowsBefore RowsAfter ColsBefore ColsAfter
#>   <chr>        <chr>  <lgl>                 <int>     <int>      <int>     <int>
#> 1 pbc baselin… FAIL   FALSE                   418      2051          5         9
#> 2 closest lab… FAIL   FALSE                   418       418          6        10
#> 3 demographic… FAIL   FALSE                     6         7          3         5
#> # ℹ 9 more variables: ExpectedColsAdded <int>, ActualColsAdded <int>,
#> #   MatchedKeys <int>, LeftUniqueKeys <int>, MatchRate <dbl>,
#> #   DuplicateKeyGroups <int>, UnresolvedDupVars <int>, KeyHarmonization <chr>,
#> #   Note <chr>
```

An empty table here would mean every merge in the pipeline passed
cleanly; any row that appears is a merge that still needs review before
analysis.
