# Compare two versions of a dataset

Compare an old dataset and a new dataset using one or more key
variables. This function identifies record-level, variable-level, and
cell-level changes between dataset versions. It is useful when reviewing
updated data extracts, revised REDCap exports, cleaned spreadsheet
versions, or vendor-delivered dataset updates.

## Usage

``` r
CompareDatasets(OldData, NewData, keys, Keys = lifecycle::deprecated())
```

## Arguments

- OldData:

  A data frame representing the earlier dataset version.

- NewData:

  A data frame representing the newer dataset version.

- keys:

  Character vector of key variables used to align records across the two
  datasets. Multiple keys are supported, such as
  `c("study_id", "TimePoint")`.

- Keys:

  **Deprecated** (since 19.15.0). Use `keys` instead.

## Value

A list with dataset comparison results, including:

- SummaryText:

  A plain-text summary of the dataset comparison.

- Summary:

  One-row tibble with core comparison metrics.

- Fingerprint:

  Tibble comparing rows, columns, and unique key combinations.

- KeyTypes:

  Tibble showing key variable classes before coercion.

- Checks:

  Tibble summarizing comparison checks and pass/warning/fail status.

- StructureChanges:

  Tibble of variables added to or removed from NewData, using normalized
  variable names.

- AddedVariables:

  Tibble of variables present in NewData but not OldData, using
  normalized variable names.

- RemovedVariables:

  Tibble of variables present in OldData but not NewData, using
  normalized variable names.

- AddedRecords:

  Tibble of key combinations present in NewData but not OldData.

- RemovedRecords:

  Tibble of key combinations present in OldData but not NewData.

- DuplicateKeys:

  List containing duplicated key rows from OldData and NewData.

- ComparisonKeys:

  List describing matching keys, compared keys, and keys excluded from
  cell comparison due to duplicate key combinations.

- NameRepairAudit:

  Tibble describing variables whose raw names differ after removing
  tibble-style `...number` suffixes.

- ComparisonVariableMap:

  Tibble mapping normalized variable names to the raw OldData and
  NewData names used for cell-level comparison.

- ClassAudit:

  Tibble comparing variable classes for common non-key variables.

- ModifiedValues:

  Long-format tibble of cell-level value changes.

- VariableChangeSummary:

  Tibble summarizing changes by variable.

- TopChangedVariables:

  Top changed variables by number of modified values.

- SuspiciousChanges:

  Tibble of high-change-rate or class-change variables.

## Details

The function always returns both detailed cell-level changes and
summary-level outputs. Cell-level changes are stored in long format,
with one row per changed value.

Key variables are coerced to character internally before comparison
because IDs are often stored as numeric in one file and character in
another. The original key classes are preserved in the `KeyTypes`
output.

Variable names are also audited for common tibble/readxl name-repair
suffixes such as `...372`. These suffixes are ignored when identifying
variables added or removed, while raw variable names are still preserved
in the output.

If duplicate key combinations are detected, the function still runs, but
cell-level comparison is performed only for key combinations that are
unique in both datasets. Duplicate keys are returned separately and
flagged in `Checks`.

## Reading the result

Start with `Checks`, the traffic-light table. Everything that changed
between the two extracts is classified, counted, and explained there, so
a release can be signed off - or stopped - without reading any raw data.
`Summary` carries the same story as headline metrics, one row per
measure.

From there, `ModifiedValues` drills down to the individual cells that
moved, with the old and new value side by side and each change
classified; `VariableChangeSummary` rolls those up per variable, which
is usually what gets circulated to collaborators. `SummaryText` is a
plain-text version of the whole comparison, for a log or an email.

## Examples

``` r
# \donttest{
data(SampleData)

# A revised extract: corrected ages, lost values, dropped rows, a new column
df_OldData <- dplyr::mutate(
  SampleData,
  ParticipantID = seq_len(nrow(SampleData)),
  .before = 1
)

set.seed(4)
df_NewData <- df_OldData
df_NewData$age[1:5] <- df_NewData$age[1:5] + 1
df_NewData$MMP7[c(2, 9, 14)] <- NA
df_NewData$NewBiomarker <- stats::rnorm(nrow(df_NewData))
df_NewData <- df_NewData[-c(3, 7), ]

comparison <- CompareDatasets(
  OldData = df_OldData,
  NewData = df_NewData,
  keys = "ParticipantID"
)

ShowTable <- function(x, height = "300px") {
  htmltools::browsable(htmltools::HTML(as.character(
    FreezeTableHeader(x, height = height, full_width = TRUE)
  )))
}

# The traffic-light checks
ShowTable(comparison$Checks, height = "360px")


 Check 
```
