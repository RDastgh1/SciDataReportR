# Compare two versions of a dataset

Compare an old dataset and a new dataset using one or more key
variables. This function identifies record-level, variable-level, and
cell-level changes between dataset versions. It is useful when reviewing
updated data extracts, revised REDCap exports, cleaned spreadsheet
versions, or vendor-delivered dataset updates.

## Usage

``` r
CompareDatasets(OldData, NewData, Keys)
```

## Arguments

- OldData:

  A data frame representing the earlier dataset version.

- NewData:

  A data frame representing the newer dataset version.

- Keys:

  Character vector of key variables used to align records across the two
  datasets. Multiple keys are supported, such as
  `c("study_id", "TimePoint")`.

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
