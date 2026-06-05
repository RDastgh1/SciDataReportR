# Validate a merged dataset against its source datasets

Audit a merged dataset against the two source datasets used to create
it. This function is designed to catch common merge problems before
analysis, including incompatible key types, duplicate key combinations,
unexpected overlapping variables, unresolved `.x` / `.y` variables, and
value conflicts between duplicated variables.

## Usage

``` r
ValidateMerge(LeftData, RightData, MergedData, Keys)
```

## Arguments

- LeftData:

  A data frame used as one source for the merge.

- RightData:

  A data frame used as the other source for the merge.

- MergedData:

  The merged data frame to audit.

- Keys:

  Character vector of key variables intended to define the merge.
  Multiple keys are supported, such as `c("study_id", "TimePoint")`.

## Value

A list with merge validation results, including:

- SummaryText:

  A plain-text summary of the merge audit.

- ReadyForAnalysis:

  Logical value indicating whether major merge-integrity issues were
  absent.

- Summary:

  One-row tibble with core merge metrics.

- Fingerprint:

  Tibble comparing rows, columns, and unique key combinations across
  datasets.

- KeyTypes:

  Tibble showing key variable classes before coercion.

- Checks:

  Tibble summarizing validation checks and pass/warning/fail status.

- SuggestedActions:

  Tibble with suggested review steps.

- Relationship:

  Tibble describing key uniqueness in LeftData, RightData, and
  MergedData.

- IDCoverage:

  List containing Matching, LeftOnly, and RightOnly key combinations.

- DuplicateKeys:

  List containing duplicated key rows from Left, Right, and Merged
  datasets.

- OverlappingVariables:

  Tibble of variables present in both source datasets but not listed as
  keys.

- PotentialMergeRisk:

  Tibble of overlap variables that could have affected an unspecified
  dplyr join.

- JoinAudit:

  Tibble showing variables present in both source datasets and whether
  they were specified as merge keys.

- OverlapAudit:

  Tibble auditing all source variables by source presence and key
  status.

- UnresolvedDuplicateVariables:

  Tibble of `.x` / `.y` duplicate variable pairs still present in
  MergedData.

- DuplicateVariables:

  Tibble with agreement, conflict counts, and classes for duplicated
  variables.

- SuspiciousConflicts:

  Subset of duplicated variables with low agreement or mismatched
  classes.

- VariableConflicts:

  Long-format tibble of record-level value conflicts.

## Details

This function does not perform a merge. It assumes the user has already
created a merged dataset and wants to check whether the merge can be
trusted.
