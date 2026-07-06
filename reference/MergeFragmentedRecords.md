# Merge fragmented records into a single observation

Merges multiple rows representing fragments of the same observation
(e.g., participant visit, assessment session, study encounter, or
questionnaire administration) into a single row by selecting the first
non-missing value within each variable.

## Usage

``` r
MergeFragmentedRecords(
  data,
  id_var = "subject",
  date_var = "date",
  session_var = "session",
  keep_session = TRUE,
  session_name = "first_session",
  n_rows_name = "n_rows_collapsed",
  arrange_desc_session = FALSE,
  empty_strings_to_na = TRUE,
  df = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing fragmented records.

- id_var:

  Character string specifying the participant identifier variable.
  Default is `"subject"`.

- date_var:

  Character string specifying the visit or assessment date variable.
  Default is `"date"`.

- session_var:

  Character string specifying the session identifier variable used to
  order fragmented records. Default is `"session"`.

- keep_session:

  Logical. If `TRUE`, the first session value encountered within each
  group is retained. Default is `TRUE`.

- session_name:

  Character string specifying the name of the retained session variable.
  Default is `"first_session"`.

- n_rows_name:

  Character string specifying the name of the variable recording the
  number of rows merged. Default is `"n_rows_collapsed"`.

- arrange_desc_session:

  Logical. If `TRUE`, records are ordered by descending session number
  before merging. Default is `FALSE`.

- empty_strings_to_na:

  Logical. If `TRUE`, empty character strings are converted to missing
  values prior to merging. Default is `TRUE`.

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A data frame containing one row per unique combination of `id_var` and
`date_var`.

Additional variables may include:

- n_rows_collapsed:

  Number of fragmented rows merged.

- first_session:

  First session value retained, if `keep_session = TRUE`.

## Details

This function is useful when data collection is interrupted and
restarted, resulting in multiple partial records for the same
participant and date. Examples include tablet-based cognitive testing,
mobile applications, REDCap surveys, wearable device uploads, and
electronic assessments where internet connectivity or software issues
may split a single visit across multiple records.

Rows are grouped by the combination of `id_var` and `date_var`. Within
each group, observations are ordered by `session_var`, and the first
non-missing value encountered for each variable is retained.

For each variable within a participant-date grouping, the function
returns the first non-missing value after sorting by session order.

For example:


    subject   date         session   Stroop   Trails
    101       2024-01-01   1         50       NA
    101       2024-01-01   2         NA       100

    becomes

    subject   date         Stroop   Trails
    101       2024-01-01   50       100

No attempt is made to resolve conflicting non-missing values across
sessions. If multiple non-missing values exist for the same variable,
the first value encountered after sorting is retained.
