# Create Summary Table

Generate a descriptive summary table for specified variables in a
dataset.

## Usage

``` r
CreateSummaryTable(
  data,
  variables = NULL,
  digits = 2,
  Relabel = TRUE,
  Ordinal = FALSE,
  ScrollBoxHeight = "700px",
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  numdecimals = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset containing the variables of interest.

- variables:

  A character vector specifying the variables for which summary
  statistics will be calculated.

- digits:

  Number of decimal places to round the summary statistics.

- Relabel:

  Logical, indicating whether to use variable labels as column headers.

- Ordinal:

  Logical, indicating whether ordinal variables should be included in
  the summary.

- ScrollBoxHeight:

  Height of the scroll box for displaying the table.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- numdecimals:

  **Deprecated** (since 19.15.0). Use `digits` instead.

## Value

A formatted HTML table displaying summary statistics.
