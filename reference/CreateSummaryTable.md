# Create Summary Table

Generate a descriptive summary table for specified variables in a
dataset.

## Usage

``` r
CreateSummaryTable(data, Ordinal = FALSE, ScrollBoxHeight = 400)
```

## Arguments

- Data:

  The dataset containing the variables of interest.

- Variables:

  A character vector specifying the variables for which summary
  statistics will be calculated.

- numdecimals:

  Number of decimal places to round the summary statistics.

- Relabel:

  Logical, indicating whether to use variable labels as column headers.

- Ordinal:

  Logical, indicating whether ordinal variables should be included in
  the summary.

- ScrollBoxHeight:

  Height of the scroll box for displaying the table.

## Value

A formatted HTML table displaying summary statistics.
