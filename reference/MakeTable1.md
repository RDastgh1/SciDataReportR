# Create Summary Table using gtsummary

This function is a wrapper around
[`gtsummary::tbl_summary`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
that ensures continuous variables are treated as continuous.

## Usage

``` r
MakeTable1(
  DataFrame,
  Variables = NULL,
  TreatOrdinalAs = "Continuous",
  AutoDetectDistribution = FALSE,
  IncludeMissing = "ifany"
)
```

## Arguments

- DataFrame:

  The dataframe to create the summary table from.

- Variables:

  Optional. A character vector specifying the names of variables to
  include in the summary table. If NULL, all variables are included.

- TreatOrdinalAs:

  Character. Specifies how ordinal variables should be treated. Can be
  "Continuous", "Categorical", or "Both".

- AutoDetectDistribution:

  Logical. If TRUE, the function will attempt to automatically detect
  the distribution of variables. Default is FALSE.

- IncludeMissing:

  Character matching gtsummary criteria. Can be "no", "ifany", or
  "always". Default is "ifany"

## Value

A summary table created using gtsummary.
