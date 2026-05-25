# Calculate robust M-scores for numeric variables

Calculate median/MAD-based M-scores for selected numeric variables and
return both the transformed data and the parameters needed to review or
reuse the transformation.

## Usage

``` r
CreateMScoreObject(
  df,
  variables = NULL,
  names_prefix = "M_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE,
  constant = 1.4826
)
```

## Arguments

- df:

  A data frame.

- variables:

  Character vector of numeric variables to transform. If `NULL`, numeric
  variables are detected with
  [`getNumVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getNumVars.md).

- names_prefix:

  Prefix for generated M-score columns.

- RetainLabels:

  Logical; keep variable labels when possible.

- RenameLabels:

  Logical; rename generated labels when labels are retained.

- center:

  Logical; subtract the median before scaling.

- scale:

  Logical; divide by the median absolute deviation.

- constant:

  Scaling constant passed to
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html).

## Value

An object of class `"MScoreObj"`, a list with:

- `MScores`: data frame of M-score variables only

- `DataWithM`: original `df` plus M-score variables

- `Parameters`: data frame with `Variable`, `N`, `Median`, and `MAD`

- `Center`: logical flag used

- `Scale`: logical flag used

- `Constant`: MAD scaling constant used
