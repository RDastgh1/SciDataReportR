# Calculate M-scores using median and MAD

Calculate robust standardized scores using the median as the center and
median absolute deviation as the scale. M-scores are useful for skewed,
heavy-tailed, or outlier-prone variables where standard z-scores may be
overly influenced by the mean and standard deviation.

## Usage

``` r
CalcMScore(
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

  Data frame with variables to robust-standardize.

- variables:

  Character vector of variable names. If `NULL`, uses
  `SciDataReportR::getNumVars(df, Ordinal = FALSE)`.

- names_prefix:

  Prefix to prepend to variable names. Default is `"M_"`.

- RetainLabels:

  Logical; if `TRUE` and `Hmisc` is available, copy labels.

- RenameLabels:

  Logical; if `TRUE`, apply the same prefix to labels.

- center:

  Logical; if `TRUE`, subtract the median.

- scale:

  Logical; if `TRUE`, divide by the MAD.

- constant:

  Numeric scaling constant passed to
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html). Default is
  `1.4826`.

## Value

An object of class `"MScoreObj"`, a list with:

- `MScores`: data frame of M-score variables only

- `DataWithM`: original `df` plus M-score variables

- `Parameters`: data frame with `Variable`, `N`, `Median`, and `MAD`

- `Center`: logical flag used

- `Scale`: logical flag used

- `Constant`: MAD scaling constant used

## Details

An M-score is calculated as:

`(x - median(x)) / MAD(x)`

By default, [`stats::mad()`](https://rdrr.io/r/stats/mad.html) uses a
scaling constant of `1.4826`, which makes MAD approximately comparable
to the standard deviation when data are normally distributed.
