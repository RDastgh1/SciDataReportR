# Calculate Z-scores (or standardized scores) and return data + parameters

Calculate Z-scores (or standardized scores) and return data + parameters

## Usage

``` r
CalcZScore(
  df,
  variables = NULL,
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- df:

  Data frame with variables to standardize.

- variables:

  Character vector of variable names. If NULL, uses
  SciDataReportR::getNumVars(df).

- names_prefix:

  Prefix to prepend to variable names (default "Z\_").

- RetainLabels:

  Logical; if TRUE and Hmisc is available, copy labels.

- RenameLabels:

  Logical; if TRUE, apply the same prefix to labels.

- center:

  Logical; if TRUE, subtract the mean.

- scale:

  Logical; if TRUE, divide by the SD.

## Value

An object of class "ZScoreObj", a list with:

- ZScores: data frame of standardized variables only

- DataWithZ: original df + standardized variables

- Parameters: data frame with Variable, N, Mean, SD

- Center: logical flag used

- Scale: logical flag used
