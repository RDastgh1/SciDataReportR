# Project standardized scores onto new data using external parameters

Project standardized scores onto new data using external parameters

## Usage

``` r
Project_ZScore(
  df,
  variables = NULL,
  parameters,
  ParameterInputType = c("df_parameter", "ZScoreObj", "ExternalDataframe"),
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- df:

  Data frame on which to project scores.

- variables:

  Character vector; if NULL, project onto all variables for which
  parameters exist and that are present in df.

- parameters:

  Source of parameters, interpreted by ParameterInputType:

  - "df_parameter": a data frame with cols Variable, N, Mean, SD

  - "ZScoreObj": output object from CalcZScore()

  - "ExternalDataframe": a raw reference data frame; parameters are
    estimated via CalcZScore() on that frame.

- ParameterInputType:

  One of "df_parameter", "ZScoreObj", "ExternalDataframe".

- names_prefix:

  Prefix for projected variable names.

- RetainLabels:

  Logical; if TRUE and Hmisc available, copy labels from df to new
  variables.

- RenameLabels:

  Logical; if TRUE, prefix labels the same way as names.

- center:

  Logical; used when ParameterInputType is "df_parameter" or
  "ExternalDataframe". Ignored for "ZScoreObj" (it uses stored flags).

- scale:

  Logical; same logic as `center`.

## Value

List with same structure as CalcZScore():

- ZScores: projected standardized variables only

- DataWithZ: df with projected scores appended

- Parameters: parameter data frame actually used for projection

- Center, Scale: flags used
