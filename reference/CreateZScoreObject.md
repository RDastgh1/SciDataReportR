# Calculate Z-scores (or standardized scores) and return data + parameters

Standardizes each variable to a common scale and, critically, returns
the constants used to do it so the identical transformation can be
replayed on other data later.

`CalcZScore()` has been superseded by `CreateZScoreObject()`. It remains
available as a backwards-compatible alias and returns the same reusable
Z-score object.

## Usage

``` r
CreateZScoreObject(
  data,
  variables = NULL,
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE,
  df = lifecycle::deprecated()
)

CalcZScore(...)
```

## Arguments

- data:

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

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

- ...:

  Arguments passed to `CreateZScoreObject()`.

## Value

An object of class "ZScoreObj", a list with:

- ZScores: data frame of standardized variables only

- DataWithZ: original df + standardized variables

- Parameters: data frame with Variable, N, Mean, SD

- Center: logical flag used

- Scale: logical flag used

## The equation

For each variable, every value is expressed as its distance from that
variable's mean, measured in standard deviations:

\$\$z_i = \frac{x_i - \bar{x}}{s}\$\$

where \\\bar{x}\\ and \\s\\ are the mean and standard deviation of the
variable, both computed with `na.rm = TRUE`. The result has mean 0 and
standard deviation 1, so `z = -2` means the same thing - two standard
deviations below average - no matter which variable it came from.

`center = FALSE` drops the \\-\bar{x}\\ term and `scale = FALSE` drops
the division by \\s\\, which is how the "Center Only" and "Scale Only"
options elsewhere in the package are expressed.

## Why the parameters are returned

The \\\bar{x}\\ and \\s\\ used for each variable are stored in
`Parameters`, and this is the whole point of returning an object rather
than just a standardized data frame.

A z-score is only interpretable relative to the sample it was computed
from. If a validation cohort, a follow-up visit, or a new site is
standardized against *its own* mean and SD, then every group is centered
at 0 by construction and any real difference between them is scaled
away - the clinical cohort looks exactly like the reference cohort.
Worse, the resulting columns are not comparable to the original ones
even though they carry the same names.

Passing this object to
[`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md)
applies the *frozen* training constants to the new data instead, so new
observations land on the original scale and group differences survive.
The same principle is why the clustering pipelines store a
`ZScoreObject` alongside the fitted model: a projected participant must
pass through exactly the scaling the model was trained on.

## What the example demonstrates

The first pair of density plots shows why standardizing helps at all:
before it, each biomarker sits at its own location, so no single axis is
meaningful for all of them and they cannot be compared or pooled; after
it, every variable is centered at 0 with a standard deviation of 1, so
they share one axis and a value of -2 means the same thing everywhere.

The second part splits the cohort in half and scores the same
participants two ways. Freezing the training mean and SD and applying
them through
[`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md)
keeps whatever offset the new sample really has. Re-standardizing the
new half against itself forces every variable back to mean 0, so the new
cohort can never differ from the training cohort no matter what the
measured values were. That difference is the information loss the stored
parameters exist to prevent.

## See also

[`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md)
to apply stored parameters to new data, and
[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md)
when the reference values should also be adjusted for covariates such as
age or education.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
vars_Biomarkers <- c("AXL", "Adiponectin", "Ferritin", "MMP7", "tau", "p_tau")

z_obj <- CreateZScoreObject(df_Labelled, variables = vars_Biomarkers)

# The stored centering and scaling constants
htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(z_obj$Parameters, full_width = TRUE)
)))


   
```
