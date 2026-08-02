# Plot Numerical Interaction Effects Matrix

This function calculates interaction effects between numerical variables
and plots them as matrices.

## Usage

``` r
PlotNumInteractionEffectsMatrix(
  data,
  predictor_vars,
  outcome_vars = NULL,
  xVarLabels = NULL,
  yVarLabels = NULL,
  interVar = NULL,
  covariates = NULL,
  Data = lifecycle::deprecated(),
  xVars = lifecycle::deprecated(),
  yVars = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  covars = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset containing the variables.

- predictor_vars:

  A character vector of the names of the x-axis numerical variables.

- outcome_vars:

  A character vector of the names of the y-axis numerical variables.
  Defaults to NULL.

- xVarLabels:

  A character vector of labels for the x-axis variables. Defaults to
  NULL.

- yVarLabels:

  A character vector of labels for the y-axis variables. Defaults to
  NULL.

- interVar:

  The interaction variable. Defaults to NULL.

- covariates:

  A character vector of the names of covariate variables. Defaults to
  NULL.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- xVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- yVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all interaction p-values at once
  (historical behavior). `"per_outcome"` corrects separately within each
  y-axis variable (`outcome_vars`).

- covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A list containing matrices, ggplot objects for visualizations, and
tables of p-values.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotNumInteractionEffectsMatrix(
  Labelled,
  predictor_vars = c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                     "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                     "Apolipoprotein_A2", "Apolipoprotein_B",
                     "Apolipoprotein_CI", "Apolipoprotein_CIII",
                     "Apolipoprotein_D", "Apolipoprotein_E"),
  outcome_vars = c("ACE_CD143_Angiotensin_Converti",
                   "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
                   "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
                   "Alpha_1_Microglobulin"),
  interVar = "age"
)
#> Joining with `by = join_by(X, Y)`
#> Joining with `by = join_by(X, Y)`

# Raw p-value interaction matrix
result$p


# FDR-adjusted interaction matrix
result$p_FDR
```
