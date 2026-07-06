# Plot Interaction Effects Matrix

Creates a heatmap visualization of interaction effects between
continuous variables, showing whether interactions result in slope
reversals or maintain the same direction.

## Usage

``` r
PlotInteractionEffectsMatrix(
  data,
  interVar = NULL,
  outcome_vars = NULL,
  predictor_vars = NULL,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = FALSE,
  Data = lifecycle::deprecated(),
  outcomeVars = lifecycle::deprecated(),
  predictorVars = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome"),
  covars = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing the variables to be analyzed

- interVar:

  Character string specifying the interaction variable (moderator). Can
  be categorical or continuous.

- outcome_vars:

  Character vector of outcome variable names (displayed on rows)

- predictor_vars:

  Character vector of predictor variable names (displayed on columns)

- covariates:

  Character vector of covariate names to include in the models

- Relabel:

  Logical indicating whether to use variable labels if available
  (default: TRUE)

- Ordinal:

  Logical indicating whether to treat ordered factors as numeric
  (default: FALSE)

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- outcomeVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- predictorVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all interaction p-values at once
  (historical behavior). `"per_outcome"` corrects separately within each
  outcome variable (`outcome_vars`).

- covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A list containing:

- Unadjusted:

  List with unadjusted results including:

  - C: Matrix of interaction coefficients

  - S: Matrix of slope direction indicators (1 = same, -1 = reversed)

  - P: Matrix of p-values

  - D: Matrix of interaction coefficient signs

  - Slope1: Matrix of slopes at low values of interVar (mean - 1SD for
    continuous, reference group for categorical)

  - Slope2: Matrix of slopes at high values of interVar (mean + 1SD for
    continuous, comparison group for categorical)

  - plot: ggplot object of the heatmap

  - pvaltable: P-value table in wide format

- FDRCorrected:

  List with FDR-corrected results (same structure as Unadjusted)

- Relabel:

  Logical indicating whether relabeling was applied

- Covariates:

  Character vector of covariates used

- interVar:

  The interaction variable name

- raw_data:

  Processed data frame with all calculated values

## Details

The function fits linear models of the form: outcome ~ predictor \*
interVar + covariates for each combination of outcome and predictor
variables.

For continuous interaction variables, slopes are calculated at mean -
1SD and mean + 1SD. For categorical interaction variables, slopes are
calculated for each category.

The function then determines whether the interaction causes a slope
reversal (opposite signs) or maintains the same direction (same signs)
for the predictor-outcome relationship.

Color coding:

- Blue gradient: Significant interaction with slopes in the same
  direction

- Red gradient: Significant interaction with slope reversal

- White: Non-significant interaction (p \> 0.05)

- Grey: Missing data or model could not be fit

Darker colors indicate higher significance:

- \*\*\*: p \<= 0.001 (darkest)

- \*\*: p \<= 0.01 (medium)

- \*: p \<= 0.05 (light)

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

outcomes <- c("age", "ACE_CD143_Angiotensin_Converti",
              "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
              "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
              "Alpha_1_Microglobulin")
predictors <- c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
                "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E")

# With a categorical interaction variable (Diagnosis)
results <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "Diagnosis",
  outcome_vars = outcomes,
  predictor_vars = predictors
)
#> Warning: Ignoring unknown aesthetics: text

# With a continuous interaction variable (age)
results_cont <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "age",
  outcome_vars = setdiff(outcomes, "age"),
  predictor_vars = predictors
)
#> Warning: Ignoring unknown aesthetics: text

# Display the plot (colored tiles mark significant interactions)
results$Unadjusted$plot

# }
```
