# Plot Interaction Effects Matrix

Creates a heatmap visualization of interaction effects between
continuous variables, showing whether interactions result in slope
reversals or maintain the same direction.

## Usage

``` r
PlotInteractionEffectsMatrix(
  Data,
  interVar = NULL,
  outcomeVars = NULL,
  predictorVars = NULL,
  covars = NULL,
  Relabel = TRUE,
  Ordinal = FALSE
)
```

## Arguments

- Data:

  A data frame containing the variables to be analyzed

- interVar:

  Character string specifying the interaction variable (moderator). Can
  be categorical or continuous.

- outcomeVars:

  Character vector of outcome variable names (displayed on rows)

- predictorVars:

  Character vector of predictor variable names (displayed on columns)

- covars:

  Character vector of covariate names to include in the models

- Relabel:

  Logical indicating whether to use variable labels if available
  (default: TRUE)

- Ordinal:

  Logical indicating whether to treat ordered factors as numeric
  (default: FALSE)

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

- \*\*\*: p ≤ 0.001 (darkest)

- \*\*: p ≤ 0.01 (medium)

- \*: p ≤ 0.05 (light)

## Examples

``` r
if (FALSE) { # \dontrun{
# With categorical interaction variable
results <- PlotInteractionEffectsMatrix(
  Data = mydata,
  interVar = "HIV_status",
  outcomeVars = sleep_vars,
  predictorVars = metabolites,
  covars = c("age", "sex")
)

# With continuous interaction variable
results <- PlotInteractionEffectsMatrix(
  Data = mydata,
  interVar = "age",
  outcomeVars = sleep_vars,
  predictorVars = metabolites,
  covars = c("sex", "BMI")
)

# Display the plot
print(results$Unadjusted$plot)
} # }
```
