# Plot Categorical Interaction Effects Matrix

This function calculates and visualizes the interaction effects between
categorical variables.

## Usage

``` r
PlotCatInteractionEffectsMatrix(
  Data,
  xVars,
  yVars = NULL,
  xVarLabels = NULL,
  yVarLabels = NULL,
  interVar
)
```

## Arguments

- Data:

  The dataset containing the variables.

- xVars:

  A character vector of the names of the x-axis categorical variables.

- yVars:

  A character vector of the names of the y-axis categorical variables.
  Defaults to NULL, in which case it takes the same values as xVars.

- xVarLabels:

  A character vector of labels for the x-axis variables. Defaults to
  NULL, in which case it takes the same values as xVars.

- yVarLabels:

  A character vector of labels for the y-axis variables. Defaults to
  NULL, in which case it takes the same values as yVars.

- interVar:

  The name of the interaction variable.

## Value

A list containing matrices of interaction coefficients, p-values, ggplot
objects for visualizations, and tables of FDR-corrected p-values.
