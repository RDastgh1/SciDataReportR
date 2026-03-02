# Plot Numerical Interaction Effects Matrix

This function calculates interaction effects between numerical variables
and plots them as matrices.

## Usage

``` r
PlotNumInteractionEffectsMatrix(
  Data,
  xVars,
  yVars = NULL,
  xVarLabels = NULL,
  yVarLabels = NULL,
  interVar = NULL,
  covars = NULL
)
```

## Arguments

- Data:

  The dataset containing the variables.

- xVars:

  A character vector of the names of the x-axis numerical variables.

- yVars:

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

- covars:

  A character vector of the names of covariate variables. Defaults to
  NULL.

## Value

A list containing matrices, ggplot objects for visualizations, and
tables of p-values.
