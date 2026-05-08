# PlotMiningMatrix

Generate a matrix of statistical relationships between variables.

## Usage

``` r
PlotMiningMatrix(
  Data,
  OutcomeVars,
  PredictorVars = NULL,
  Covariates = NULL,
  Relabel = TRUE,
  Parametric = TRUE
)
```

## Arguments

- Data:

  A data frame.

- OutcomeVars:

  Outcome variables.

- PredictorVars:

  Predictor variables. If NULL, uses OutcomeVars.

- Covariates:

  Optional covariates (reserved for future use).

- Relabel:

  Use labels instead of names.

- Parametric:

  Use parametric tests.

## Value

List with tables and plots.
