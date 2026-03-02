# PlotMiningMatrix

Generate a matrix of statistical relationships between outcome and
predictor variables. Returns unadjusted and FDR-adjusted results plus
corresponding dot-matrix plots.

## Usage

``` r
PlotMiningMatrix(
  Data,
  OutcomeVars,
  PredictorVars,
  Covariates = NULL,
  Relabel = TRUE,
  Parametric = TRUE
)
```

## Arguments

- Data:

  A data frame.

- OutcomeVars:

  Outcome variables (typically numeric, e.g., biomarker-by-batch
  columns).

- PredictorVars:

  Predictor variables (covariates to scan).

- Covariates:

  Optional adjustment covariates.

- Relabel:

  Whether to use variable labels instead of names.

- Parametric:

  If TRUE uses Pearson; if FALSE uses Spearman and non-parametric ANOVA
  where relevant.

## Value

List with Unadjusted (PvalTable, plot), FDRCorrected (PvalTable, plot),
method, Relabel, Covariates.
