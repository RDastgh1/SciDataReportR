# Plot Chi-Square Tests for Categorical Associations (optionally stratified by covariates)

Conducts Chi-square tests between sets of categorical variables and
visualizes the results. NOTE: Chi-square tests do not natively "adjust"
for covariates. If `covars` are provided, this function can (optionally)
run tests *within strata* (each combination of covariate levels), and
combine p-values across strata (Fisher's method) for a single summary
p-value per pair. If you need true covariate adjustment, use
regression-based models (logistic/multinomial).

## Usage

``` r
PlotChiSqCovar(
  Data,
  xVars,
  yVars,
  covars = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  min_n = 4
)
```

## Arguments

- Data:

  A data.frame containing the dataset.

- xVars:

  Character vector of x-axis categorical variables.

- yVars:

  Character vector of y-axis categorical variables. If NULL, uses xVars.

- covars:

  Optional character vector of covariate variables used for
  stratification (not adjustment).

- Relabel:

  Logical; whether to use variable labels (sjlabelled) in the plot.

- Ordinal:

  Logical; included for backward compatibility (currently unused here).

- Stratify:

  Logical; if TRUE and covars provided, run chi-square within covariate
  strata and combine p-values with Fisher's method. If FALSE, covars are
  ignored.

- MinExpected:

  Minimum expected cell count threshold for chi-square validity warning.
  (used to annotate; test still runs).

## Value

A list with:

- p:

  ggplot for unadjusted p-values

- pvaltable:

  wide table of unadjusted p-values

- p_FDR:

  ggplot for FDR-adjusted p-values

- pvaltable_FDR:

  wide table of FDR-adjusted p-values

- details:

  long table with diagnostics (n, warnings, strata info)
