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
  data,
  predictor_vars,
  outcome_vars,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  min_n = 4,
  Data = lifecycle::deprecated(),
  xVars = lifecycle::deprecated(),
  yVars = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome"),
  covars = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data.frame containing the dataset.

- predictor_vars:

  Character vector of x-axis categorical variables.

- outcome_vars:

  Character vector of y-axis categorical variables. If NULL, uses xVars.

- covariates:

  Optional character vector of covariate variables used for
  stratification (not adjustment).

- Relabel:

  Logical; whether to use variable labels (sjlabelled) in the plot.

- Ordinal:

  Logical; included for backward compatibility (currently unused here).

- min_n:

  Minimum number of complete observations required for a tested
  association.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- xVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- yVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all p-values at once (historical behavior).
  `"per_outcome"` corrects separately within each outcome: outcomes are
  the y-axis variables (`outcome_vars`).

- covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

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

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotChiSqCovar(
  Labelled,
  predictor_vars = c("Diagnosis", "Genotype"),
  outcome_vars = c("Diagnosis", "Genotype")
)
#> Warning: There were 3 warnings in `dplyr::summarise()`.
#> The first warning was:
#> ℹ In argument: `test = list(tryCatch(stats::chisq.test(XVal, YVal), error =
#>   function(e) NULL))`.
#> ℹ In group 2: `XVar = "Diagnosis"`, `YVar = "Genotype"`.
#> Caused by warning in `stats::chisq.test()`:
#> ! Chi-squared approximation may be incorrect
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 2 remaining warnings.

result$p
```
