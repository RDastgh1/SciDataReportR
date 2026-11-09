# Partial Regression Plot

Generate a partial regression plot for a specified independent and
dependent variable while adjusting for covariates. In addition to the
figure, key parameters such as the correlation method, whether
relabeling was used, the covariates, R-squared, p-value, and sample size
are returned.

## Usage

``` r
PlotPartialRegressionScatter(
  data,
  IndepVar,
  DepVar,
  covariates = NULL,
  Relabel = TRUE,
  DataFrame = lifecycle::deprecated(),
  Covariates = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  The dataset to use.

- IndepVar:

  A string specifying the independent variable.

- DepVar:

  A string specifying the dependent variable.

- covariates:

  A character vector of covariate names for adjustment. Defaults to
  NULL.

- Relabel:

  Logical indicating whether to use labelled names from the data (using
  sjlabelled::get_label). Defaults to TRUE.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

- Pthresh:

  Significance threshold, default `0.05`. An explicitly supplied value
  takes precedence over a legacy significance threshold argument.

- PVar:

  P-value column used for significance decisions. Defaults to raw
  `"PValue"`; select `"PValueFDR"` to use Benjamini-Hochberg adjusted
  values. Explicit legacy p-column names remain accepted. Existing,
  explicitly named FDR plots and adjusted-value columns are retained.

  Explicit `Pthresh` takes precedence over `Alpha`, `FDRAlpha`,
  `label_q`, or `adjusted_significance_threshold`. Pairwise
  `star_p = "adjusted"` retains its legacy adjustment method unless
  `PVar` is supplied; `PVar = "PValueFDR"` always selects BH.
  `star_p = "none"` suppresses legacy pairwise stars. Stars retain
  conventional p-value tiers within the requested cutoff.

## Value

A list containing:

- plot:

  A ggplot2 object representing the partial regression plot.

- method:

  The correlation method (as provided).

- Relabel:

  Logical; whether relabeling was applied.

- Covariates:

  The vector of covariates.

- r2:

  The R-squared of the partial regression model.

- p_value:

  The p-value for the independent variable coefficient.

- n:

  The sample size (number of complete cases).

- equation:

  The regression equation string.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotPartialRegressionScatter(
  Labelled,
  IndepVar = "age",
  DepVar = "AXL",
  covariates = "Adiponectin"
)

result$plot
#> `geom_smooth()` using formula = 'y ~ x'
```
