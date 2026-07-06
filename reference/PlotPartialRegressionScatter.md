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
  Covariates = lifecycle::deprecated()
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

## Examples

``` r
data(SampleData)

result <- PlotPartialRegressionScatter(
  SampleData,
  IndepVar = "age",
  DepVar = "AXL",
  covariates = "Adiponectin"
)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the SciDataReportR package.
#>   Please report the issue at
#>   <https://github.com/RDastgh1/SciDataReportR/issues>.

result$plot
#> `geom_smooth()` using formula = 'y ~ x'
```
