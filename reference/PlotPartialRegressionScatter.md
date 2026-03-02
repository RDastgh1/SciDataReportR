# Partial Regression Plot

Generate a partial regression plot for a specified independent and
dependent variable while adjusting for covariates. In addition to the
figure, key parameters such as the correlation method, whether
relabeling was used, the covariates, R², p-value, and sample size are
returned.

## Usage

``` r
PlotPartialRegressionScatter(
  DataFrame,
  IndepVar,
  DepVar,
  Covariates = NULL,
  Relabel = TRUE
)
```

## Arguments

- DataFrame:

  The dataset to use.

- IndepVar:

  A string specifying the independent variable.

- DepVar:

  A string specifying the dependent variable.

- Covariates:

  A character vector of covariate names for adjustment. Defaults to
  NULL.

- Relabel:

  Logical indicating whether to use labelled names from the data (using
  sjlabelled::get_label). Defaults to TRUE.

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

  The R² of the partial regression model.

- p_value:

  The p-value for the independent variable coefficient.

- n:

  The sample size (number of complete cases).

- equation:

  The regression equation string.
