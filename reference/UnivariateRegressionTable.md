# Univariate Regression Table

Creates a list of univariate regression tables with variable labels and
standardized coefficients (if specified).

## Usage

``` r
UnivariateRegressionTable(
  Data,
  OutcomeVars,
  PredictorVars,
  Covars = NULL,
  Standardize = FALSE
)
```

## Arguments

- Data:

  Dataframe containing the variables

- OutcomeVars:

  Character vector of outcome variable names

- PredictorVars:

  Character vector of predictor variable names

- Covars:

  Character vector of covariate variable names (default: NULL)

- Standardize:

  Logical indicating whether to standardize numeric variables (default:
  FALSE)

## Value

A list containing:

- FormattedTable: A merged table with formatted regression results

- LargeTable: A merged table with unformatted regression results

- ModelSummaries: A list of lm model summaries
