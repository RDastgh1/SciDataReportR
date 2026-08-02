# Univariate Regression Table

Creates a list of univariate regression tables with variable labels and
standardized coefficients (if specified).

`UnivariateRegressionTable()` was renamed to
`MakeUnivariateRegressionTable()` in SciDataReportR 20.5.0 to match the
package's `Make*` naming convention. It remains available as a
backwards-compatible synonym.

## Usage

``` r
MakeUnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  Data = lifecycle::deprecated(),
  OutcomeVars = lifecycle::deprecated(),
  PredictorVars = lifecycle::deprecated(),
  Covars = lifecycle::deprecated()
)

UnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE,
  Data = lifecycle::deprecated(),
  OutcomeVars = lifecycle::deprecated(),
  PredictorVars = lifecycle::deprecated(),
  Covars = lifecycle::deprecated()
)
```

## Arguments

- data:

  Dataframe containing the variables

- outcome_vars:

  Character vector of outcome variable names

- predictor_vars:

  Character vector of predictor variable names

- covariates:

  Character vector of covariate variable names (default: NULL)

- Standardize:

  Logical indicating whether to standardize numeric variables (default:
  FALSE)

- Method:

  Character. Regression method to use. `"auto"` detects linear
  regression for numeric outcomes and logistic regression for two-level
  outcomes. `"lm"` and `"logistic"` force one model family for all
  outcomes.

- LogisticExponentiate:

  Logical. If `TRUE`, logistic regression estimates are exponentiated
  and reported as odds ratios.

- ReturnModels:

  Logical. If `TRUE`, return fitted model objects in `ModelSummaries`.
  Default is `FALSE` to keep large screening runs lighter.

- Relabel:

  Logical; if TRUE (default), display attached variable labels.

- TreatOrdinalAs:

  How ordinal outcomes and predictors are handled.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- OutcomeVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- PredictorVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- Covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A list containing:

- FormattedTable: A `gt` table with formatted regression results

- LargeTable: A `gt` table with unformatted regression results

- Results: A tidy dataframe with one row per estimated term. Columns:
  `Outcome`, `OutcomeLabel`, `OutcomeFamily`, `EffectType`, `Predictor`,
  `PredictorLabel`, `Term`, `Level`, `TermLabel`, `N`, `Estimate`,
  `StdError`, `ConfLow`, `ConfHigh`, `PValue`, `Significant`, and
  `ReferenceValue`. This dataframe can be filtered and passed directly
  to
  [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).

- ModelSummaries: A list of fitted model objects when
  `ReturnModels = TRUE`, otherwise `NULL`

- Metadata: Outcome families and analysis settings
