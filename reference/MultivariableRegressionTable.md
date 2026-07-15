# Multivariable regression table

Fit one multivariable regression model per outcome and return a stable,
label-aware regression object for downstream tables, diagnostics, and
plots.

## Usage

``` r
MultivariableRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = TRUE,
  Relabel = TRUE,
  FDR = TRUE,
  FDRAlpha = 0.05,
  Method = c("lm", "ridge", "lasso", "elasticnet"),
  CVFolds = 10,
  Lambda = c("lambda.min", "lambda.1se"),
  Seed = 123,
  MissingDataStrategy = c("drop_sparse_impute", "impute", "complete_cases",
    "drop_sparse_complete_cases"),
  MaxMissingPredictor = 0.3,
  ImputeMethod = c("median_mode"),
  MinCompleteCases = NULL,
  Data = lifecycle::deprecated(),
  OutcomeVars = lifecycle::deprecated(),
  PredictorVars = lifecycle::deprecated(),
  Covars = lifecycle::deprecated()
)
```

## Arguments

- data:

  Data frame containing outcomes, predictors, and covariates.

- outcome_vars:

  Character vector of outcome variable names.

- predictor_vars:

  Character vector of predictor variable names.

- covariates:

  Optional character vector of covariate variable names. Covariates are
  treated as mandatory adjustments: for penalized methods (`"ridge"`,
  `"lasso"`, `"elasticnet"`) they are exempted from the penalty
  (`penalty.factor = 0`), so they are never shrunk or selected out of
  the model.

- Standardize:

  Logical. If `TRUE`, ordinary models are fit on standardized continuous
  variables for the primary estimate. Standardized coefficients are
  always calculated separately regardless of this setting.

- Relabel:

  Logical. If `TRUE`, use variable labels from `sjlabelled` when
  available.

- FDR:

  Logical. If `TRUE`, calculate FDR-adjusted p-values for ordinary
  regression terms.

- FDRAlpha:

  Numeric FDR threshold retained in metadata.

- Method:

  Regression method. One of `"lm"`, `"ridge"`, `"lasso"`, or
  `"elasticnet"`.

- CVFolds:

  Number of cross-validation folds for penalized models.

- Lambda:

  Lambda selection rule for penalized models. One of `"lambda.min"` or
  `"lambda.1se"`.

- Seed:

  Random seed used for deterministic cross-validation folds.

- MissingDataStrategy:

  Missing-data handling strategy. The default, `"drop_sparse_impute"`,
  drops sparse predictors and covariates, then imputes remaining
  predictor missingness.

- MaxMissingPredictor:

  Maximum allowed missingness proportion for predictors and covariates
  before they are dropped by sparse-drop strategies. Default is `0.30`.

- ImputeMethod:

  Imputation method for predictor/covariate missingness. Currently
  `"median_mode"`: median for numeric variables and mode for factor,
  character, and logical variables.

- MinCompleteCases:

  Optional minimum number of modeling rows required after missing-data
  handling.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- OutcomeVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- PredictorVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- Covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A named list with stable components: `Models`, `FormattedTable`,
`LargeTable`, `RegressionMatrix`, `VariableImportanceMatrix`,
`Predictions`, `Diagnostics`, `ModelSummary`, `Multicollinearity`,
`Plots`, and `Metadata`. `FormattedTable` is a report-facing `gt` table
grouped by outcome, matching the style of
[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md):
predictor rows only, a combined `Estimate (95% CI)` cell, and bold
significant p-values. `LargeTable` is a data frame holding the full
per-term detail (including covariate rows) for programmatic use, plus an
`Aliased` flag marking perfectly collinear terms the model dropped.
`ModelSummary` reports per-outcome `Converged`, `SeparationDetected`,
and `AliasedTermCount`. For ordinary (`"lm"`) logistic fits,
quasi-complete separation is detected (fitted probabilities pinned at
0/1, exploded standardized coefficients, or non-convergence); the
affected model's estimates are blanked (`NA`) and `Converged` is set to
`FALSE` so unreliable coefficients do not propagate into tables or
plots. `ModelSummary` also carries an omnibus model test per outcome
(`ModelStat`, `ModelStatType`, `ModelPValue`): an F-test for linear
models and a likelihood-ratio test for logistic models (`NA` for
penalized fits, which have no valid classical omnibus test). `Plots`
contains ggplot objects built from the stored result tables and
predictions without refitting models; the coefficient heatmap uses
robust, clamped fill limits so a single extreme value cannot dominate
the scale, and each outcome column is annotated at the top with its
omnibus p-value (ordinary models) or cross-validated deviance explained
(penalized models) to discourage interpreting coefficients from a model
that is not significant overall.

## Examples

``` r
# \donttest{
data(SampleData)

result <- MultivariableRegressionTable(
  SampleData,
  outcome_vars = "AXL",
  predictor_vars = c("Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin"),
  covariates = "age"
)

# Display the regression coefficient matrix plot
result$Plots$RegressionMatrix

# }
```
