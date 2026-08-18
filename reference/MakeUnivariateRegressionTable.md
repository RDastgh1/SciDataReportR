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

## Details

Fits one model per outcome-predictor pair and collects every result into
a single table. "Univariate" means each predictor is tested on its own:
for three outcomes and eight predictors this is 24 separate models, not
one model with eight terms. That is the screening step - finding which
associations exist at all - and it is deliberately different from
[`MultivariableRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MultivariableRegressionTable.md),
which puts all the predictors in one model and reports each one adjusted
for the others.

The model family is chosen per outcome. With `Method = "auto"` (the
default), numeric outcomes get linear regression and two-level outcomes
get logistic regression, so a mixed set of outcomes can be screened in
one call and each is still modeled correctly. `"lm"` or `"logistic"`
force one family for everything.

## Three views of the same results

The return value carries the same estimates in three shapes, for three
different jobs:

- `FormattedTable` - a wide `gt` table with one spanner per outcome,
  confidence intervals, stars, and significant results emphasized.

- `LargeTable` - the same wide layout with the individual model
  statistics, when the raw numbers matter more than the presentation.

- `Results` - one tidy row per estimated term, with `Estimate`,
  `StdError`, confidence bounds, `PValue`, and the labels. This is the
  one to work with programmatically: filter it, correct it with
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md),
  or pass it straight to
  [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).

Categorical predictors contribute one row per non-reference level, with
the level in `Level` and the level being compared against in
`ReferenceValue`.

## Options that change the estimates

`Standardize = TRUE` z-scores the numeric variables first, so
coefficients come back in standard deviations per standard deviation.
Use it when predictors are on scales that are not comparable to each
other - it is what makes "which predictor matters more" a meaningful
question.

`covariates` adds the named variables to every model, turning a screen
of raw associations into a screen of associations adjusted for them.
This is usually the difference between "these biomarkers track with the
outcome" and "these biomarkers track with the outcome beyond age and
sex".

`LogisticExponentiate = TRUE` (the default) reports logistic results as
odds ratios rather than log-odds, so the null value in the table is 1
rather than 0. `ReturnModels = TRUE` keeps the fitted model objects for
residual checks; it is off by default because a large screen would
otherwise hold hundreds of models in memory.

Running many models at once means many p-values. Nothing here corrects
for that automatically - the family is yours to define - so pass
`Results$PValue` through
[`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md)
once you have decided what the family is.

## See also

[`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md)
to visualize `Results`,
[`MultivariableRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MultivariableRegressionTable.md)
for mutually adjusted models, and
[`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md)
for multiple-comparison correction.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
# Make the reference example deterministic: age is visibly associated with
# AXL, so the displayed table includes a clear significant result.
ExampleData <- Labelled
set.seed(20260810)
ExampleData$AXL <- as.numeric(scale(ExampleData$age)) * 2 +
  stats::rnorm(nrow(ExampleData), sd = 0.25)
attr(ExampleData$AXL, "label") <- sjlabelled::get_label(Labelled$AXL)

vars_Outcomes <- c("AXL", "tau", "p_tau")
vars_Predictors <- c("age", "sex", "Adiponectin", "Cortisol")

# Three outcomes by four predictors: twelve separate models, one table
urt <- MakeUnivariateRegressionTable(
  data = ExampleData,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors
)

# Format 1: the report-ready table
urt$FormattedTable


  

```
