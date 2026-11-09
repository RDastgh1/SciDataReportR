# Cox Proportional Hazards Regression Table

Fits one Cox proportional hazards model per event-predictor pair and
returns report-ready tables plus a tidy results dataframe compatible
with
[`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).

## Usage

``` r
MakeCoxRegressionTable(
  data,
  time_var,
  event_vars,
  predictor_vars,
  covariates = NULL,
  design = NULL,
  reference_levels = NULL,
  Standardize = FALSE,
  FDR = TRUE,
  FDRAlpha = 0.05,
  IncludeTrendSymbol = FALSE,
  CheckPH = TRUE,
  ReturnModels = FALSE,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  Data frame containing the time, event, predictor, and covariate
  variables. When `design` is supplied, `data` provides the
  variable-label contract and the design variables provide the modeling
  data.

- time_var:

  Character scalar naming the follow-up-time variable.

- event_vars:

  Character vector naming event indicators.

- predictor_vars:

  Character vector naming predictors.

- covariates:

  Optional character vector naming adjustment variables.

- design:

  Optional `survey.design` whose variables contain all model variables.
  Models are fit with
  [`survey::svycoxph()`](https://rdrr.io/pkg/survey/man/svycoxph.html)
  when supplied.

- reference_levels:

  Optional named list or named character vector giving the reference
  level for categorical predictors. Unspecified predictors retain their
  first factor level.

- Standardize:

  Logical. If `TRUE`, standardize numeric predictors and covariates
  within each model's complete-case analysis set.

- FDR:

  Logical. If `TRUE`, add `FDR`, adjusted across all returned rows using
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).

- FDRAlpha:

  Legacy significance-threshold alias, retained in metadata. An
  explicitly supplied `Pthresh` takes precedence; `PVar` selects the
  p-value source.

- IncludeTrendSymbol:

  Logical; add a plus sign for `0.05 < p <= 0.10`.

- CheckPH:

  Logical. If `TRUE`, calculate the global
  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
  p-value for ordinary Cox models. Survey models return `NA`.

- ReturnModels:

  Logical. If `TRUE`, retain fitted models in `ModelSummaries`.

- Relabel:

  Logical. If `TRUE`, use attached variable labels.

- TreatOrdinalAs:

  How ordered predictors and covariates are handled.

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

- `FormattedTable`: a report-facing `gt` table.

- `LargeTable`: a detailed `gt` table.

- `Results`: one row per predictor coefficient, with the tidy regression
  fields used by
  [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  plus `Events`, `Concordance`, `PH_PValue`, and (when requested) `FDR`.

- `ModelSummaries`: fitted models when `ReturnModels = TRUE`, otherwise
  `NULL`.

- `Metadata`: event coding, modeling engine, and analysis settings.

## Details

Each model has the form
`Surv(time_var, event) ~ predictor + covariates`. Without `design`,
models are fit with
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). When
`design` is a `survey.design`, the corresponding design variables are
used and models are fit with
[`survey::svycoxph()`](https://rdrr.io/pkg/survey/man/svycoxph.html).
`N` and `Events` are unweighted analytic row counts in both cases.

Numeric events must be coded 0/1. Logical events use `TRUE` as the
event. For two-level factors, the second factor level is the event.
Character event variables are deliberately rejected so the event
ordering is never implicit.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## See also

[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md),
[`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md),
[`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md)

## Examples

``` r
lung <- survival::lung
lung$event <- lung$status == 2
attr(lung$age, "label") <- "Age"

cox_results <- MakeCoxRegressionTable(
  data = lung,
  time_var = "time",
  event_vars = "event",
  predictor_vars = c("age", "sex")
)
cox_results$FormattedTable


  

```
