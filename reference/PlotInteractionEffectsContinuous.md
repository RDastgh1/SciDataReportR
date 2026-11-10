# Plot Single Interaction Effect

Creates a scatter plot with regression lines showing the interaction
between a predictor and outcome variable, moderated by either a
continuous or categorical variable.

## Usage

``` r
PlotInteractionEffectsContinuous(
  data,
  interVar = NULL,
  outcome_var = NULL,
  predictor_var = NULL,
  covariates = NULL,
  n_lines = 3,
  alpha = 0.6,
  point_size = 2,
  Data = lifecycle::deprecated(),
  outcomeVar = lifecycle::deprecated(),
  predictorVar = lifecycle::deprecated(),
  covars = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue",
  adjust = TRUE,
  standardize = FALSE,
  covariate_interactions = FALSE,
  TreatOrdinalAs = "Exclude",
  interaction_results = NULL,
  source_formula = NULL
)
```

## Arguments

- data:

  A data frame containing the variables to be analyzed

- interVar:

  Character string specifying the interaction variable (moderator)

- outcome_var:

  Character string specifying the outcome variable

- predictor_var:

  Character string specifying the predictor variable

- covariates:

  Character vector of covariate names to include in the model

- n_lines:

  For continuous moderators, number of lines to plot (default: 3 for
  low/med/high)

- alpha:

  Transparency level for points (default: 0.6)

- point_size:

  Size of points (default: 2)

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- outcomeVar:

  **Deprecated** (since 19.15.0). Use `outcome_var` instead.

- predictorVar:

  **Deprecated** (since 19.15.0). Use `predictor_var` instead.

- covars:

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

- adjust:

  Display model-based covariate-adjusted partial residuals.

- standardize:

  Use pooled model-sample SD units for both axes.

- covariate_interactions:

  Allow moderator-specific covariate effects.

- TreatOrdinalAs:

  Ordinal treatment passed to PlotInteractionEffectsMatrix.

- interaction_results:

  Optional source interaction row, retaining its original
  multiple-testing q-value in the annotation.

- source_formula:

  Optional stored matrix formula; validated against the supplied model
  specification and used for adjusted displays.

## Value

A ggplot object showing the interaction effect

## Details

For categorical moderators, separate regression lines are plotted for
each category. For continuous moderators, regression lines are plotted
at the mean and +/- 1 SD.

Model-based partial residuals preserve the fitted interaction slopes.
Ribbons show 95% confidence intervals for the displayed model component.

Variable labels are used if available in the data frame.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axis titles and legend
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Categorical moderator: one regression line per Diagnosis group
PlotInteractionEffectsContinuous(
  Labelled,
  interVar = "Diagnosis",
  outcome_var = "AXL",
  predictor_var = "age"
)


# Continuous moderator: lines at mean and +/- 1 SD of the moderator
PlotInteractionEffectsContinuous(
  Labelled,
  interVar = "Adiponectin",
  outcome_var = "AXL",
  predictor_var = "age"
)
```
