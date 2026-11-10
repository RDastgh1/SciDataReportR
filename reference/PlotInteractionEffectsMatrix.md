# Plot Interaction Effects Matrix

Creates a heatmap visualization of interaction effects between
continuous variables, showing whether interactions result in slope
reversals or maintain the same direction.

## Usage

``` r
PlotInteractionEffectsMatrix(
  data,
  interVar = NULL,
  outcome_vars = NULL,
  predictor_vars = NULL,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = lifecycle::deprecated(),
  TreatOrdinalAs = "Exclude",
  Data = lifecycle::deprecated(),
  outcomeVars = lifecycle::deprecated(),
  predictorVars = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  covars = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue",
  covariate_interactions = FALSE
)
```

## Arguments

- data:

  A data frame containing the variables to be analyzed

- interVar:

  Character string specifying the interaction variable (moderator). Can
  be categorical or continuous.

- outcome_vars:

  Character vector of outcome variable names (displayed on rows)

- predictor_vars:

  Character vector of predictor variable names (displayed on columns)

- covariates:

  Character vector of covariate names to include in the models

- Relabel:

  Logical indicating whether to use variable labels if available
  (default: TRUE)

- Ordinal:

  **Deprecated** (since 20.20.0). Use `TreatOrdinalAs` instead.

- TreatOrdinalAs:

  How ordinal variables are handled. This interaction matrix accepts
  `"Exclude"`, `"Continuous"`, or `"Categorical"`. Categorical ordinal
  outcomes are not supported by the linear models.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- outcomeVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- predictorVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- fdr_scope:

  One of `"matrix"` (default), `"per_outcome"`, or `"per_predictor"`,
  passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` adjusts all interaction p-values together. `"per_outcome"`
  adjusts across predictor columns separately within each outcome row.
  `"per_predictor"` adjusts across outcome rows separately within each
  predictor column.

- covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

- Pthresh:

  Interaction significance cutoff, default 0.05. Cells are displayed
  when their raw (unadjusted plot) or BH-adjusted (FDR plot) interaction
  p-value is less than or equal to this cutoff.

- PVar:

  P-value column used for `Results$Significant`, default `"PValue"`.
  `"PValueFDR"` uses BH-adjusted values. The named raw and FDR matrix
  plots always use their respective p-value columns.

- covariate_interactions:

  Logical, default FALSE. Interact every covariate with the moderator,
  allowing moderator-specific covariate effects.

## Value

A list containing:

- Unadjusted:

  List with unadjusted results including:

  - C: Matrix of interaction coefficients

  - S: Matrix of slope direction indicators (1 = same, -1 = reversed)

  - P: Matrix of p-values

  - D: Matrix of interaction coefficient signs

  - Slope1: Matrix of slopes at low values of interVar (mean - 1SD for
    continuous, reference group for categorical)

  - Slope2: Matrix of slopes at high values of interVar (mean + 1SD for
    continuous, comparison group for categorical)

  - Slope1Std, Slope2Std, SlopeDeltaStd: Common-scale standardized slope
    matrices

  - plot: ggplot object of the heatmap

  - pvaltable: P-value table in wide format

- FDRCorrected:

  List with FDR-corrected results (same structure as Unadjusted)

- Relabel:

  Logical indicating whether relabeling was applied

- Covariates:

  Character vector of covariates used

- interVar:

  The interaction variable name

- raw_data:

  Processed data frame with all calculated values

- Results:

  Long-form interaction and simple-slope inference

- Metadata:

  Model provenance, failures, scaling, and analysis choices

- SlopeDifferenceLimits:

  Shared symmetric plot limits

## Details

Models use `outcome ~ predictor * interVar + covariates` by default.
With `covariate_interactions = TRUE`, every covariate also interacts
with the moderator. For categorical moderators, point estimates of group
slopes then equal separate within-group linear fits on the same complete
cases. Their standard errors use the pooled model residual variance,
rather than separate group residual variances. The interaction test
differs from the common-covariate-effects model.

Standardized slopes multiply raw slopes by SD(predictor)/SD(outcome).
Both SDs come from all complete cases of the fitted model, pooled across
groups. For continuous moderators, simple slopes are evaluated at the
fitted-sample moderator mean minus/plus one SD. Standardization of the
interaction coefficient uses the same predictor/outcome SD ratio; it
does not multiply by the moderator SD. `SlopeDeltaStd` instead describes
the difference between the two evaluated standardized slopes.

`Results` contains one `Interaction` row and one `Slope: <level>` row
per evaluated slope. `Estimate`, `StdError`, `ConfLow`, `ConfHigh`, and
`PValue` describe raw-unit inference; `EstimateStd`, `StdErrorStd`,
`ConfLowStd`, and `ConfHighStd` use the common standardization.
Intervals are 95% t intervals from linear combinations of coefficients
and `vcov(model)`. `N` is model N; `NGroup`, `NReference`, and
`NComparison` report categorical group counts. `ModeratorValue` records
continuous moderator evaluation points.

Categorical moderators with more than two levels use a nested-model
joint F-test for all predictor-by-moderator terms. All level slopes are
returned. The plotted pair is the pair with the largest absolute slope
difference; factor-level order sets reference then comparison, with ties
resolved by the first pair. `ReferenceLevel`, `ComparisonLevel`, and
`LargestLevelDifference` identify this contrast. The interaction row's
estimate and CI describe that selected contrast, while `PValue` is the
joint F-test; `ContrastPValue` reports the contrast's t-test separately.
These contrast intervals are unadjusted for selecting the largest pair.
For two levels, the contrast is comparison minus reference. Factor
coding does not change the simple slopes or these contrasts.

The matrix fill is `SlopeDeltaStd`: comparison minus reference (high
minus low for continuous moderators), with a single diverging scale.
Positive values are blue, negative red, and zero white. Shared symmetric
limits are returned in `SlopeDifferenceLimits`. Only interaction
p-values \<= `Pthresh` are displayed: raw p-values in `Unadjusted`,
BH-adjusted p-values in `FDRCorrected`, regardless of `PVar`.
Nonsignificant cells are blank; failed models are grey and recorded with
warnings in `Metadata$Failed`. A no-interactions note is displayed when
no cell passes the threshold. Hatching marks strict sign reversals (zero
is not a reversal); stars describe interaction significance, not
significance of either slope.

Existing raw-unit matrices, p-value tables, aliases, and `raw_data`
columns are retained. `Slope1Std`, `Slope2Std`, and `SlopeDeltaStd`
matrices are added to both result branches. For multilevel moderators
these matrices describe the selected pair. `raw_data` retains
`StrengthChange` as a descriptive magnitude change; it no longer
controls fill. `StrengthLimits` is retained for compatibility and is
distinct from the new plotted limits. Simple-slope BH corrections form a
separate family from interaction tests, using the same `fdr_scope`.
`Metadata$PooledScaling` records each model's SDs.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

outcomes <- c("age", "ACE_CD143_Angiotensin_Converti",
              "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
              "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
              "Alpha_1_Microglobulin")
predictors <- c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
                "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E")

# With a categorical interaction variable (Diagnosis)
results <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "Diagnosis",
  outcome_vars = outcomes,
  predictor_vars = predictors
)

# With a continuous interaction variable (age)
results_cont <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "age",
  outcome_vars = setdiff(outcomes, "age"),
  predictor_vars = predictors
)

# Raw p-value interaction matrix
results$Unadjusted$plot


# FDR-adjusted interaction matrix
results$FDRCorrected$plot

# }
```
