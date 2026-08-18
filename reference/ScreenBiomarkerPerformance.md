# Screen biomarker performance across outcomes

Applies
[`EvaluateBiomarkerPerformance()`](https://rdastgh1.github.io/SciDataReportR/reference/EvaluateBiomarkerPerformance.md)
across many candidate biomarkers and one or more binary or continuous
outcomes. Each biomarker-outcome pair uses all complete observations
available for that pair and the requested covariates. The function
returns comparison tables and screening plots, including an
interactive-ready heatmap whose cells contain hover text.

## Usage

``` r
ScreenBiomarkerPerformance(
  data,
  outcome_vars,
  biomarker_vars,
  covariates = NULL,
  PositiveLevel = NULL,
  OutcomeType = c("auto", "binary", "continuous"),
  Validation = c("none", "bootstrap", "cross_validation"),
  BootstrapR = 500,
  CVFolds = 10,
  CIBootstrapR = 200,
  CILevel = 0.95,
  HeatmapMetric = "AdjustedAUC",
  Seed = 123,
  Relabel = TRUE,
  codebook = NULL
)
```

## Arguments

- data:

  A data frame.

- outcome_vars:

  Character vector of outcome variable names.

- biomarker_vars:

  Character vector of biomarker variable names.

- covariates:

  Optional character vector of covariate variable names.

- PositiveLevel:

  Optional positive level. Supply one value for all binary outcomes or a
  named character vector with names matching `outcome_vars`.

- OutcomeType:

  One of `"auto"`, `"binary"`, or `"continuous"`, applied to all
  outcomes unless auto-detection is used.

- Validation:

  One of `"none"`, `"bootstrap"`, or `"cross_validation"`.

- BootstrapR:

  Number of bootstrap resamples for internal validation.

- CVFolds:

  Number of cross-validation folds.

- CIBootstrapR:

  Number of bootstrap resamples for performance confidence intervals.

- CILevel:

  Confidence level. Default is `0.95`.

- HeatmapMetric:

  Performance metric shown by the heatmap fill. Default is
  `"AdjustedAUC"`. Common alternatives are `"AUC"`, `"DeltaAUC"`,
  `"AdjustedR2"`, and `"DeltaR2"`.

- Seed:

  Random seed used for resampling.

- Relabel:

  Logical indicating whether labels should be used when available.

- codebook:

  Optional data frame with columns `Variable` and `Label`.

## Value

A named list with `PerformanceTable`, `RegressionTable`,
`ThresholdTable`, `FailureTable`, `Evaluations`, `Plots`, and
`Metadata`. For binary outcomes, `Plots` also includes `BiomarkerPanels`
(raw outcome-stratified distributions for continuous biomarkers) and
`ROCFacets` (adjusted ROC curves), each annotated with pair-specific
performance.

## Examples

``` r
if (FALSE) { # \dontrun{
screen <- ScreenBiomarkerPerformance(
  data = df,
  outcome_vars = c("DiseaseCohort", "Progression"),
  biomarker_vars = c("NfL", "GFAP", "Acrocyanosis"),
  covariates = c("Age", "Sex")
)

screen$PerformanceTable
screen$Plots$Heatmap
screen$Plots$Heatmap %>% add_biomarker_values()
} # }
```
