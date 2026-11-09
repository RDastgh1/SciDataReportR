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
  codebook = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  Pthresh = 0.05,
  PVar = "PValue"
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

- cluster_rows, cluster_columns:

  Logical; cluster heatmap rows and columns.

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

A named list with `PerformanceTable`, `RegressionTable`,
`ThresholdTable`, `FailureTable`, `Evaluations`, `Plots`, and
`Metadata`. For binary outcomes, `Plots` also includes `BiomarkerPanels`
(raw outcome-stratified distributions for continuous biomarkers) and
`ROCFacets` (adjusted ROC curves), each annotated with pair-specific
performance.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

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
