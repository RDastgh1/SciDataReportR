# Evaluate biomarker performance

Evaluates a continuous or categorical biomarker against a binary or
continuous outcome. Binary analyses use ordinary logistic regression.
Models with separation, non-convergence, or aliased coefficients return
stable unavailable metrics.

## Usage

``` r
EvaluateBiomarkerPerformance(
  data,
  outcome_var,
  biomarker_var,
  covariates = NULL,
  PositiveLevel = NULL,
  OutcomeType = c("auto", "binary", "continuous"),
  ThresholdMethod = c("youden", "sensitivity", "specificity", "custom"),
  ThresholdValue = NULL,
  RawThresholdValue = NULL,
  ProbabilityThresholdValue = NULL,
  Validation = c("none", "bootstrap", "cross_validation"),
  BootstrapR = 500,
  CVFolds = 10,
  CIBootstrapR = 500,
  CILevel = 0.95,
  CalibrationGroups = 10,
  Seed = 123,
  Relabel = TRUE,
  codebook = NULL,
  Verbose = TRUE
)
```

## Arguments

- data:

  A data frame.

- outcome_var:

  One outcome variable name.

- biomarker_var:

  One biomarker variable name.

- covariates:

  Optional covariate variable names.

- PositiveLevel:

  Positive binary outcome level, or `NULL` to use the second observed
  level.

- OutcomeType:

  One of `"auto"`, `"binary"`, or `"continuous"`.

- ThresholdMethod:

  One of `"youden"`, `"sensitivity"`, `"specificity"`, or `"custom"`.

- ThresholdValue:

  Target sensitivity or specificity.

- RawThresholdValue:

  Custom raw-biomarker threshold.

- ProbabilityThresholdValue:

  Custom predicted-probability threshold.

- Validation:

  One of `"none"`, `"bootstrap"`, or `"cross_validation"`.

- BootstrapR:

  Number of bootstrap optimism-correction resamples.

- CVFolds:

  Number of cross-validation folds.

- CIBootstrapR:

  Number of bootstrap confidence-interval resamples.

- CILevel:

  Confidence level.

- CalibrationGroups:

  Maximum grouped-calibration bins.

- Seed:

  Random seed.

- Relabel:

  Use codebook labels, then label attributes, for presentation.

- codebook:

  Optional data frame with `Variable` and `Label`.

- Verbose:

  Print positive-level information.

## Value

A stable named list containing models, performance, thresholds,
predictions, calibration, validation, plots, and metadata.

## Examples

``` r
if (FALSE) EvaluateBiomarkerPerformance(df, "DiseaseCohort", "NfL", c("Age", "Sex")) # \dontrun{}
```
