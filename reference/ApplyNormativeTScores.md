# Apply a normative T-score model to new data

Applies a previously fitted normative regression model to new data and
computes predicted values, z-scores, and T-scores using the same
preprocessing settings used during model development.

## Usage

``` r
ApplyNormativeTScores(df, normative_obj, score_prefix = "Norm")
```

## Arguments

- df:

  A data frame containing the test variable, count variable, and all
  predictors required by the normative model.

- normative_obj:

  A list returned by
  [`CreateNormativeTScores()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScores.md).

- score_prefix:

  A character string prefix used when naming output columns. Defaults to
  `"Norm"`.

## Value

A tibble containing the original data plus scored columns:

- Raw:

  The raw input score.

- Scaled:

  The transformed analysis-scale score.

- Predicted:

  The predicted score from the normative model.

- Z:

  The z-score.

- T:

  The T-score.

## Details

This function is designed to work with the output of
[`CreateNormativeTScores()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScores.md).
It uses the saved model and preprocessing settings to score new
observations consistently.

## Examples

``` r
df <- tibble::tibble(
  Group = c(
    rep("Reference", 8),
    rep("Clinical", 2)
  ),
  Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
  Education = factor(c(
    "College", "College", "Graduate", "Graduate",
    "College", "Graduate", "College", "Graduate",
    "College", "Graduate"
  )),
  Sex = factor(c(
    "F", "M", "F", "M", "F", "M", "F", "M", "F", "M"
  )),
  Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
  TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49) * 1000
)

norm_obj <- CreateNormativeTScores(
  df = df,
  test_var = "TrailsA",
  count_var = "Visit",
  covariates = c("Age", "Education", "Sex"),
  reference_var = "Group",
  reference_value = "Reference",
  include_practice_effect = TRUE,
  reverse_score = TRUE,
  convert_seconds = TRUE,
  log_transform = TRUE,
  return_plots = FALSE
)

scored_df <- ApplyNormativeTScores(
  df = df,
  normative_obj = norm_obj
)
```
