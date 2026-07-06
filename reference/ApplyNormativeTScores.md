# Apply a normative T-score model to new data

Applies a previously fitted normative regression model to new data and
computes predicted values, z-scores, and T-scores using the same
preprocessing settings used during model development.

## Usage

``` r
ApplyNormativeTScores(
  data,
  normative_obj,
  score_prefix = "Norm",
  df = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing the test variable, count variable, and all
  predictors required by the normative model.

- normative_obj:

  A list returned by
  [`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md).

- score_prefix:

  A character string prefix used when naming output columns. Defaults to
  `"Norm"`.

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A tibble containing the original data plus scored columns, each named
using the `score_prefix` string as a prefix:

- `{score_prefix}Raw`:

  The raw input score.

- `{score_prefix}Scaled`:

  The transformed analysis-scale score.

- `{score_prefix}Predicted`:

  The predicted score from the normative model.

- `{score_prefix}Z`:

  The z-score.

- `{score_prefix}T`:

  The T-score.

## Details

This function is designed to work with the output of
[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md).
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

norm_obj <- CreateNormativeTScoreModel(
  data = df,
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
  data = df,
  normative_obj = norm_obj
)
```
