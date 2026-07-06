# Create normative T-scores from a regression model

Compatibility alias for
[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md).
Prefer
[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md)
in new code because this function fits a reusable normative T-score
model.

## Usage

``` r
CreateNormativeTScores(...)
```

## Arguments

- ...:

  Arguments passed to
  [`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md).

## Value

The same normative model object returned by
[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md).

## See also

[`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md)
for the canonical function and full examples.

## Examples

``` r
df <- tibble::tibble(
  Group = c(rep("Reference", 8), rep("Clinical", 2)),
  Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
  Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
  TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49)
)

out <- CreateNormativeTScores(
  data = df,
  test_var = "TrailsA",
  count_var = "Visit",
  covariates = "Age",
  reference_var = "Group",
  reference_value = "Reference",
  return_plots = TRUE
)
out$plots$tscore
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
