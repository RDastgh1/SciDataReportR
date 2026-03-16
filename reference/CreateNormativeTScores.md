# Create normative T-scores from a regression model

Fits a normative regression model in a user-defined reference subgroup
and uses the model residual standard deviation to convert observed
scores into z-scores and T-scores. This is useful for creating
demographically adjusted cognitive norms with optional practice effect
adjustment and optional preprocessing such as unit conversion, log
transformation, and reverse scoring.

## Usage

``` r
CreateNormativeTScores(
  df,
  test_var,
  count_var,
  covariates,
  reference_var,
  reference_value,
  include_practice_effect = FALSE,
  baseline_count_value = 1,
  reverse_score = FALSE,
  convert_seconds = FALSE,
  seconds_divisor = 1000,
  log_transform = TRUE,
  codebook = NULL,
  return_plots = TRUE
)
```

## Arguments

- df:

  A data frame containing the test variable, count variable, reference
  group variable, and covariates.

- test_var:

  A character string naming the raw test score variable.

- count_var:

  A character string naming the visit count or practice count variable.

- covariates:

  A character vector of covariate column names to include in the
  normative model.

- reference_var:

  A character string naming the variable used to define the normative
  reference group.

- reference_value:

  The value of `reference_var` that defines the normative reference
  group.

- include_practice_effect:

  Logical. If `TRUE`, `count_var` is included as a predictor and the
  model is fit using all available visits in the reference group. If
  `FALSE`, the model is fit only on rows where
  `count_var == baseline_count_value`.

- baseline_count_value:

  The value of `count_var` used to define the baseline visit when
  `include_practice_effect = FALSE`. Defaults to `1`.

- reverse_score:

  Logical. If `TRUE`, the analysis-scale score is multiplied by `-1` so
  that higher values reflect better performance.

- convert_seconds:

  Logical. If `TRUE`, the raw score is divided by `seconds_divisor`
  before further processing.

- seconds_divisor:

  Numeric divisor used when `convert_seconds = TRUE`. Defaults to
  `1000`.

- log_transform:

  Logical. If `TRUE`, applies
  [`log10()`](https://rdrr.io/r/base/Log.html) to the analysis score
  after optional unit conversion and before optional reverse scoring.

- codebook:

  Optional data frame with columns `Variable` and `Label`. If supplied,
  plot labels use variable labels when available.

- return_plots:

  Logical. If `TRUE`, returns a list of diagnostic plots.

## Value

A list with the following elements:

- data:

  A tibble containing the original data plus `NormRaw`, `NormScaled`,
  `NormPredicted`, `NormZ`, and `NormT`.

- model:

  The fitted `lm` object.

- model_summary:

  A tibble of coefficient estimates.

- model_fit:

  A one-row tibble containing model fit statistics.

- training_data:

  The rows used to fit the normative model.

- plots:

  A named list of ggplot objects when `return_plots = TRUE`.

- settings:

  A list of preprocessing and modeling settings used.

## Details

The reference group is defined by `reference_var == reference_value`.
When `include_practice_effect = FALSE`, the normative model is fit only
on rows where `count_var == baseline_count_value`. When
`include_practice_effect = TRUE`, `count_var` is added to the model and
all available visits in the reference group are used.

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

out <- CreateNormativeTScores(
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

out$data
#> # A tibble: 10 × 12
#>    Group       Age Education Sex   Visit TrailsA NormRaw NormCount NormScaled
#>    <chr>     <dbl> <fct>     <fct> <dbl>   <dbl>   <dbl>     <dbl>      <dbl>
#>  1 Reference    30 College   F         1   35000   35000         1      -1.54
#>  2 Reference    34 College   M         1   38000   38000         1      -1.58
#>  3 Reference    38 Graduate  F         1   40000   40000         1      -1.60
#>  4 Reference    42 Graduate  M         1   43000   43000         1      -1.63
#>  5 Reference    46 College   F         2   36000   36000         2      -1.56
#>  6 Reference    50 Graduate  M         2   39000   39000         2      -1.59
#>  7 Reference    54 College   F         2   41000   41000         2      -1.61
#>  8 Reference    58 Graduate  M         2   44000   44000         2      -1.64
#>  9 Clinical     40 College   F         1   47000   47000         1      -1.67
#> 10 Clinical     52 Graduate  M         2   49000   49000         2      -1.69
#> # ℹ 3 more variables: NormPredicted <dbl>, NormZ <dbl>, NormT <dbl>
out$model
#> 
#> Call:
#> stats::lm(formula = model_formula, data = training_data)
#> 
#> Coefficients:
#>       (Intercept)                Age  EducationGraduate               SexM  
#>        -1.4384063         -0.0068747         -0.0002786         -0.0055008  
#>             Visit  
#>         0.0989387  
#> 
```
