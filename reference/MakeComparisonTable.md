# Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts

Create a label-aware comparison table using
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
with optional global hypothesis tests, covariate-adjusted tests, effect
sizes, and pairwise comparisons.

## Usage

``` r
MakeComparisonTable(
  DataFrame,
  CompVariable = NULL,
  Variables,
  ...,
  Covariates = NULL,
  ValueDigits = 2,
  pDigits = 3,
  AddEffectSize = FALSE,
  EffectSizeDigits = 2,
  AddPairwise = FALSE,
  PairwiseMethod = "bonferroni",
  Parametric = TRUE,
  ParametricDisplay = NULL,
  IncludeOverallN = FALSE,
  IncludeMissing = FALSE,
  suppress_warnings = FALSE,
  Referent = NULL,
  IncludeOverallStats = FALSE,
  ShowPositiveBinaryOnLabel = TRUE,
  CatMethod = c("auto", "chisq", "fisher"),
  MultiCatAdjusted = c("multinomial_LR", "none"),
  ShowNotes = c("auto", "always", "never"),
  NotesPosition = c("last", "after_test", "before_pairwise")
)
```

## Arguments

- DataFrame:

  A data frame.

- CompVariable:

  Character scalar naming the grouping variable.

- Variables:

  Character vector of variables to summarize.

- ...:

  Optional additional variable names supplied individually.

- Covariates:

  Optional character vector of covariates for adjusted models.

- ValueDigits:

  Number of digits for descriptive statistics.

- pDigits:

  Number of digits for p-values.

- AddEffectSize:

  Logical; add effect-size columns.

- EffectSizeDigits:

  Number of digits for effect sizes.

- AddPairwise:

  Logical; add pairwise comparison columns.

- PairwiseMethod:

  P-value adjustment method. Use `"none"` for no adjustment.

- Parametric:

  Logical; use parametric tests for continuous outcomes.

- ParametricDisplay:

  Logical; display continuous summaries as mean (SD). If `FALSE`,
  display median [IQR](https://rdrr.io/r/stats/IQR.html). Defaults to
  `Parametric`.

- IncludeOverallN:

  Logical; add N column.

- IncludeMissing:

  Logical; include missing rows in summaries.

- suppress_warnings:

  Logical; suppress selected gtsummary warnings.

- Referent:

  Optional reference group for pairwise comparisons.

- IncludeOverallStats:

  Logical; add overall summary column.

- ShowPositiveBinaryOnLabel:

  Logical; for binary variables, show only the positive level where
  identifiable.

- CatMethod:

  Categorical test method. One of `"auto"`, `"chisq"`, `"fisher"`.

- MultiCatAdjusted:

  Adjusted multicategory method. Currently `"multinomial_LR"` or
  `"none"`.

- ShowNotes:

  Whether to show Notes column. One of `"auto"`, `"always"`, or
  `"never"`.

- NotesPosition:

  Notes column position. One of `"last"`, `"after_test"`, or
  `"before_pairwise"`.

## Value

A `gtsummary` object.

## Details

Continuous variables are numeric variables with more than two unique
non-missing values. Numeric variables with exactly two unique values are
treated as dichotomous categorical variables.

When covariates are supplied, continuous outcomes are tested using
ANCOVA with Type II tests. If `Parametric = FALSE`, robust HC3
covariance is used for the group-level Wald test and adjusted pairwise
comparisons.

Binary categorical outcomes with covariates are tested using logistic
regression likelihood-ratio tests. Multicategory categorical outcomes
with covariates are tested using multinomial likelihood-ratio tests.

Pairwise comparisons preserve non-standard group labels and variable
names.

## Examples

``` r
if (FALSE) { # \dontrun{
MakeComparisonTable(
  DataFrame = mtcars,
  CompVariable = "am",
  Variables = c("mpg", "hp", "wt"),
  AddEffectSize = TRUE,
  AddPairwise = TRUE
)
} # }
```
