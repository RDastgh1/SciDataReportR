# Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts

Create a label-aware comparison table using
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
with optional global hypothesis tests, covariate-adjusted tests, effect
sizes, and pairwise comparisons.

## Usage

``` r
MakeComparisonTable(
  data,
  group_var = NULL,
  variables,
  ...,
  covariates = NULL,
  value_digits = 2,
  p_digits = 3,
  AddEffectSize = FALSE,
  effect_size_digits = 2,
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
  NotesPosition = c("last", "after_test", "before_pairwise"),
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  DataFrame = lifecycle::deprecated(),
  CompVariable = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  Covariates = lifecycle::deprecated(),
  ValueDigits = lifecycle::deprecated(),
  pDigits = lifecycle::deprecated(),
  EffectSizeDigits = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame.

- group_var:

  Character scalar naming the grouping variable.

- variables:

  Character vector of variables to summarize.

- ...:

  Optional additional variable names supplied individually.

- covariates:

  Optional character vector of covariates for adjusted models.

- value_digits:

  Number of digits for descriptive statistics.

- p_digits:

  Number of digits for p-values.

- AddEffectSize:

  Logical; add effect-size columns.

- effect_size_digits:

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

  Whether to show the Analysis notes column. One of `"auto"`,
  `"always"`, or `"never"`.

- NotesPosition:

  Analysis notes column position. One of `"last"`, `"after_test"`, or
  `"before_pairwise"`.

- Relabel:

  Logical; if TRUE (default), use attached variable labels.

- TreatOrdinalAs:

  How ordinal variables are treated: `"Categorical"`, `"Continuous"`,
  `"Both"`, or `"Exclude"`.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- CompVariable:

  **Deprecated** (since 19.15.0). Use `group_var` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

- ValueDigits:

  **Deprecated** (since 19.15.0). Use `value_digits` instead.

- pDigits:

  **Deprecated** (since 19.15.0). Use `p_digits` instead.

- EffectSizeDigits:

  **Deprecated** (since 19.15.0). Use `effect_size_digits` instead.

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

For parametric continuous outcomes, effect sizes are reported as
absolute Cohen's d for two-group comparisons and Cohen's f for omnibus
comparisons involving more than two groups. With covariates, two-group
Cohen's d uses the estimated marginal mean difference divided by the
model residual standard deviation, and multi-group Cohen's f is
calculated from the Type II ANCOVA group effect. These effect-size
scales are not numerically equivalent.

Binary categorical outcomes with covariates are tested using logistic
regression likelihood-ratio tests. Multicategory categorical outcomes
with covariates are tested using multinomial likelihood-ratio tests.

Pairwise comparisons preserve non-standard group labels and variable
names.

## Choosing the options

The defaults produce one test per variable, chosen from that variable's
type: a t-test for continuous variables, a chi-squared test for
categorical ones. The remaining arguments each answer a specific
question the default table cannot.

**`AddEffectSize`** adds the standardized magnitude beside the p-value.
A p-value says whether a difference is detectable, not whether it is
large, and effect sizes are what make two significant rows comparable to
each other.

**`AddPairwise`** matters as soon as there are more than two groups. The
omnibus test only says "these groups are not all the same"; it never
says which pair differs. Pairwise contrasts answer that, corrected
across the contrasts so that hunting through them does not inflate the
error rate. `PairwiseMethod` selects the correction - Bonferroni by
default, `"fdr"` when there are many contrasts, `"none"` for exploratory
work. `Referent` compares every group against one reference level
instead of against each other, which is usually what a control group is
for.

**`covariates`** replaces the simple test with a model-based one that
holds the named variables constant. An unadjusted group difference in a
biomarker may only reflect that the groups differ in age or sex;
adjusting reports the group difference that remains.

**`Parametric = FALSE`** switches continuous comparisons to rank-based
tests (Wilcoxon, or Kruskal-Wallis for more than two groups), which is
the right choice for skewed measures or small, uneven groups.

**`IncludeOverallN`** and **`IncludeMissing`** add an overall column and
explicit missing-value counts, for a table that has to stand on its own
in a manuscript.

## References

This function wraps gtsummary. Please cite:

Sjoberg, D. D., Whiting, K., Curry, M., Lavery, J. A., & Larmarange, J.
(2021). Reproducible summary tables with the gtsummary package. *The R
Journal*, 13(1), 570-580.
[doi:10.32614/RJ-2021-053](https://doi.org/10.32614/RJ-2021-053)

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

vars_Compare <- c("age", "sex", "AXL", "Adiponectin", "tau", "p_tau")

# Two groups, default tests
MakeComparisonTable(
  data = Labelled,
  group_var = "Diagnosis",
  variables = vars_Compare
)


  
Comparison by Diagnosis (values: mean (SD)).

  
Characteristic
```
