# Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts

Create a label-friendly comparison table (via **gtsummary**) summarizing
variables across groups, with optional global hypothesis tests,
covariate-adjusted tests, effect sizes, and optional pairwise
comparisons.

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
  MultiCatAdjusted = c("multinomial_LR", "none")
)
```

## Arguments

- DataFrame:

  A data frame.

- CompVariable:

  Grouping variable name (character scalar). If NULL or missing and
  `IncludeOverallStats = TRUE`, an overall-only table is returned.

- Variables:

  Character vector of variables to include. You may also pass additional
  variable names as unnamed arguments via `...` (convenience).

- ...:

  Optional additional variable names (character scalars). Useful when
  you accidentally typed `Variables="A", "B"` instead of
  `Variables=c("A","B")`.

- Covariates:

  Optional covariate names (character vector) used for adjusted models.
  Any covariates appearing in `Variables` are removed from the displayed
  table with a warning.

- ValueDigits:

  Digits for summary statistics (default 2).

- pDigits:

  Digits to display formatted p-values using
  [`format.pval()`](https://rdrr.io/r/base/format.pval.html).

- AddEffectSize:

  Add effect sizes column (default FALSE).

- EffectSizeDigits:

  Digits for effect sizes (default 2).

- AddPairwise:

  Add pairwise p-value columns (default FALSE).

- PairwiseMethod:

  P-adjustment method (default "bonferroni"). Use "none" for no
  adjustment. Must be one of "none" or
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).

- Parametric:

  If TRUE, use parametric tests where relevant. If FALSE and covariates
  are present, robust ANCOVA is used for continuous outcomes and robust
  covariance is used for continuous pairwise.

- ParametricDisplay:

  If TRUE show mean (SD); if FALSE show median
  [IQR](https://rdrr.io/r/stats/IQR.html). Defaults to `Parametric`.

- IncludeOverallN:

  If TRUE add overall N column.

- IncludeMissing:

  If TRUE include missing/unknown category rows in summaries.

- suppress_warnings:

  If TRUE suppress gtsummary warnings.

- Referent:

  Optional reference level for pairwise contrasts (character scalar). If
  set, pairwise comparisons are against this referent.

- IncludeOverallStats:

  If TRUE return overall-only summary (ignores grouping).

- ShowPositiveBinaryOnLabel:

  If TRUE, display only the “positive” level for binary categorical
  variables.

- CatMethod:

  Categorical test method: "chisq", "fisher", or "auto" (default).

- MultiCatAdjusted:

  Adjusted multi-category method when covariates are present: currently
  supports "multinomial_LR" (default) or "none".

## Value

A
[`gtsummary::tbl_summary`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
object with added columns in `table_body`.

## Details

This function is designed for publication tables where you want the
**Test** column to reflect the method used for the reported p-value(s),
while still leveraging gtsummary's formatting, labeling, and table
structure.

### Overview

`MakeComparisonTable()` is a wrapper around
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
that:

- creates group-wise descriptive summaries,

- computes global p-values (unadjusted or covariate-adjusted),

- optionally computes effect sizes,

- optionally computes pairwise p-values (with multiplicity control and
  optional referent),

- adds a **Notes** column to explain common edge cases (e.g.,
  complete-case filtering dropping a group level, causing missing
  adjusted pairwise contrasts).

The function aims to use gtsummary capabilities wherever possible
(summary, formatting, styling, captioning) and computes additional
statistics not provided directly by gtsummary (adjusted global tests,
adjusted pairwise, effect sizes).

### Variable typing rules

- Continuous: numeric variables with \> 2 unique non-missing values.

- Dichotomous numeric: numeric variables with exactly 2 unique
  non-missing values (treated as categorical).

- Categorical: factors, characters, logicals, and non-numeric variables
  (or dichotomous numeric).

### Global tests (p-values)

Global p-values are computed per variable, using complete-case data for
the required fields.

#### Continuous outcomes

- No covariates:

  - `Parametric = TRUE`:

    - 2 groups: Welch t-test
      ([`stats::t.test()`](https://rdrr.io/r/stats/t.test.html), unequal
      variances)

    - 3+ groups: one-way ANOVA
      ([`stats::aov()`](https://rdrr.io/r/stats/aov.html))

  - `Parametric = FALSE`:

    - 2 groups: Wilcoxon rank-sum
      ([`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html))

    - 3+ groups: Kruskal-Wallis
      ([`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html))

- With covariates:

  - `Parametric = TRUE`: ANCOVA via linear model
    ([`stats::lm()`](https://rdrr.io/r/stats/lm.html)) with Type II test
    for the grouping term (car::Anova(type = 2)).

  - `Parametric = FALSE`: robust ANCOVA via
    [`stats::lm()`](https://rdrr.io/r/stats/lm.html) plus HC3 robust
    covariance (sandwich::vcovHC(type = "HC3")) and a Wald F-test for
    the grouping term
    ([`car::linearHypothesis()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)).

#### Categorical outcomes (unadjusted)

- `CatMethod = "chisq"`: Pearson chi-squared with `correct = FALSE`
  ([`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)).

- `CatMethod = "fisher"`: Fisher's exact test
  ([`stats::fisher.test()`](https://rdrr.io/r/stats/fisher.test.html));
  for RxC tables, uses simulated p-value (B = 1e4).

- `CatMethod = "auto"` (default): uses chi-squared unless expected
  counts are small; then Fisher.

#### Categorical outcomes (adjusted; covariates present)

- Binary outcome: logistic regression likelihood ratio test (LR) using
  stats::glm(family = binomial) and stats::drop1(test = "Chisq").

- Multi-category outcome (3+ levels): multinomial LR (default) via
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) and
  an LR comparison of models with and without the grouping term.
  Controlled by `MultiCatAdjusted`, which defaults to
  `"multinomial_LR"`.

### Pairwise comparisons

Pairwise columns are added when `AddPairwise = TRUE`.

#### Which pairs are compared

- If `Referent` is `NULL`: all pairwise group comparisons.

- If `Referent` is set: all groups are compared against the referent
  (treatment-vs-control).

#### Multiplicity control

- `PairwiseMethod` defaults to `"bonferroni"` and may be any value in
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).

- Use `"none"` for no adjustment.

#### Continuous outcomes

- No covariates:

  - Parametric: stats::pairwise.t.test(pool.sd = FALSE)

  - Nonparametric:
    [`stats::pairwise.wilcox.test()`](https://rdrr.io/r/stats/pairwise.wilcox.test.html)

- With covariates:

  - Uses adjusted means via
    [`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
    on the ANCOVA model.

  - If `Parametric = FALSE`, the same ANCOVA model is used, but pairwise
    inference uses the HC3 robust covariance matrix supplied to emmeans
    (so adjusted pairwise still works).

#### Categorical outcomes

- No covariates: pairwise chi-squared or Fisher (per `CatMethod`) on the
  2-group subset.

- With covariates:

  - Binary outcome: for each pair, fit logistic models with and without
    group and use LR p-value.

  - Multi-category outcome: for each pair, fit multinomial models with
    and without group and use LR p-value.

  - Edge case: in a given pairwise subset, a multi-category outcome may
    collapse to 2 levels; in that case, the function automatically
    switches to logistic LR for that pair.

### Effect sizes

Effect sizes are provided when `AddEffectSize = TRUE`.

- Continuous, 2 groups, unadjusted: absolute Cohen's d (\|d\|) via
  [`effectsize::cohens_d()`](https://easystats.github.io/effectsize/reference/cohens_d.html).

- Continuous, 3+ groups, unadjusted parametric: eta-squared (η²) via
  [`effectsize::eta_squared()`](https://easystats.github.io/effectsize/reference/eta_squared.html).

- Continuous, adjusted: partial eta-squared (partial η²) from Type II
  ANOVA table.

- Continuous, nonparametric: epsilon-squared approximation from
  Kruskal-Wallis.

- Categorical: Cramer's V (uses DescTools if available; otherwise a
  chi-squared-based approximation).

### Captions

The caption describes:

- display statistic for continuous variables (mean SD vs median IQR),

- whether covariates were used (and how continuous outcomes were
  tested),

- categorical test selection (`CatMethod`),

- adjusted multi-category method (`MultiCatAdjusted`),

- pairwise inclusion and multiplicity control (`PairwiseMethod`).

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  Cluster = factor(sample(0:3, 120, TRUE)),
  Age = rnorm(120, 50, 10),
  Sex = factor(sample(c("F","M"), 120, TRUE)),
  wrat3 = rnorm(120, 95, 12),
  Education_years = rnorm(120, 14, 2),
  Race = factor(sample(c("A","B","D"), 120, TRUE)),
  Hypertension = factor(sample(c("No","Yes"), 120, TRUE))
)

MakeComparisonTable(
  df, "Cluster",
  Variables = c("Education_years", "Race", "Hypertension"),
  Covariates = c("Age","Sex","wrat3"),
  AddPairwise = TRUE,
  PairwiseMethod = "none",
  Parametric = FALSE
)
} # }
```
