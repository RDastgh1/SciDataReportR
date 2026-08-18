# Create a merged gtsummary table by faceting comparisons across multiple categorical variables

Generates a series of comparison tables using
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
for each categorical variable (facet) in the provided list and merges
them side-by-side using
[`gtsummary::tbl_merge()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_merge.html).
This function extends the functionality of
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
by automatically detecting which facet variables are categorical (factor
or character) and producing a faceted summary of how the main comparison
variable (e.g., Cluster, TreatmentArm) differs across multiple
categorical dimensions such as Race, Sex, or HIV status.

## Usage

``` r
MakeFacetCatComparisonTable(
  data,
  FacetVariables,
  variables,
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
  CompFun = MakeComparisonTable,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  ...,
  DataFrame = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  Covariates = lifecycle::deprecated(),
  ValueDigits = lifecycle::deprecated(),
  pDigits = lifecycle::deprecated(),
  EffectSizeDigits = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing all variables to be analyzed.

- FacetVariables:

  A character vector of variable names to facet by. The function
  automatically selects those that are categorical (`factor` or
  `character`).

- variables:

  A character string naming the variable(s) being compared (e.g.,
  "Cluster").

- covariates:

  Optional character vector of covariate names to adjust for.

- value_digits:

  Number of decimal digits to display for numeric values (default = 2).

- p_digits:

  Number of decimal digits to display for p-values (default = 3).

- AddEffectSize:

  Logical; if TRUE, include effect sizes (default = FALSE).

- effect_size_digits:

  Decimal digits for effect size values (default = 2).

- AddPairwise:

  Logical; if TRUE, include pairwise comparisons (default = FALSE).

- PairwiseMethod:

  Method for pairwise comparison p-value adjustment (default =
  "bonferroni").

- Parametric:

  Logical; if TRUE, use parametric tests (default = TRUE).

- ParametricDisplay:

  Optional vector specifying which statistics to display for parametric
  tests.

- IncludeOverallN:

  Logical; if TRUE, adds overall N to the table (default = FALSE).

- IncludeMissing:

  Logical; if TRUE, includes missing categories (default = FALSE).

- suppress_warnings:

  Logical; suppress internal warnings (default = FALSE).

- Referent:

  Optional string specifying the referent category for binary or
  categorical comparisons.

- IncludeOverallStats:

  Logical; if TRUE, adds overall descriptive statistics (default =
  FALSE).

- ShowPositiveBinaryOnLabel:

  Logical; if TRUE, labels binary variables with positive outcome
  (default = TRUE).

- CompFun:

  Comparison function to apply; defaults to `MakeComparisonTable`.

- Relabel:

  Logical; if TRUE (default), use attached variable labels.

- TreatOrdinalAs:

  How ordinal variables are treated in each table.

- ...:

  Additional arguments passed to the comparison function.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

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

A `gtsummary` table created by merging each facet's
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
output side-by-side using
[`gtsummary::tbl_merge()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_merge.html).
Each facet variable is labeled with its own tab spanner header for
clarity.

## Details

MakeFacetCatComparisonTable

Use this when the same set of measures needs to be compared across
several different groupings at once - by diagnosis, and by sex, and by
APOE status - and those comparisons belong in one table rather than
three.

Running
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
once per grouping produces tables that are each correct but cannot be
read against each other: the variables are listed three times, and
spotting that a biomarker separates diagnosis groups but not sexes means
looking back and forth between pages. Merging them puts one row per
variable and one column block per grouping, so that comparison is a
left-to-right read.

Facet variables that are not categorical are dropped automatically, so a
vector of candidate groupings can be passed without pre-filtering it.

Every argument that controls the individual tables - `AddPairwise`,
`covariates`, `Parametric`, and the rest - is passed through to each
facet, so the blocks stay consistent with one another.

## See also

[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
for a single grouping, and
[`MakeTable1()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeTable1.md)
for a plain descriptive table with no grouping at all.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

vars_Compare <- c("age", "AXL", "Adiponectin", "tau", "p_tau")

# The same five measures compared across two groupings at once
MakeFacetCatComparisonTable(
  data = Labelled,
  FacetVariables = c("Diagnosis", "sex"),
  variables = vars_Compare
)


  
Comparison by Diagnosis (values: mean (SD)).

  

Characteristic
```
