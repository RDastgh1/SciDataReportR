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
  DataFrame,
  FacetVariables,
  Variables,
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
  CompFun = MakeComparisonTable,
  ...
)
```

## Arguments

- DataFrame:

  A data frame containing all variables to be analyzed.

- FacetVariables:

  A character vector of variable names to facet by. The function
  automatically selects those that are categorical (`factor` or
  `character`).

- Variables:

  A character string naming the variable(s) being compared (e.g.,
  "Cluster").

- Covariates:

  Optional character vector of covariate names to adjust for.

- ValueDigits:

  Number of decimal digits to display for numeric values (default = 2).

- pDigits:

  Number of decimal digits to display for p-values (default = 3).

- AddEffectSize:

  Logical; if TRUE, include effect sizes (default = FALSE).

- EffectSizeDigits:

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

- ...:

  Additional arguments passed to the comparison function.

## Value

A `gtsummary` table created by merging each facet’s
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
output side-by-side using
[`gtsummary::tbl_merge()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_merge.html).
Each facet variable is labeled with its own tab spanner header for
clarity.

## Details

MakeFacetCatComparisonTable
