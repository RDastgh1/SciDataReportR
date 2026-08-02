# Create Summary Table using gtsummary

This function is a wrapper around
[`gtsummary::tbl_summary`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
that ensures continuous variables are treated as continuous.

## Usage

``` r
MakeTable1(
  data,
  variables = NULL,
  TreatOrdinalAs = "Categorical",
  Relabel = TRUE,
  AutoDetectDistribution = FALSE,
  IncludeMissing = "ifany",
  DataFrame = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataframe to create the summary table from.

- variables:

  Optional. A character vector specifying the names of variables to
  include in the summary table. If NULL, all variables are included.

- TreatOrdinalAs:

  Character. Specifies how ordinal variables should be treated. Can be
  "Continuous", "Categorical", or "Both".

- Relabel:

  Logical; if TRUE (default), use attached variable labels.

- AutoDetectDistribution:

  Logical. If TRUE, the function will attempt to automatically detect
  the distribution of variables. Default is FALSE.

- IncludeMissing:

  Character matching gtsummary criteria. Can be "no", "ifany", or
  "always". Default is "ifany"

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A summary table created using gtsummary.

## References

This function wraps gtsummary. Please cite:

Sjoberg, D. D., Whiting, K., Curry, M., Lavery, J. A., & Larmarange, J.
(2021). Reproducible summary tables with the gtsummary package. *The R
Journal*, 13(1), 570-580.
[doi:10.32614/RJ-2021-053](https://doi.org/10.32614/RJ-2021-053)

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels so the table is publication-ready
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A Table 1 across 16 mixed-type variables
vars <- c(
  "Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin",
  "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin", "Apolipoprotein_A1",
  "Apolipoprotein_B", "C_Reactive_Protein", "Cortisol", "Cystatin_C",
  "Ferritin", "Insulin", "Leptin"
)

MakeTable1(Labelled, variables = vars)


  

Characteristic
```
