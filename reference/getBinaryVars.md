# Identify Binary Variables

This function identifies and returns a list of binary variables in a
dataframe. Binary variables are defined as having exactly two unique
values or levels. The function supports options for handling ordinal
factors and revalued data.

## Usage

``` r
getBinaryVars(
  data,
  Ordinal = TRUE,
  Revalued = TRUE,
  DataFrame = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe to analyze for binary variables.

- Ordinal:

  Logical. If TRUE, ordinal factors are included in the search for
  binary variables. Default is TRUE.

- Revalued:

  Logical. If TRUE, the function checks factors and their levels;
  otherwise, it checks for variables with two unique values.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing the names of binary variables.

## See also

[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md)
to fix which level counts as positive, and
[`getCatVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getCatVars.md)
/
[`getNumVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getNumVars.md)
for the other partitions.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Two-level factors, which are the ones that can be modelled as 0/1
vars_Binary <- getBinaryVars(Labelled)
vars_Binary
#> [1] "Diagnosis" "sex"      

# `Revalued = FALSE` looks for any column with two distinct values instead
# of two factor levels, for frames that have not been through RevalueData().
getBinaryVars(SampleData, Revalued = FALSE)
#> [1] "Diagnosis" "sex"      

# Which level each one is scored against
createBinaryMapping(Labelled, vars_Binary)
#>    Variable     Label PositiveLevel NegativeLevel
#> 1 Diagnosis Diagnosis      Impaired       Control
#> 2       sex       Sex          Male        Female
```
