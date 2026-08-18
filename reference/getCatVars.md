# Get Categorical Variables

Extracts categorical variables from a data frame.

## Usage

``` r
getCatVars(data, Ordinal = TRUE, DataFrame = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame from which to extract categorical variables.

- Ordinal:

  Logical, indicating whether to include ordinal variables.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing the names of categorical variables.

## See also

[`getNumVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getNumVars.md)
and
[`getBinaryVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getBinaryVars.md)
for the other partitions.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Every factor in the frame
getCatVars(Labelled)
#> [1] "Diagnosis" "sex"       "Genotype" 

# `Ordinal = FALSE` drops ordered factors, which is what you want when the
# ordered variables are going to be analyzed on their numeric scale instead.
getCatVars(Labelled, Ordinal = FALSE)
#> [1] "Diagnosis" "sex"       "Genotype" 

# Only meaningful after RevalueData(): in the raw extract `sex` is still a
# bare 0/1 numeric column and is not detected as categorical.
getCatVars(SampleData)
#> character(0)
```
