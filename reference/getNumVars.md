# Get Numeric Variables

Extracts numeric variables from a data frame.

## Usage

``` r
getNumVars(data, Ordinal = FALSE, DataFrame = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame from which to extract numeric variables.

- Ordinal:

  Logical, indicating whether to include ordinal variables.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing the names of numeric variables.

## Details

The variable-selector helpers exist so that a variable set is derived
from the data rather than typed out and left to drift. `getNumVars()`,
[`getCatVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getCatVars.md),
and
[`getBinaryVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getBinaryVars.md)
partition a labelled data frame the way the analysis functions expect,
and each takes an `Ordinal` argument because an ordered factor can
legitimately count as either continuous or categorical depending on the
analysis.

Run them on a frame that has already been through
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md).
Before relabelling, a 0/1 diagnosis is still numeric and would be picked
up here as continuous.

## See also

[`getCatVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getCatVars.md),
[`getBinaryVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getBinaryVars.md),
and
[`ConvertOrdinalToNumeric()`](https://rdastgh1.github.io/SciDataReportR/reference/ConvertOrdinalToNumeric.md)
for the ordinal policy these share.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# 128 numeric columns: age plus the biomarker panel
vars_Numeric <- getNumVars(Labelled)
length(vars_Numeric)
#> [1] 128
utils::head(vars_Numeric)
#> [1] "age"                             "ACE_CD143_Angiotensin_Converti" 
#> [3] "ACTH_Adrenocorticotropic_Hormon" "AXL"                            
#> [5] "Adiponectin"                     "Alpha_1_Antichymotrypsin"       

# Ordinal variables are excluded by default; include them when they should
# be modelled on their numeric scale.
length(getNumVars(Labelled, Ordinal = TRUE))
#> [1] 128

# The point of deriving the set: it feeds straight into an analysis without
# a hand-typed vector that can fall out of date.
MakeTable1(Labelled, variables = utils::head(vars_Numeric, 5))


  

Characteristic
```
