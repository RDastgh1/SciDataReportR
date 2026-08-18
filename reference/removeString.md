# Remove Strings from a Vector

This function removes strings from a vector that are present in another
vector.

## Usage

``` r
removeString(Orig, Remove)
```

## Arguments

- Orig:

  The original vector containing strings.

- Remove:

  The vector of strings to be removed from the original vector.

## Value

A vector containing the strings from the original vector that were not
present in the removal vector.

## Details

A set-difference helper for variable vectors, equivalent to
`setdiff(Orig, Remove)` but keeping duplicates in `Orig`. It exists for
the common case of taking an auto-detected variable set and dropping the
few columns that should not be analyzed - identifiers, dates, the
outcome itself - without retyping the ones that should.

Names in `Remove` that do not appear in `Orig` are ignored silently, so
a stale exclusion list will not error; it will simply stop having an
effect.

## See also

[`getNumVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getNumVars.md),
[`getCatVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getCatVars.md),
and
[`getBinaryVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getBinaryVars.md),
which produce the vectors this usually trims.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Every numeric column, including ones that are not really measurements
vars_Numeric <- getNumVars(Labelled)
length(vars_Numeric)
#> [1] 128

# Drop the ones that should not be screened as biomarkers
vars_Biomarkers <- removeString(vars_Numeric, c("age", "tau", "p_tau"))
length(vars_Biomarkers)
#> [1] 125
utils::head(vars_Biomarkers, 4)
#> [1] "ACE_CD143_Angiotensin_Converti"  "ACTH_Adrenocorticotropic_Hormon"
#> [3] "AXL"                             "Adiponectin"                    

# Names that are not present are ignored, so an over-broad exclusion list
# is harmless
removeString(c("age", "sex", "AXL"), c("sex", "not_a_column"))
#> [1] "age" "AXL"

# Unlike setdiff(), duplicates in the original are kept
removeString(c("a", "a", "b", "c"), "b")
#> [1] "a" "a" "c"
setdiff(c("a", "a", "b", "c"), "b")
#> [1] "a" "c"
```
