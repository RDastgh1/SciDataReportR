# Example Dataset: SampleVariableTypes

An example of a modified VariableTypes file to be used to Revalue
SampleData.

## Usage

``` r
SampleVariableTypes
```

## Format

A data frame with 11 columns and 138 rows:

- Variable:

  Variable name, as listed in column name of data file.

- Label:

  Variable label, as desired in tables and figures.

- Type:

  How the variable should be treated, e.g., continuous (double) or
  categorical.

- Category:

  Optional for future filtering: the category type this variable belongs
  to.

- Recode:

  0 or 1, indicating whether or not this variable should be recoded.

- Code:

  If the variable should be recoded, how it should be recoded.

- Notes:

  Optional notes for the variable.

- Exclude:

  Optional for filtering: whether or not this variable should be
  excluded.

- Subcategory:

  Optional additional categories for the variable.

- Include:

  Optional for filtering: whether or not this variable should be
  included.

- MissingCode:

  Optional: if a value should be considered NA.

## Source

Exported from CreateVariableTypes(SampleData, "SampleData.csv") and
modified in Excel.

## Details

This is what an edited codebook looks like in practice. It began as the
output of
[`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md)
run on
[SampleData](https://rdastgh1.github.io/SciDataReportR/reference/SampleData.md),
and was then filled in by hand: labels written for the columns that
needed them, `sex` marked for recoding with its `Code` mapping, a
missing code recorded where one applies, and `Category` / `Subcategory`
/ `Include` used to group variables for later selection.

Most rows are left mostly blank, and that is normal. Of 138 variables
only a handful carry a recode or a missing code; the rest simply need a
readable label. A codebook is a place to record the exceptions, not a
form to complete.

Passing it with
[SampleData](https://rdastgh1.github.io/SciDataReportR/reference/SampleData.md)
to
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
is the first step of nearly every example in this package.

## See also

[SampleData](https://rdastgh1.github.io/SciDataReportR/reference/SampleData.md)
for the data it describes,
[`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md)
to generate one for your own data, and
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
to apply it.

## Examples

``` r
data(SampleVariableTypes)

dim(SampleVariableTypes)
#> [1] 138  11
table(SampleVariableTypes$Type)
#> 
#> Categorical      Double 
#>           8         130 

# \donttest{
# The whole codebook: one row per variable in SampleData
htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(
    SampleVariableTypes,
    height = "420px", full_width = TRUE
  )
)))


 Variable 
```
