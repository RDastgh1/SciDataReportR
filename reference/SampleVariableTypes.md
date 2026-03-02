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

## Examples

``` r
data(SampleVariableTypes)
head(SampleVariableTypes)
#>    Variable          Label        Type Category Recode                 Code
#> 1 Diagnosis      Diagnosis Categorical     <NA>   <NA>                 <NA>
#> 2       age            Age      Double     <NA>   <NA>                 <NA>
#> 3       sex            Sex Categorical     <NA>    yes 0 = Female; 1 = Male
#> 4     Group          Group Categorical     <NA>   <NA>                 <NA>
#> 5      Race           Race Categorical     <NA>   <NA>                 <NA>
#> 6   Marital Marital Status Categorical     <NA>   <NA>                 <NA>
#>   Notes Exclude Subcategory Include MissingCode
#> 1    NA      NA                   1          NA
#> 2    NA      NA                   1         999
#> 3    NA      NA                   1          NA
#> 4    NA      NA                  NA          NA
#> 5    NA      NA                  NA          NA
#> 6    NA      NA                  NA          NA
```
