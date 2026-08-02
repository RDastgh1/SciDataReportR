# Add a new variable to a codebook

Adds one variable entry to a codebook while preserving its existing
schema. In addition to the standard codebook fields, named values
supplied through `...` populate user-defined columns. A new named `...`
column is added to the codebook (with `NA` for existing rows) and
produces a warning.

## Usage

``` r
AddToCodebook(
  codebook,
  VariableName,
  VariableLabel = NA,
  VariableType = NA,
  VariableCategory = NA,
  VariableRecode = NA,
  VariableCode = NA,
  VariableExclude = NA,
  VariableNotes = NA,
  CB = lifecycle::deprecated(),
  ...
)
```

## Arguments

- codebook:

  A data frame representing the codebook. It must contain a `Variable`
  column.

- VariableName:

  A single, non-missing, non-empty character variable name. It must not
  already appear in `codebook$Variable`.

- VariableLabel:

  A single label for the variable. Defaults to `VariableName` when `NA`.

- VariableType:

  A single value for the `Type` column.

- VariableCategory:

  A single value for the `Category` column.

- VariableRecode:

  A single value for the `Recode` column.

- VariableCode:

  A single value for the `Code` column.

- VariableExclude:

  A single value for the `Exclude` column.

- VariableNotes:

  A single value for the `Notes` column.

- CB:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

- ...:

  Named, single atomic values for user-defined codebook columns. Names
  matching existing columns populate them. New names create a column and
  warn. Standard fields (`Variable`, `Label`, `Type`, `Category`,
  `Recode`, `Code`, `Exclude`, and `Notes`) must be supplied through
  their corresponding formal arguments.

## Value

A data frame representing the updated codebook with the new variable
added.

## Details

The function uses the supplied codebook as the source of truth for coded
metadata. For every supplied field except `Variable`, `Label`, and
`Notes`, it warns when a non-missing value has not previously appeared
in that column or has a different storage type. These warnings are
advisory: the row is still added and R may promote the column type while
binding the new row.

## Examples

``` r
codebook <- data.frame(
  Variable = "sex", Label = "Sex assigned at birth", Type = "Categorical",
  Category = "Demographics", Recode = 1, Code = "0 = Female; 1 = Male",
  Exclude = FALSE, Notes = NA_character_, Domain = "Clinical"
)

# Populate an existing user-defined column.
AddToCodebook(
  codebook, "age", "Age at enrollment", "Double", "Demographics",
  Domain = "Clinical"
)
#> Warning: `Type` has not previously appeared in this column.
#>   Variable                 Label        Type     Category Recode
#> 1      sex Sex assigned at birth Categorical Demographics      1
#> 2      age     Age at enrollment      Double Demographics     NA
#>                   Code Exclude Notes   Domain
#> 1 0 = Female; 1 = Male   FALSE  <NA> Clinical
#> 2                 <NA>      NA  <NA> Clinical

# This adds `Source` and warns because it is a new column.
AddToCodebook(codebook, "site", Source = "REDCap")
#> Warning: `Source` is not an existing codebook column; adding it with NA for existing rows.
#>   Variable                 Label        Type     Category Recode
#> 1      sex Sex assigned at birth Categorical Demographics      1
#> 2     site                  site        <NA>         <NA>     NA
#>                   Code Exclude Notes   Domain Source
#> 1 0 = Female; 1 = Male   FALSE  <NA> Clinical   <NA>
#> 2                 <NA>      NA  <NA>     <NA> REDCap

# These values are retained but warn because they conflict with the old schema.
AddToCodebook(codebook, "participant_id", VariableRecode = 4,
              VariableExclude = 3)
#> Warning: `Recode` has not previously appeared in this column.
#> Warning: `Exclude` has not previously appeared in this column and has a different storage type or class from existing values.
#>         Variable                 Label        Type     Category Recode
#> 1            sex Sex assigned at birth Categorical Demographics      1
#> 2 participant_id        participant_id        <NA>         <NA>      4
#>                   Code Exclude Notes   Domain
#> 1 0 = Female; 1 = Male       0  <NA> Clinical
#> 2                 <NA>       3  <NA>     <NA>
```
