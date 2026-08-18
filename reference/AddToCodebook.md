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
# \donttest{
# An existing codebook, including the user-defined `Domain` column
codebook <- data.frame(
  Variable = c("sex", "visit", "mmse"),
  Label = c("Sex assigned at birth", "Study visit", "MMSE total score"),
  Type = c("Categorical", "Categorical", "Double"),
  Category = c("Demographics", "Design", "Cognition"),
  Recode = c(1, 0, 0),
  Code = c("0 = Female; 1 = Male", "1 = Baseline; 2 = Follow-up", NA),
  Exclude = c(FALSE, FALSE, FALSE),
  Notes = NA_character_,
  Domain = c("Clinical", "Study", "Clinical")
)

# Each call returns the updated codebook, so entries chain
codebook <- AddToCodebook(
  codebook, "age", "Age at enrollment", "Double", "Demographics",
  Domain = "Clinical"
)

# `Source` is a new column: created with a warning, back-filled with NA
codebook <- AddToCodebook(
  codebook, "site", "Enrolling site", "Categorical", "Design",
  Source = "REDCap"
)
#> Warning: `Source` is not an existing codebook column; adding it with NA for existing rows.

# Kept, but warned about: they conflict with the established schema
codebook <- AddToCodebook(
  codebook, "participant_id", "Participant identifier",
  VariableRecode = 4, VariableExclude = 3
)
#> Warning: `Recode` has not previously appeared in this column.
#> Warning: `Exclude` has not previously appeared in this column and has a different storage type or class from existing values.

# The result as a formatted table
htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(codebook, height = "320px", full_width = TRUE)
)))


 Variable 
```
