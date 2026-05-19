# Interactive codebook harmonization dashboard

Launch a Shiny dashboard for reviewing and harmonizing multiple
codebooks before deterministic merging with MergeCodebooks().

## Usage

``` r
CodebookMergeApp(
  codebooks,
  VariableCol = "Variable",
  auto_type_mapping = TRUE,
  ignore_columns = NULL
)
```

## Arguments

- codebooks:

  Named list of codebook data frames.

- VariableCol:

  Name of variable identifier column.

- auto_type_mapping:

  Logical; normalize common type synonyms.

- ignore_columns:

  Optional metadata columns to ignore.

## Value

Launches a Shiny app.
