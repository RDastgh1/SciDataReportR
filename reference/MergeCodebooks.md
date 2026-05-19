# Merge multiple codebooks using harmonization rules

Deterministically merge multiple codebooks using optional harmonization
rules generated from
[`CodebookMergeApp()`](https://rdastgh1.github.io/SciDataReportR/reference/CodebookMergeApp.md).

## Usage

``` r
MergeCodebooks(
  codebooks,
  Rules = NULL,
  VariableCol = "Variable",
  warn = TRUE,
  strict = FALSE
)
```

## Arguments

- codebooks:

  Named list of codebook data frames.

- Rules:

  Optional harmonization rules generated from
  [`CodebookMergeApp()`](https://rdastgh1.github.io/SciDataReportR/reference/CodebookMergeApp.md).

- VariableCol:

  Name of variable identifier column.

- warn:

  Logical; emit warnings.

- strict:

  Logical; stop on unresolved conflicts.

## Value

A list containing:

- Codebook:

  Merged harmonized codebook

- ConflictReport:

  Detected conflicts

- AppliedRules:

  Applied rules

## Details

If the supplied rules match the current codebooks and conflict
structure, warnings are suppressed because the harmonization has already
been reviewed.

Warnings are only emitted when:

- new conflicts appear

- codebooks differ from those used to generate rules

- conflict structure changes

## Examples

``` r
if (FALSE) { # \dontrun{

MergedCB <- MergeCodebooks(

  codebooks = list(
    Study1 = cb1,
    Study2 = cb2
  ),

  Rules = MergeRules
)

} # }
```
