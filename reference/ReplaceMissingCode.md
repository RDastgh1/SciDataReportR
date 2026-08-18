# Replace Missing Codes with NA

This function replaces specified missing codes in a data frame with `NA`
values based on a given variable codebook.

## Usage

``` r
ReplaceMissingCode(
  data,
  codebook,
  DataFrame = lifecycle::deprecated(),
  VariableCodebook = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing the data.

- codebook:

  A data frame containing the variable codebook. It must have columns
  `Variable` and `MissingCode`. `Variable` names the column in `data`;
  `MissingCode` holds the code, or several codes separated by commas or
  semicolons, to replace with `NA`. Rows with a missing or blank
  `MissingCode` are skipped.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- VariableCodebook:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

A data frame with specified missing codes replaced by `NA`.

## Details

Study databases rarely leave missing values empty. They record *why* the
value is missing, using sentinel numbers outside the plausible range -
`-7` for "not applicable", `-8` for "refused", `-9` for "not asked",
`999` for "unknown". Left in place, those codes are read as real
measurements: a mean age climbs into the hundreds, and a variable coded
`-7`/`-8` gets a negative mean that no one notices until much later.

This converts every code listed in the codebook to `NA`, so the values
are excluded from statistics and counted as missing.

## Several markers for one variable

A variable usually has more than one missing code, since the whole point
of sentinels is to distinguish reasons. List them all in a single
`MissingCode` cell, separated by commas or semicolons:

    Variable   MissingCode
    age        -7, -8, -9
    income     999; 9999

Whitespace around the separators is ignored, so `-7,-8` and `-7, -8`
behave identically. When a codebook is assembled in R rather than read
from a spreadsheet, `MissingCode` may also be a list-column of vectors.

Codes are matched on the column's own scale: numerically for numeric
columns, and as text for character and factor columns, so markers like
`"Refused"` or `"Unknown"` work. When such a marker is a factor level,
the level is removed along with the values, so it does not linger as an
empty category in later tables.

Variables whose `MissingCode` is blank are left untouched, and a
codebook naming variables that are not in `data` warns rather than
failing - a codebook usually describes more columns than any one
analysis selects.

## See also

[`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md),
which generates a codebook with a `MissingCode` column ready to fill in,
and
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md),
which applies missing codes as part of the full relabelling workflow.

## Examples

``` r
# `age` uses three different codes to record three different reasons
df <- data.frame(
  id     = 1:6,
  age    = c(34, 999, 52, -7, 41, -8),
  score  = c(10, -9, -9, 15, 20, 12),
  status = c("Active", "Unknown", "Active", "Withdrawn", "Unknown", "Active")
)

# Several markers for one variable go in one cell
codebook <- data.frame(
  Variable    = c("age", "score", "status"),
  MissingCode = c("999, -7, -8", "-9", "Unknown")
)

# Before: the sentinels are averaged in as if they were ages
df
#>   id age score    status
#> 1  1  34    10    Active
#> 2  2 999    -9   Unknown
#> 3  3  52    -9    Active
#> 4  4  -7    15 Withdrawn
#> 5  5  41    20   Unknown
#> 6  6  -8    12    Active
mean(df$age)
#> [1] 185.1667

# After: every listed code becomes NA
cleaned <- ReplaceMissingCode(df, codebook)
cleaned
#>   id age score    status
#> 1  1  34    10    Active
#> 2  2  NA    NA      <NA>
#> 3  3  52    NA    Active
#> 4  4  NA    15 Withdrawn
#> 5  5  41    20      <NA>
#> 6  6  NA    12    Active
mean(cleaned$age, na.rm = TRUE)
#> [1] 42.33333
colSums(is.na(cleaned))
#>     id    age  score status 
#>      0      3      2      2 

# Blank codes are skipped
ReplaceMissingCode(
  df,
  data.frame(Variable = c("id", "age"), MissingCode = c(NA, "999, -7, -8"))
)
#>   id age score    status
#> 1  1  34    10    Active
#> 2  2  NA    -9   Unknown
#> 3  3  52    -9    Active
#> 4  4  NA    15 Withdrawn
#> 5  5  41    20   Unknown
#> 6  6  NA    12    Active
```
