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

  A data frame containing the variable codebook. It should have columns
  `Variable` and `MissingCode`, where `Variable` specifies the variable
  name in the data frame and `MissingCode` specifies the code to be
  replaced with `NA`.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- VariableCodebook:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

A data frame with specified missing codes replaced by `NA`.

## Examples

``` r
# A data frame that uses sentinel values to encode missingness
df <- data.frame(
  id    = 1:6,
  age   = c(34, 999, 52, 999, 41, 29),
  score = c(10, -9, -9, 15, 20, 12)
)

# A codebook mapping each variable to its missing code(s)
codebook <- data.frame(
  Variable    = c("age", "score"),
  MissingCode = c("999", "-9")
)

# Before: sentinel codes still present
df
#>   id age score
#> 1  1  34    10
#> 2  2 999    -9
#> 3  3  52    -9
#> 4  4 999    15
#> 5  5  41    20
#> 6  6  29    12

# After: sentinel codes replaced with NA
ReplaceMissingCode(df, codebook)
#>   id age score
#> 1  1  34    10
#> 2  2  NA    NA
#> 3  3  52    NA
#> 4  4  NA    15
#> 5  5  41    20
#> 6  6  29    12
```
