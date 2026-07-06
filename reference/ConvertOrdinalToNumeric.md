# Convert ordinal variables to numeric

Convert ordinal variables in a dataframe to numeric if they contain
numeric values in their character representation.

## Usage

``` r
ConvertOrdinalToNumeric(
  data,
  variables = NULL,
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataframe containing the variables.

- variables:

  A character vector specifying the names of variables to consider. If
  NULL, all columns of the dataframe will be considered.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

The dataframe with ordinal variables potentially converted to numeric.

## Examples

``` r
# An ordered factor with numeric levels, and one with non-numeric levels
df <- data.frame(
  id     = 1:5,
  likert = factor(c("1", "2", "3", "2", "1"),
                  levels = c("1", "2", "3"), ordered = TRUE),
  grade  = factor(c("A", "B", "A", "C", "B"),
                  levels = c("A", "B", "C"), ordered = TRUE)
)

out <- ConvertOrdinalToNumeric(df)

# likert becomes numeric; grade stays an ordered factor (levels are not numeric)
sapply(out, class)
#> $id
#> [1] "integer"
#> 
#> $likert
#> [1] "numeric"
#> 
#> $grade
#> [1] "ordered" "factor" 
#> 
```
