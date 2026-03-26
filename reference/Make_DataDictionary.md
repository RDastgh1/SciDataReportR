# Create a data dictionary for a data frame

This function generates a data dictionary for a given data frame,
including variable names, labels, types, missing value statistics,
summary statistics, and other relevant information. Use this when you
want a stable, label-aware summary table across mixed variable types.

## Usage

``` r
Make_DataDictionary(DataFrame, numdecimals = 2)
```

## Arguments

- DataFrame:

  A data frame.

- numdecimals:

  Number of decimals to display for numeric variables.

## Value

A tibble with one row per variable and stable summary columns.

## Examples

``` r
df <- tibble::tibble(
  age = c(20, 25, 30),
  group = factor(c("A", "B", "A"))
)
Make_DataDictionary(df)
#> # A tibble: 2 × 14
#>   Variable Label Type    `Ordered Factor` n_missing complete_rate n_unique
#>   <chr>    <chr> <chr>   <chr>            <chr>     <chr>         <chr>   
#> 1 age      ""    numeric " "              0         1             3       
#> 2 group    ""    factor  "FALSE"          0         1             2       
#> # ℹ 7 more variables: `Top Counts` <chr>, Mean <chr>, Median <chr>, SD <chr>,
#> #   Min <chr>, Max <chr>, Histogram <chr>
```
