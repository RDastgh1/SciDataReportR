# Create a data dictionary for a data frame

This function generates a stable data dictionary for a data frame using
[`codebook::skim_codebook()`](https://rubenarslan.github.io/codebook/reference/skim_codebook.html).
It is designed to work even when skim output omits type-specific summary
columns, such as numeric summaries for data frames without numeric
variables.

`Make_DataDictionary()` has been superseded by `MakeDataDictionary()`.
It remains available as a backwards-compatible alias and returns the
same data dictionary.

## Usage

``` r
MakeDataDictionary(
  data,
  digits = 2,
  DataFrame = lifecycle::deprecated(),
  numdecimals = lifecycle::deprecated()
)

Make_DataDictionary(...)
```

## Arguments

- data:

  A data frame.

- digits:

  Number of decimals to display for numeric variables.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- numdecimals:

  **Deprecated** (since 19.15.0). Use `digits` instead.

- ...:

  Arguments passed to `MakeDataDictionary()`.

## Value

A data frame with one row per variable and stable summary columns.

## Examples

``` r
df <- tibble::tibble(
  group = factor(c("A", "B", "A")),
  status = c("yes", "no", "yes")
)
if (requireNamespace("codebook", quietly = TRUE)) {
  MakeDataDictionary(df)
}
#> # A tibble: 2 × 14
#>   Variable Label Type      `Ordered Factor` n_missing complete_rate n_unique
#>   <chr>    <chr> <chr>     <chr>            <chr>     <chr>         <chr>   
#> 1 group    ""    factor    "FALSE"          0         1             2       
#> 2 status   ""    character " "              0         1             2       
#> # ℹ 7 more variables: `Top Counts` <chr>, Mean <chr>, Median <chr>, SD <chr>,
#> #   Min <chr>, Max <chr>, Histogram <chr>
```
