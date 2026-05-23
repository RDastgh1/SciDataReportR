# Prepare numeric data safely for analysis

Safely coerces selected variables in a data frame to numeric format
while preserving column names and replacing non-finite values (`Inf`,
`-Inf`, `NaN`) with `NA`.

## Usage

``` r
PrepNumericData(df, variables = names(df))
```

## Arguments

- df:

  A data frame.

- variables:

  Character vector of variable names to process. Defaults to all
  columns.

## Value

A data frame with selected variables converted to numeric and non-finite
values replaced with `NA`.

## Details

This function is useful for preparing real-world biomedical, omics, and
survey datasets for downstream statistical analyses such as
correlations, PCA, clustering, regression, and factor analysis.

Character and factor variables are converted using
`as.numeric(as.character(x))` to avoid incorrect factor level coercion.

## Examples

``` r
df <- data.frame(
  x = c("1", "2", "3"),
  y = c(1, 2, Inf),
  z = factor(c("4", "5", "6"))
)

PrepNumericData(df)
#>   x  y z
#> 1 1  1 4
#> 2 2  2 5
#> 3 3 NA 6
```
