# Create a Forest Plot from Univariate Regression Tables

This function generates a forest plot from a list of formatted
univariate regression tables.

## Usage

``` r
plotForestFromTable(UnivariateRegressionTables, pSize = 2)
```

## Arguments

- UnivariateRegressionTables:

  A list containing formatted regression tables and styling information.
  Expected structure:

  - `FormattedTable$tbls`: A list of tables, each containing a
    `table_body` dataframe.

  - `LargeTable$table_styling$header`: A dataframe with a `label` and
    `spanning_header` column for headers.

- pSize:

  Numeric. Size of the points in the plot. Default is 2.

## Value

A ggplot object representing the forest plot.

## Examples

``` r
# Example usage:
# plotForestFromTable(UnivariateRegressionTables)
```
