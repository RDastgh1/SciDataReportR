# Create facet labels for ggplot2 based on variable labels in a data frame

This function takes a data frame containing variable labels and creates
facet labels suitable for use with ggplot2 facet functions.

## Usage

``` r
createFacetLabels(data, DataFrame = lifecycle::deprecated())
```

## Arguments

- data:

  A data frame containing variable labels.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing facet labels.
