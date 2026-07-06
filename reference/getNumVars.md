# Get Numeric Variables

Extracts numeric variables from a data frame.

## Usage

``` r
getNumVars(data, Ordinal = FALSE, DataFrame = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame from which to extract numeric variables.

- Ordinal:

  Logical, indicating whether to include ordinal variables.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing the names of numeric variables.
