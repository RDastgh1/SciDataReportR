# Get Categorical Variables

Extracts categorical variables from a data frame.

## Usage

``` r
getCatVars(data, Ordinal = TRUE, DataFrame = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame from which to extract categorical variables.

- Ordinal:

  Logical, indicating whether to include ordinal variables.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A character vector containing the names of categorical variables.
