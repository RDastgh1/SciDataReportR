# Reverse Levels of Categorical Factors

This function reverses the levels of specified categorical variables in
a given dataframe.

## Usage

``` r
reverseFactorLevels(data, variables, df = lifecycle::deprecated())
```

## Arguments

- data:

  A dataframe containing the categorical variables to be reversed.

- variables:

  A character vector of column names in the dataframe to reverse levels.
  These columns must be categorical factors.

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A dataframe with the levels of the specified factors reversed. Columns
not specified in `variables` remain unchanged.
