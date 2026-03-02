# Reverse Levels of Categorical Factors

This function reverses the levels of specified categorical variables in
a given dataframe.

## Usage

``` r
reverseFactorLevels(df, variables)
```

## Arguments

- df:

  A dataframe containing the categorical variables to be reversed.

- variables:

  A character vector of column names in the dataframe to reverse levels.
  These columns must be categorical factors.

## Value

A dataframe with the levels of the specified factors reversed. Columns
not specified in `variables` remain unchanged.
