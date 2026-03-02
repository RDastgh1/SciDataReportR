# Identify Binary Variables

This function identifies and returns a list of binary variables in a
dataframe. Binary variables are defined as having exactly two unique
values or levels. The function supports options for handling ordinal
factors and revalued data.

## Usage

``` r
getBinaryVars(DataFrame, Ordinal = TRUE, Revalued = TRUE)
```

## Arguments

- DataFrame:

  A dataframe to analyze for binary variables.

- Ordinal:

  Logical. If TRUE, ordinal factors are included in the search for
  binary variables. Default is TRUE.

- Revalued:

  Logical. If TRUE, the function checks factors and their levels;
  otherwise, it checks for variables with two unique values.

## Value

A character vector containing the names of binary variables.
