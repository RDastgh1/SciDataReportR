# Create a data dictionary for a data frame

This function generates a data dictionary for a given data frame,
including variable names, labels, types, missing value statistics,
summary statistics, and other relevant information.

## Usage

``` r
Make_DataDictionary(DataFrame, numdecimals = 2)
```

## Arguments

- DataFrame:

  The data frame for which the data dictionary is to be created.

- numdecimals:

  Number of decimals to display for numeric variables (default: 2).

## Value

A data frame representing the data dictionary.
