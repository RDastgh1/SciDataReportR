# Create a formatted data dictionary table

This function generates a formatted data dictionary table using the
specified data frame. The table includes variable names, labels, types,
and additional formatting based on variable types.

## Usage

``` r
FormattedDataDictionary(DataFrame, numdecimals = 2)
```

## Arguments

- DataFrame:

  The data frame for which the data dictionary is to be created.

- numdecimals:

  Number of decimals to display for numeric variables (default: 2).

## Value

A formatted data dictionary table (gt object).

## Details

This function requires the `gt` package. If not installed, the function
will return an error.
