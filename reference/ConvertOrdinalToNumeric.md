# ConvertOrdinalToNumeric

Convert ordinal variables in a dataframe to numeric if they contain
numeric values in their character representation.

## Usage

``` r
ConvertOrdinalToNumeric(Data, Variables = NULL)
```

## Arguments

- Data:

  The dataframe containing the variables.

- Variables:

  A character vector specifying the names of variables to consider. If
  NULL, all columns of the dataframe will be considered.

## Value

The dataframe with ordinal variables potentially converted to numeric.
