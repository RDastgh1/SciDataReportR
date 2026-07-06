# Create a formatted data dictionary table

This function generates a formatted data dictionary table using the
specified data frame. The table includes variable names, labels, types,
and additional formatting based on variable types.

## Usage

``` r
FormattedDataDictionary(
  data,
  digits = 2,
  DataFrame = lifecycle::deprecated(),
  numdecimals = lifecycle::deprecated()
)
```

## Arguments

- data:

  The data frame for which the data dictionary is to be created.

- digits:

  Number of decimals to display for numeric variables (default: 2).

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- numdecimals:

  **Deprecated** (since 19.15.0). Use `digits` instead.

## Value

A formatted data dictionary table (gt object).

## Details

This function requires the `gt` package. If not installed, the function
will return an error.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach variable labels and factor levels first so the dictionary shows
# human-readable labels and correct variable types.
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Formatted data dictionary for a subset of variables
FormattedDataDictionary(
  Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin")]
)


  

Variable1,2
```
