# Revalue Factors

This function revalues factor variables in a dataset according to the
specifications provided in a codebook.

## Usage

``` r
ReValueFactors(DatatoRevalue, VarTypes)
```

## Arguments

- DatatoRevalue:

  The dataset to be revalued.

- VarTypes:

  A data frame containing information about the variables and how they
  should be revalued. It should have columns: Variable (variable names),
  Recode (yes/no for recoding), and Code (the revalue codes separated by
  "=" and ","). DO NOT USE COMMAS ANYWHERE ELSE IN THIS COLUMN.

## Value

The revalued dataset.
