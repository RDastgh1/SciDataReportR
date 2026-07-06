# Revalue Factors

This function revalues factor variables in a dataset according to the
specifications provided in a codebook.

## Usage

``` r
ReValueFactors(
  data,
  codebook,
  DatatoRevalue = lifecycle::deprecated(),
  VarTypes = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset to be revalued.

- codebook:

  A data frame containing information about the variables and how they
  should be revalued. It should have columns: Variable (variable names),
  Recode (yes/no for recoding), and Code (the revalue codes separated by
  "=" and ","). DO NOT USE COMMAS ANYWHERE ELSE IN THIS COLUMN.

- DatatoRevalue:

  **Deprecated** (since 19.15.0). Use `data` instead.

- VarTypes:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

The revalued dataset.
