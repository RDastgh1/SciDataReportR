# Revalue Data

Revalues variables in a dataset using a VarTypes codebook.

## Usage

``` r
RevalueData(DatatoRevalue, VarTypes, missingVal = -999, splitchar = ";")
```

## Arguments

- DatatoRevalue:

  A data.frame or tibble to be revalued.

- VarTypes:

  A data.frame with columns: Variable, Recode, Code, Type, Label,
  MissingCode. Only Variable is required. (Backward compatible: if
  MissingCode is absent/NA, will fall back to Missing.)

- missingVal:

  Default value to treat as missing when VarTypes\$MissingCode is absent
  or NA.

- splitchar:

  Separator used in VarTypes\$Code between pairs (default ";").

## Value

A list with: RevaluedData (data), warninglist (character), recodedvars
(character), not_in_data (character).
