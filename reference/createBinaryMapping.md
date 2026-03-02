# Create a Mapping Table for Binary Variables

Identifies binary variables and returns a deterministic mapping for 0/1
coding.

- Factors with *explicit order* (ordered = TRUE) do NOT use heuristics;
  the highest (last) level is Positive.

- Logicals map to Negative = "FALSE", Positive = "TRUE".

- Numeric 0/1 (or any 2-value numeric) maps Positive to the numeric
  maximum.

- Characters / unordered factors use minimal heuristics (no race/PWH/sex
  terms).

## Usage

``` r
createBinaryMapping(Data, CatVars, prefer = NULL)
```

## Arguments

- Data:

  A dataframe.

- CatVars:

  Character vector of candidate binary variables.

- prefer:

  Optional named character vector of explicit positive levels, e.g.,
  c(STATUS = "PWH", Smoker = "Yes"). This overrides other rules.

## Value

A data.frame with columns: Variable, Label, PositiveLevel,
NegativeLevel.
