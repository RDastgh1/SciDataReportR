# Revalue Data

Revalues variables in a dataset using a VarTypes codebook.

## Usage

``` r
RevalueData(
  data,
  codebook,
  missingVal = -999,
  splitchar = ";",
  on_error = c("stop", "warn"),
  DatatoRevalue = lifecycle::deprecated(),
  VarTypes = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data.frame or tibble to be revalued.

- codebook:

  A data.frame with columns: Variable, Recode, Code, Type, Label,
  MissingCode. Only Variable is required. (Backward compatible: if
  MissingCode is absent/NA, will fall back to Missing.)

- missingVal:

  Default value to treat as missing when VarTypes\$MissingCode is absent
  or NA.

- splitchar:

  Separator used in VarTypes\$Code between pairs (default ";").

- on_error:

  Whether to stop at the first variable-level error (the default) or
  continue and record errors in the returned object.

- DatatoRevalue:

  **Deprecated** (since 19.15.0). Use `data` instead.

- VarTypes:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

A list with: RevaluedData (data), warninglist (character), recodedvars
(character), not_in_data (character), and errors (data frame with
`Variable` and `Error` columns). In the default `on_error = "stop"`
mode, an error names the offending variable and preserves the underlying
message.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Before: no variable labels and coded factors are still raw
# (e.g. sex is stored as 0/1 with no label)
sjlabelled::get_label(SampleData$age)   # NULL
#> NULL
class(SampleData$sex)                    # "integer"
#> [1] "integer"

# Revalue using the codebook: attach labels and recode coded factors
revalued <- RevalueData(SampleData, SampleVariableTypes)
Labelled <- revalued$RevaluedData

# After: labels are attached and sex is a labelled factor
sjlabelled::get_label(Labelled$age)     # "Age"
#> [1] "Age"
levels(Labelled$sex)                     # "Female" "Male"
#> [1] "Female" "Male"  

# Variables that were recoded, and any codebook entries not found in the data
revalued$recodedvars
#> [1] "sex"
revalued$not_in_data
#> [1] "Group"   "Race"    "Marital" "Height"  "Income"  "Smokes"  "Smoker" 
```
