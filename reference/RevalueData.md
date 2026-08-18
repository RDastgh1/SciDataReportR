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

## What changes

In the raw extract `sex` is a bare 0/1 column and nothing carries a
label. Afterwards it is a factor with real levels, and the labels are
attached for every downstream table and plot to pick up automatically -
which is what makes output readable without renaming anything by hand.

Anything the codebook could not act on is reported in `warninglist`
rather than silently skipped, so a mistyped variable name or an
unparseable `Code` string is visible instead of quietly leaving a
variable unrecoded.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Before: no labels, and sex is stored as 0/1
sjlabelled::get_label(SampleData$age)   # NULL
#> NULL
class(SampleData$sex)                    # "integer"
#> [1] "integer"

# Revalue using the codebook
revalued <- RevalueData(SampleData, SampleVariableTypes)
Labelled <- revalued$RevaluedData

# After: labels attached, sex is a labelled factor
sjlabelled::get_label(Labelled$age)     # "Age"
#> [1] "Age"
levels(Labelled$sex)                     # "Female" "Male"
#> [1] "Female" "Male"  

# Recoded variables, and codebook entries not found in the data
revalued$recodedvars
#> [1] "sex"
revalued$not_in_data
#> [1] "Group"   "Race"    "Marital" "Height"  "Income"  "Smokes"  "Smoker" 

# \donttest{
# A side-by-side view of what changed
vars_Show <- c("Diagnosis", "age", "sex", "Genotype", "AXL")

df_Before <- utils::head(SampleData[, vars_Show], 6)
df_After <- utils::head(Labelled[, vars_Show], 6)

ShowTable <- function(x, caption) {
  htmltools::browsable(htmltools::HTML(as.character(
    kableExtra::kable_styling(
      knitr::kable(x, format = "html", caption = caption, row.names = FALSE),
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
  )))
}

ShowTable(df_Before, "Before: raw codes as imported")

Before: raw codes as imported
 
 Diagnosis 
```
