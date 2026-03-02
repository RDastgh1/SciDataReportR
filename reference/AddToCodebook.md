# Add a new variable to a codebook

This function adds a new variable entry to an existing codebook.

## Usage

``` r
AddToCodebook(
  CB,
  VariableName,
  VariableLabel = NA,
  VariableType = NA,
  VariableCategory = NA,
  VariableRecode = NA,
  VariableCode = NA,
  VariableExclude = NA,
  VariableNotes = NA
)
```

## Arguments

- CB:

  A data frame representing the codebook.

- VariableName:

  The name of the variable.

- VariableLabel:

  The label of the variable (default: same as VariableName).

- VariableType:

  The type of the variable (e.g., numeric, categorical).

- VariableCategory:

  The category to which the variable belongs.

- VariableRecode:

  The recode information for the variable.

- VariableCode:

  The code associated with the variable.

- VariableExclude:

  A flag indicating whether the variable should be excluded.

- VariableNotes:

  Any additional notes or comments about the variable.

## Value

A data frame representing the updated codebook with the new variable
added.

## Examples

``` r
# Create an empty codebook
codebook <- data.frame(Variable = character(0), Label = character(0),
                       Type = character(0), Category = character(0),
                       Recode = character(0), Code = character(0),
                       Exclude = logical(0), Notes = character(0))

# Add a new variable to the codebook
codebook <- AddToCodebook(codebook, "Age", "Age of participants", "numeric", "Demographics")
```
