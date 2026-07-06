# Update an existing codebook based on a given dataframe

This function updates a codebook by adding missing variables, removing
variables that no longer exist in the dataframe (if specified), and
optionally replacing outdated labels.

## Usage

``` r
UpdateCodebook(
  data,
  codebook,
  RemoveMissing = TRUE,
  ReplaceLabels = FALSE,
  Dataframe = lifecycle::deprecated(),
  Codebook = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe for which the codebook needs to be updated.

- codebook:

  A dataframe representing the existing codebook with at least a
  'Variable' column.

- RemoveMissing:

  Logical; if TRUE, removes variables from the codebook that are not in
  the dataframe.

- ReplaceLabels:

  Logical; if TRUE, replaces outdated labels in the codebook with new
  ones from the dataframe.

- Dataframe:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Codebook:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

A list containing:

- `UpdatedCodebook`: The updated codebook dataframe.

- `NewVariables`: Variables present in the dataframe but missing from
  the original codebook.

- `NotExistingVariables`: Variables in the codebook that are not present
  in the dataframe.

- `MismatchedLabels`: A dataframe of variables with mismatched labels
  between the codebook and the dataframe.
