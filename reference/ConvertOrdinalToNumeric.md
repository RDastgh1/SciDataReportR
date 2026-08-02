# Prepare ordinal variables for analysis

Applies a consistent ordinal-treatment policy to selected variables.
Ordinal score mappings recorded by
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
are used when available; otherwise ordered-factor ranks are used.

## Usage

``` r
ConvertOrdinalToNumeric(
  data,
  variables = NULL,
  TreatOrdinalAs = c("Continuous", "Categorical", "Both", "Exclude"),
  Relabel = TRUE,
  ReturnMetadata = FALSE,
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  The data frame containing the variables.

- variables:

  Character vector of variables to consider. If `NULL`, all columns are
  considered.

- TreatOrdinalAs:

  How ordinal variables are handled: `"Continuous"`, `"Categorical"`,
  `"Both"`, or `"Exclude"`.

- Relabel:

  Logical; when `TreatOrdinalAs = "Both"`, apply descriptive labels to
  the derived categorical and continuous variables.

- ReturnMetadata:

  Logical; if `FALSE` (default), return only the transformed data frame.
  If `TRUE`, return a list containing the data, selected variables,
  ordinal variables, variable map, and treatment.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A transformed data frame, or a metadata list when
`ReturnMetadata = TRUE`.

## Examples

``` r
df <- data.frame(
  id = 1:5,
  likert = factor(c("1", "2", "3", "2", "1"),
                  levels = c("1", "2", "3"), ordered = TRUE),
  grade = factor(c("A", "B", "A", "C", "B"),
                 levels = c("A", "B", "C"), ordered = TRUE)
)

# Convert ordinal ranks to numeric scores.
ConvertOrdinalToNumeric(df)
#>   id likert grade
#> 1  1      1     1
#> 2  2      2     2
#> 3  3      3     1
#> 4  4      2     3
#> 5  5      1     2

# Keep both categorical and continuous versions for a summary table.
ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Both", ReturnMetadata = TRUE)
#> $data
#>   id likert grade .scidr_ordinal_categorical_likert
#> 1  1      1     A                                 1
#> 2  2      2     B                                 2
#> 3  3      3     A                                 3
#> 4  4      2     C                                 2
#> 5  5      1     B                                 1
#>   .scidr_ordinal_continuous_likert .scidr_ordinal_categorical_grade
#> 1                                1                                A
#> 2                                2                                B
#> 3                                3                                A
#> 4                                2                                C
#> 5                                1                                B
#>   .scidr_ordinal_continuous_grade
#> 1                               1
#> 2                               2
#> 3                               1
#> 4                               3
#> 5                               2
#> 
#> $variables
#> [1] "id"                                ".scidr_ordinal_categorical_likert"
#> [3] ".scidr_ordinal_continuous_likert"  ".scidr_ordinal_categorical_grade" 
#> [5] ".scidr_ordinal_continuous_grade"  
#> 
#> $ordinal_variables
#> [1] "likert" "grade" 
#> 
#> $variable_map
#> $variable_map$id
#> [1] "id"
#> 
#> $variable_map$likert
#> [1] ".scidr_ordinal_categorical_likert" ".scidr_ordinal_continuous_likert" 
#> 
#> $variable_map$grade
#> [1] ".scidr_ordinal_categorical_grade" ".scidr_ordinal_continuous_grade" 
#> 
#> 
#> $treatment
#> [1] "Both"
#> 
```
