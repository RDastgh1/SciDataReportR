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

## The four policies

`"Continuous"` replaces each level with its rank. `"Categorical"` and
`"Exclude"` both leave the values alone; they differ in whether the
variable is offered to the analysis at all. `"Both"` produces one column
of each.

With `ReturnMetadata = TRUE` the derived column names come back
alongside the data, so downstream functions know which column carries
which treatment without having to guess from the naming convention.

## Examples

``` r
# \donttest{
df <- data.frame(
  id = 1:6,
  Severity = factor(
    c("None", "Mild", "Severe", "Mild", "Moderate", "None"),
    levels = c("None", "Mild", "Moderate", "Severe"), ordered = TRUE
  ),
  Education = factor(
    c("HighSchool", "College", "Graduate", "College", "Graduate", "HighSchool"),
    levels = c("HighSchool", "College", "Graduate"), ordered = TRUE
  )
)

# The four policies, applied to the same ordinal variable
df_Both <- ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Both")

df_Policies <- data.frame(
  id = df$id,
  Original = as.character(df$Severity),
  Continuous = ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Continuous")$Severity,
  Categorical = as.character(
    ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Categorical")$Severity
  ),
  Both_Categorical = as.character(df_Both$.scidr_ordinal_categorical_Severity),
  Both_Continuous = df_Both$.scidr_ordinal_continuous_Severity
)

htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(df_Policies, full_width = TRUE)
)))


 id 
```
