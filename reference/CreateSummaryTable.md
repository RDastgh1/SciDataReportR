# Create Summary Table

Generate a descriptive summary table for specified variables in a
dataset.

## Usage

``` r
CreateSummaryTable(
  data,
  variables = NULL,
  digits = 2,
  Relabel = TRUE,
  Ordinal = lifecycle::deprecated(),
  TreatOrdinalAs = "Categorical",
  ScrollBoxHeight = "700px",
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  numdecimals = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset containing the variables of interest.

- variables:

  A character vector specifying the variables for which summary
  statistics will be calculated.

- digits:

  Number of decimal places to round the summary statistics.

- Relabel:

  Logical, indicating whether to use variable labels as column headers.

- Ordinal:

  Deprecated logical compatibility option; use `TreatOrdinalAs` instead.

- TreatOrdinalAs:

  How ordinal variables are handled. This numeric descriptive table
  accepts `"Continuous"` or `"Exclude"`.

- ScrollBoxHeight:

  Height of the scroll box for displaying the table.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- numdecimals:

  **Deprecated** (since 19.15.0). Use `digits` instead.

## Value

A formatted HTML table displaying summary statistics.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Distribution diagnostics for six biomarkers
summary_table <- CreateSummaryTable(
  data = df_Labelled,
  variables = c("AXL", "Adiponectin", "Ferritin", "MMP7", "tau", "p_tau"),
  digits = 2,
  ScrollBoxHeight = "320px"
)
#> Warning: no DISPLAY variable so Tk is not available

# `browsable()` is what renders the HTML on this page; printing suffices
# in Quarto or R Markdown
htmltools::browsable(htmltools::HTML(as.character(summary_table)))

Descriptive Summary Table. IQR, Skewness, and Kurtosis are highlighted in yellow if they are indicative of a non-normal distribution. Pct.Valid is highlighted in red if over 30% of data is missing
 
   
```
