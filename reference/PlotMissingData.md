# Plot Missing Data

Visualize missing data patterns with variables as rows and observations
as columns. Optional hover variables can be included to facilitate
quality control workflows when converting the plot to an interactive
Plotly figure.

## Usage

``` r
PlotMissingData(
  data,
  variables = NULL,
  HoverVars = NULL,
  Relabel = TRUE,
  show_perc = TRUE,
  show_perc_var = TRUE,
  cluster = FALSE,
  DataFrame = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame.

- variables:

  Character vector of variables to visualize. If NULL, all columns
  except HoverVars are used.

- HoverVars:

  Optional character vector of columns to include in hover text. Useful
  for participant IDs, visit names, dates, sites, etc.

- Relabel:

  Logical. If TRUE, variable labels are used when available.

- show_perc:

  Logical. If TRUE, overall missingness percentages are shown in the
  legend.

- show_perc_var:

  Logical. If TRUE, variable-specific missingness percentages are
  appended to y-axis labels.

- cluster:

  Logical. If TRUE, variables are clustered by missingness pattern.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A ggplot object.

## Examples

``` r
data(SampleData)

# SampleData has real missingness in several assays
vars <- c("age", "AXL", "Angiotensinogen", "BMP_6", "IL_6",
          "Fetuin_A", "NT_proBNP", "ENA_78")

PlotMissingData(
  SampleData,
  variables = vars,
  HoverVars = "Diagnosis"
)

```
