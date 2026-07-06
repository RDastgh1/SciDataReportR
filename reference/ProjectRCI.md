# Project a trained RCI object onto new data

Apply previously learned regression-based Reliable Change Index (RCI)
models to a new or expanded longitudinal dataset WITHOUT refitting
models.

## Usage

``` r
ProjectRCI(
  data,
  Object,
  id_var = NULL,
  Data = lifecycle::deprecated(),
  ID = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame.

- Object:

  A SciDataReportR_RCI object created with
  [`CreateRCIObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateRCIObject.md).

- id_var:

  Optional ID column override.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- ID:

  **Deprecated** (since 19.15.0). Use `id_var` instead.

## Value

A projected SciDataReportR_RCI object.

## Details

Supports both wide and long data structures and returns:

- tidy longitudinal RCI values

- original data merged with projected RCI outputs

- publication-ready plots

## Examples

``` r
if (FALSE) { # \dontrun{
# NOTE: This example is not run because ProjectRCI() currently errors with
# "object 'VariableTable' not found" - VariableTable is used but never
# defined in the function body (tracked bug, to be fixed separately).
rci_data <- data.frame(
  id = rep(1:30, each = 2),
  visit = rep(c("Baseline", "Followup"), 30),
  Score = round(rnorm(60, mean = 50, sd = 10), 1)
)

rci <- CreateRCIObject(
  data = rci_data,
  variables = "Score",
  DataFormat = "long",
  id_var = "id",
  VisitColumn = "visit",
  BaselineVisit = "Baseline"
)

projected <- ProjectRCI(data = rci_data, Object = rci, id_var = "id")
projected$Plots$Spaghetti$Score
} # }
```
