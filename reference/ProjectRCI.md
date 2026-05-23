# Project a trained RCI object onto new data

Apply previously learned regression-based Reliable Change Index (RCI)
models to a new or expanded longitudinal dataset WITHOUT refitting
models.

## Usage

``` r
ProjectRCI(Data, Object, ID = NULL)
```

## Arguments

- Data:

  A data frame.

- Object:

  A SciDataReportR_RCI object created with
  [`CreateRCIObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateRCIObject.md).

- ID:

  Optional ID column override.

## Value

A projected SciDataReportR_RCI object.

## Details

Supports both wide and long data structures and returns:

- tidy longitudinal RCI values

- original data merged with projected RCI outputs

- publication-ready plots
