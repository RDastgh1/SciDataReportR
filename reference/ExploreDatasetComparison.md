# Explore dataset comparison results interactively

Create an interactive HTML dashboard from a
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
result object. This function is designed for data review and
quality-control workflows. It displays high-level summary cards, a
traffic-light checks table, and an expandable variable-change explorer
that shows side-by-side old and new values for modified cells.

## Usage

``` r
ExploreDatasetComparison(
  CompareObj,
  Title = "Dataset comparison explorer",
  TopN = 10
)
```

## Arguments

- CompareObj:

  A list returned by
  [`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md).

- Title:

  Character title shown at the top of the dashboard. Default is
  `"Dataset comparison explorer"`.

- TopN:

  Integer number of example variables or records to show in previews and
  expanded sections. Default is `10`.

## Value

An
[`htmltools::tagList()`](https://rstudio.github.io/htmltools/reference/tagList.html)
object containing an interactive dashboard.

## Details

This function is intended for interactive review rather than publication
tables. It returns an HTML object that can be rendered in the RStudio
Viewer, Quarto, R Markdown, or Shiny.
