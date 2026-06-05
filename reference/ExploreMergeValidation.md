# Explore merge validation results interactively

Create an interactive HTML dashboard from a
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
result object. This function is designed for merge quality-control
workflows. It displays dataset fingerprints, high-level summary cards, a
traffic-light checks table, coverage diagnostics, join-variable
auditing, and an expandable duplicate-variable conflict explorer.

## Usage

``` r
ExploreMergeValidation(
  MergeObj,
  Title = "Merge validation explorer",
  TopN = 10,
  TableHeight = 350
)
```

## Arguments

- MergeObj:

  A list returned by
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).

- Title:

  Character title shown at the top of the dashboard. Default is
  `"Merge validation explorer"`.

- TopN:

  Integer number of example variables or records to show in previews and
  expanded sections. Default is `10`.

- TableHeight:

  Height in pixels for scrollable reactable tables. Default is `350`.

## Value

An
[`htmltools::tagList()`](https://rstudio.github.io/htmltools/reference/tagList.html)
object containing an interactive dashboard.

## Details

This function is intended for interactive review rather than publication
tables. It returns an HTML object that can be rendered in the RStudio
Viewer, Quarto, R Markdown, Shiny, or saved as HTML. If needed in an
interactive console, wrap the result with
[`htmltools::browsable()`](https://rstudio.github.io/htmltools/reference/browsable.html).
