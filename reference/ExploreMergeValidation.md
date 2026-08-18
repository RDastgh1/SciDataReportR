# Explore merge validation results interactively

Create an interactive HTML dashboard from a
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
result object. This function is designed for merge quality-control
workflows. It displays a traffic-light checks table (with
rows/columns/unique-key context folded in as informational rows), a
coverage explorer, an expandable duplicate-variable conflict explorer,
and suggested actions.

## Usage

``` r
ExploreMergeValidation(
  MergeObj,
  Title = "Merge validation explorer",
  TopN = 10,
  TableHeight = 350,
  Detail = c("Compact", "Full")
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

- Detail:

  Either `"Compact"` (default) or `"Full"`. In `"Compact"` mode, the
  coverage explorer and conflicts explorer render as collapsed
  click-to-expand accordion sections labeled with their item counts. In
  `"Full"` mode, the same sections are expanded by default. In both
  modes, sections with nothing to show (no unmatched keys, no duplicated
  variables, no suggested actions) are omitted entirely.

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
The reference page uses a static preview because pkgdown cannot host the
live dashboard.

## Dashboard preview

![Merge validation dashboard showing merge checks, coverage review, and
recommended actions](figures/ExploreMergeValidation.png)

## Examples

``` r
set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Render the interactive dashboard in an HTML report, Shiny app, Viewer,
# or standalone HTML file.
dashboard <- ExploreMergeValidation(validation)
```
