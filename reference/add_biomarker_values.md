# Add values to a biomarker performance heatmap

Adds numeric cell labels to the heatmap returned by
[`ScreenBiomarkerPerformance()`](https://rdastgh1.github.io/SciDataReportR/reference/ScreenBiomarkerPerformance.md).
This helper is intended for downstream annotation of SciDataReportR
biomarker heatmaps while keeping the default heatmap uncluttered and
hover-friendly.

## Usage

``` r
add_biomarker_values(
  plot,
  value_var = "HeatmapValue",
  digits = 2,
  size = 3,
  color = "black"
)
```

## Arguments

- plot:

  A biomarker heatmap ggplot, typically
  `ScreenBiomarkerPerformance(...)$Plots$Heatmap`.

- value_var:

  Column in `plot$data` used for labels. Default is `"HeatmapValue"`.

- digits:

  Number of decimal places. Default is `2`.

- size:

  Text size. Default is `3`.

- color:

  Text color. Default is `"black"`.

## Value

A ggplot with numeric values added to each non-missing heatmap cell.

## Examples

``` r
if (FALSE) { # \dontrun{
screen$Plots$Heatmap %>% add_biomarker_values()
} # }
```
