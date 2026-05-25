# Add a Caption Explaining Star Annotations

This function adds a caption to a ggplot explaining the meaning of star
annotations (\*, \*\*, \*\*\*). It is most commonly added to correlation
heatmaps produced by
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
or to downstream plots derived from that function with
[`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md).

## Usage

``` r
geom_starcaption()
```

## Value

A [`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) object
that can be added to a ggplot, especially SciDataReportR heatmaps that
use star annotations.

## Input requirements

`geom_starcaption()` does not take an input object directly. Add it to a
ggplot with `+`, typically a plot returned inside the list created by
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
or a plot returned by
[`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md).
