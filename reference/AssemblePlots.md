# Assemble ggplot objects into a unified multi-panel figure

Combine ggplot objects or lists of ggplot objects into a consistently
styled multi-panel layout with optional legend collection, global
theming, automatic layout estimation, and suggested figure dimensions
for reporting workflows.

## Usage

``` r
AssemblePlots(
  Plots,
  ncol = NULL,
  nrow = NULL,
  AutoLayout = TRUE,
  RemoveNULL = TRUE,
  CollectLegend = TRUE,
  LegendPosition = "bottom",
  LegendRelativeSize = 0.08,
  Theme = ggplot2::theme_minimal(),
  BaseFontSize = 12,
  GlobalTheme = NULL,
  GlobalLayers = NULL,
  RemoveTitles = FALSE,
  UseNamesAsTitles = FALSE,
  Align = "hv",
  Axis = "tblr",
  Labels = NULL,
  LabelSize = 14,
  SuggestedBaseWidth = 4,
  SuggestedBaseHeight = 4,
  ReturnMetadata = FALSE
)
```

## Arguments

- Plots:

  A ggplot object or list of ggplot objects.

- ncol:

  Optional number of columns.

- nrow:

  Optional number of rows.

- AutoLayout:

  Logical; automatically estimate layout when `ncol` and `nrow` are not
  supplied.

- RemoveNULL:

  Logical; remove NULL plots before assembly.

- CollectLegend:

  Logical; combine legends into a shared legend.

- LegendPosition:

  Position of shared legend. One of `"top"`, `"bottom"`, `"left"`,
  `"right"`, or `"none"`.

- LegendRelativeSize:

  Relative size allocated to legend area.

- Theme:

  Optional global ggplot theme.

- BaseFontSize:

  Base font size applied globally.

- GlobalTheme:

  Optional additional theme applied globally.

- GlobalLayers:

  Optional list of ggplot layers/scales applied to all plots.

- RemoveTitles:

  Logical; remove plot titles globally.

- UseNamesAsTitles:

  Logical; use plot list names as titles when plot titles are missing.

- Align:

  Plot alignment passed to cowplot.

- Axis:

  Axis alignment passed to cowplot.

- Labels:

  Optional panel labels.

- LabelSize:

  Panel label font size.

- SuggestedBaseWidth:

  Base width per column used for suggested figure dimensions.

- SuggestedBaseHeight:

  Base height per row used for suggested figure dimensions.

- ReturnMetadata:

  Logical; if `TRUE`, returns plot plus metadata.

## Value

If `ReturnMetadata = FALSE`, returns a ggplot object.

If `ReturnMetadata = TRUE`, returns a list containing:

- Plot:

  Combined plot object

- nrow:

  Estimated number of rows

- ncol:

  Estimated number of columns

- SuggestedWidth:

  Suggested figure width

- SuggestedHeight:

  Suggested figure height

- NumPlots:

  Number of plots

## Details

This function is designed as a centralized visualization orchestration
layer for SciDataReportR workflows. It standardizes layout behavior
across analytical outputs while preserving full ggplot compatibility.

Internally, cowplot is used for layout assembly and legend extraction
because it provides more predictable spacing and legend behavior for
publication-style figures than patchwork.

Supports:

- Named or unnamed plot lists

- Automatic row/column estimation

- Global themes

- Global ggplot layers/scales

- Shared legends

- Automatic NULL removal

- Suggested figure dimensions

## Examples

``` r
library(ggplot2)

p1 <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point()

p2 <- ggplot(mtcars, aes(hp, wt)) +
  geom_point()

AssemblePlots(list(p1, p2))


AssemblePlots(
  list(MPG = p1, Horsepower = p2),
  UseNamesAsTitles = TRUE,
  LegendPosition = "top"
)

```
