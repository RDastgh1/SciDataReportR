# Plot a spider chart across continuous and binary variables

This function summarizes a set of variables and displays them on a
spider chart. Continuous variables are plotted as mean z-scores by
default using
[`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md),
while binary variables are plotted as percentages. It can overlay groups
on one spider chart or facet by group, optionally fill the polygons,
relabel spokes using variable labels, wrap long labels, reorder
variables to visually emphasize between-group differences, and
optionally return an interactive radar chart using plotly.

## Usage

``` r
PlotSpiderChart(
  data,
  variables,
  group_var = NULL,
  Relabel = TRUE,
  ContinuousSummary = "mean",
  ContinuousScaling = "zscore",
  Fill = FALSE,
  FillAlpha = 0.2,
  Facet = FALSE,
  VariableOrder = "input",
  VariableCategories = NULL,
  BinaryPositiveValue = 1,
  Palette = "Dark 3",
  LineSize = 1,
  PointSize = 2,
  ShowPoints = FALSE,
  LegendTitle = NULL,
  PlotTitle = NULL,
  Subtitle = NULL,
  Caption = NULL,
  AxisLabelSize = 12,
  AxisTextSize = 10,
  StripTextSize = 11,
  WrapLabels = TRUE,
  LabelWrapWidth = 22,
  LabelRadiusMultiplier = 1.22,
  PlotMarginTop = 40,
  PlotMarginRight = 120,
  PlotMarginBottom = 40,
  PlotMarginLeft = 120,
  interactive = FALSE,
  InteractiveHeight = 700,
  InteractiveWidth = NULL,
  InteractiveAxisMin = NULL,
  InteractiveAxisMax = NULL,
  tooltip_digits = 2,
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  GroupVariable = lifecycle::deprecated(),
  TooltipDigits = lifecycle::deprecated(),
  MakeInteractive = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame.

- variables:

  Character vector of variable names to plot.

- group_var:

  Optional grouping variable name. If NULL, one overall summary profile
  is plotted.

- Relabel:

  Logical; if TRUE, use variable labels when available.

- ContinuousSummary:

  Character; one of `"mean"` or `"median"`.

- ContinuousScaling:

  Character; one of `"zscore"`, `"none"`, or `"minmax"`.

- Fill:

  Logical; if TRUE, add transparent polygon fills in the static ggplot
  version and filled polygons in the interactive version.

- FillAlpha:

  Numeric transparency for fills.

- Facet:

  Logical; if TRUE and `GroupVariable` is supplied, facet by group
  instead of overlaying all groups on one spider chart. Ignored when
  `interactive = TRUE`.

- VariableOrder:

  Character; one of `"input"`, `"discrimination"`, `"hierarchical"`,
  `"greedy"`, or `"category_discrimination"`.

- VariableCategories:

  Optional character vector of categories for `Variables`. Must be the
  same length as `Variables` when supplied.

- BinaryPositiveValue:

  Optional positive value to use for non-factor binary variables.
  Defaults to `1`. For factor variables, the second factor level is
  used.

- Palette:

  Character name of an
  [`hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html) palette.
  Defaults to `"Dark 3"`.

- LineSize:

  Numeric line width for the static ggplot version.

- PointSize:

  Numeric point size for the static ggplot version.

- ShowPoints:

  Logical; if TRUE, show points at each spoke in the static ggplot
  version.

- LegendTitle:

  Optional legend title. Defaults to `GroupVariable`.

- PlotTitle:

  Optional plot title.

- Subtitle:

  Optional plot subtitle.

- Caption:

  Optional plot caption.

- AxisLabelSize:

  Numeric axis text size for spoke labels in the static ggplot version.

- AxisTextSize:

  Numeric text size for radial axis labels in the static ggplot version.

- StripTextSize:

  Numeric facet strip text size in the static ggplot version.

- WrapLabels:

  Logical; if TRUE, wrap long spoke labels.

- LabelWrapWidth:

  Numeric wrap width passed to
  [`stringr::str_wrap()`](https://stringr.tidyverse.org/reference/str_wrap.html).

- LabelRadiusMultiplier:

  Numeric multiplier controlling how far labels sit outside the spider
  in the static ggplot version.

- PlotMarginTop:

  Numeric top plot margin for the static ggplot version.

- PlotMarginRight:

  Numeric right plot margin for the static ggplot version.

- PlotMarginBottom:

  Numeric bottom plot margin for the static ggplot version.

- PlotMarginLeft:

  Numeric left plot margin for the static ggplot version.

- interactive:

  Logical; if TRUE, return an interactive plotly radar chart instead of
  a static ggplot. Default is `FALSE`.

- InteractiveHeight:

  Numeric height in pixels for the interactive widget.

- InteractiveWidth:

  Optional width passed to plotly layout. Defaults to NULL.

- InteractiveAxisMin:

  Optional numeric minimum for the interactive radial axis. If NULL,
  auto-detected from the summarized values.

- InteractiveAxisMax:

  Optional numeric maximum for the interactive radial axis. If NULL,
  auto-detected from the summarized values.

- tooltip_digits:

  Integer number of digits to show in interactive tooltips.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- GroupVariable:

  **Deprecated** (since 19.15.0). Use `group_var` instead.

- TooltipDigits:

  **Deprecated** (since 19.15.0). Use `tooltip_digits` instead.

- MakeInteractive:

  **Deprecated** (since 19.15.0). Use `interactive` instead.

## Value

A ggplot object when `interactive = FALSE`, otherwise a plotly
htmlwidget.

## Examples

``` r
data(SampleData)

# Static spider chart
PlotSpiderChart(
  SampleData,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  group_var = "Diagnosis"
)


# Interactive (plotly) version of the same chart
PlotSpiderChart(
  SampleData,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  group_var = "Diagnosis",
  interactive = TRUE
)
#> Warning: Specifying width/height in layout() is now deprecated.
#> Please specify in ggplotly() or plot_ly()

{"x":{"visdat":{"225b2bd85b64":["function () ","plotlyVisDat"]},"cur_data":"225b2bd85b64","attrs":{"225b2bd85b64":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatterpolar","mode":"lines","r":[0.047735667776953747,-0.037764082196990657,-0.066786105495183712,-0.12628285475807199,0.047735667776953747],"theta":["age","AXL","Adiponectin","Alpha_1_Antitrypsin","age"],"name":"Control","line":{"color":"#E16A86","width":2},"marker":{"color":"#E16A86","size":6},"fill":"none","fillcolor":"#E16A8633","text":["Diagnosis: Control<br>Variable: age<br>Raw name: age<br>Value:  0.05","Diagnosis: Control<br>Variable: AXL<br>Raw name: AXL<br>Value: -0.04","Diagnosis: Control<br>Variable: Adiponectin<br>Raw name: Adiponectin<br>Value: -0.07","Diagnosis: Control<br>Variable: Alpha_1_Antitrypsin<br>Raw name: Alpha_1_Antitrypsin<br>Value: -0.13","Diagnosis: Control<br>Variable: age<br>Raw name: age<br>Value:  0.05"],"hoverinfo":"text","inherit":true},"225b2bd85b64.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatterpolar","mode":"lines","r":[-0.12694540222003084,0.10042755924914015,0.17760700582235767,0.33582913023575273,-0.12694540222003084],"theta":["age","AXL","Adiponectin","Alpha_1_Antitrypsin","age"],"name":"Impaired","line":{"color":"#00AD9A","width":2},"marker":{"color":"#00AD9A","size":6},"fill":"none","fillcolor":"#00AD9A33","text":["Diagnosis: Impaired<br>Variable: age<br>Raw name: age<br>Value: -0.13","Diagnosis: Impaired<br>Variable: AXL<br>Raw name: AXL<br>Value:  0.10","Diagnosis: Impaired<br>Variable: Adiponectin<br>Raw name: Adiponectin<br>Value:  0.18","Diagnosis: Impaired<br>Variable: Alpha_1_Antitrypsin<br>Raw name: Alpha_1_Antitrypsin<br>Value:  0.34","Diagnosis: Impaired<br>Variable: age<br>Raw name: age<br>Value: -0.13"],"hoverinfo":"text","inherit":true}},"layout":{"height":700,"margin":{"b":40,"l":60,"t":25,"r":10},"title":{"text":""},"polar":{"radialaxis":{"visible":true,"range":[-1,1.1000000000000001]},"angularaxis":{"direction":"clockwise","rotation":90,"categoryorder":"array","categoryarray":["age","AXL","Adiponectin","Alpha_1_Antitrypsin"],"tickfont":{"size":12}}},"legend":{"title":{"text":"Diagnosis"},"orientation":"v"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"type":"scatterpolar","mode":"lines+markers","r":[0.047735667776953747,-0.037764082196990657,-0.066786105495183712,-0.12628285475807199,0.047735667776953747],"theta":["age","AXL","Adiponectin","Alpha_1_Antitrypsin","age"],"name":"Control","line":{"color":"#E16A86","width":2},"marker":{"color":"#E16A86","size":6,"line":{"color":"rgba(31,119,180,1)"}},"fill":"none","fillcolor":"#E16A8633","text":["Diagnosis: Control<br>Variable: age<br>Raw name: age<br>Value:  0.05","Diagnosis: Control<br>Variable: AXL<br>Raw name: AXL<br>Value: -0.04","Diagnosis: Control<br>Variable: Adiponectin<br>Raw name: Adiponectin<br>Value: -0.07","Diagnosis: Control<br>Variable: Alpha_1_Antitrypsin<br>Raw name: Alpha_1_Antitrypsin<br>Value: -0.13","Diagnosis: Control<br>Variable: age<br>Raw name: age<br>Value:  0.05"],"hoverinfo":["text","text","text","text","text"],"frame":null},{"type":"scatterpolar","mode":"lines+markers","r":[-0.12694540222003084,0.10042755924914015,0.17760700582235767,0.33582913023575273,-0.12694540222003084],"theta":["age","AXL","Adiponectin","Alpha_1_Antitrypsin","age"],"name":"Impaired","line":{"color":"#00AD9A","width":2},"marker":{"color":"#00AD9A","size":6,"line":{"color":"rgba(255,127,14,1)"}},"fill":"none","fillcolor":"#00AD9A33","text":["Diagnosis: Impaired<br>Variable: age<br>Raw name: age<br>Value: -0.13","Diagnosis: Impaired<br>Variable: AXL<br>Raw name: AXL<br>Value:  0.10","Diagnosis: Impaired<br>Variable: Adiponectin<br>Raw name: Adiponectin<br>Value:  0.18","Diagnosis: Impaired<br>Variable: Alpha_1_Antitrypsin<br>Raw name: Alpha_1_Antitrypsin<br>Value:  0.34","Diagnosis: Impaired<br>Variable: age<br>Raw name: age<br>Value: -0.13"],"hoverinfo":["text","text","text","text","text"],"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
