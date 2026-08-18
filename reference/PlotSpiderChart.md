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
  Palette = NULL,
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

  Optional character name of an
  [`hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html) palette.
  When `NULL` (the default), the SciDataReportR palette is used. Passing
  a name such as `"Dark 3"` still works exactly as before.

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
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

vars_biomarkers <- c(
  "Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10", "MMP7", "TRAIL_R3"
)
categories_biomarkers <- c(
  "Neurodegeneration", "Neurodegeneration", "Neurodegeneration",
  "Inflammation", "Inflammation", "Matrix remodeling", "TNF signaling"
)

# Input order preserves the supplied clinical/domain sequence.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "input"
)


# Discrimination puts the largest between-group differences next to each other.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "discrimination"
)


# Hierarchical order groups biomarkers with similar Diagnosis profiles.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "hierarchical"
)


# Greedy order places consecutive spokes with maximally different profiles.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "greedy"
)


# Category/discrimination order keeps domains together, then ranks the
# biomarkers within each domain by their between-group difference.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "category_discrimination",
  VariableCategories = categories_biomarkers
)


# Interactive (plotly) version of the category-aware chart.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "category_discrimination",
  VariableCategories = categories_biomarkers,
  interactive = TRUE
)

{"x":{"visdat":{"224ad7e4749":["function () ","plotlyVisDat"]},"cur_data":"224ad7e4749","attrs":{"224ad7e4749":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatterpolar","mode":"lines","r":[0.28368409270840433,-0.23869554746228758,-0.21267843354168886,-0.20034930721817537,-0.19239788626634546,-0.15054765460683786,-0.1820385545890647,0.28368409270840433],"theta":["Amyloid beta 42","Tau protein","Phosphorylated tau<br />protein","Growth-regulated alpha<br />protein","Matrix<br />metalloproteinase 10","Matrix<br />metalloproteinase 7","TNF-related<br />apoptosis-inducing<br />ligand receptor 3","Amyloid beta 42"],"name":"Control","line":{"color":"#0B1F5E","width":2},"marker":{"color":"#0B1F5E","size":6},"fill":"none","fillcolor":"#0B1F5E33","text":["Diagnosis: Control<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value:  0.28","Diagnosis: Control<br>Variable: Tau protein<br>Raw name: tau<br>Value: -0.24","Diagnosis: Control<br>Variable: Phosphorylated tau<br />protein<br>Raw name: p_tau<br>Value: -0.21","Diagnosis: Control<br>Variable: Growth-regulated alpha<br />protein<br>Raw name: GRO_alpha<br>Value: -0.20","Diagnosis: Control<br>Variable: Matrix<br />metalloproteinase 10<br>Raw name: MMP10<br>Value: -0.19","Diagnosis: Control<br>Variable: Matrix<br />metalloproteinase 7<br>Raw name: MMP7<br>Value: -0.15","Diagnosis: Control<br>Variable: TNF-related<br />apoptosis-inducing<br />ligand receptor 3<br>Raw name: TRAIL_R3<br>Value: -0.18","Diagnosis: Control<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value:  0.28"],"hoverinfo":"text","inherit":true},"224ad7e4749.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatterpolar","mode":"lines","r":[-0.70285909536709201,0.70394958065149327,0.56558440568229151,0.53279705875603434,0.51165152171929296,0.40035749906433782,0.48410252978630425,-0.70285909536709201],"theta":["Amyloid beta 42","Tau protein","Phosphorylated tau<br />protein","Growth-regulated alpha<br />protein","Matrix<br />metalloproteinase 10","Matrix<br />metalloproteinase 7","TNF-related<br />apoptosis-inducing<br />ligand receptor 3","Amyloid beta 42"],"name":"Impaired","line":{"color":"#E87400","width":2},"marker":{"color":"#E87400","size":6},"fill":"none","fillcolor":"#E8740033","text":["Diagnosis: Impaired<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value: -0.70","Diagnosis: Impaired<br>Variable: Tau protein<br>Raw name: tau<br>Value:  0.70","Diagnosis: Impaired<br>Variable: Phosphorylated tau<br />protein<br>Raw name: p_tau<br>Value:  0.57","Diagnosis: Impaired<br>Variable: Growth-regulated alpha<br />protein<br>Raw name: GRO_alpha<br>Value:  0.53","Diagnosis: Impaired<br>Variable: Matrix<br />metalloproteinase 10<br>Raw name: MMP10<br>Value:  0.51","Diagnosis: Impaired<br>Variable: Matrix<br />metalloproteinase 7<br>Raw name: MMP7<br>Value:  0.40","Diagnosis: Impaired<br>Variable: TNF-related<br />apoptosis-inducing<br />ligand receptor 3<br>Raw name: TRAIL_R3<br>Value:  0.48","Diagnosis: Impaired<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value: -0.70"],"hoverinfo":"text","inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":{"text":""},"polar":{"radialaxis":{"visible":true,"range":[-1,1.1000000000000001]},"angularaxis":{"direction":"clockwise","rotation":90,"categoryorder":"array","categoryarray":["Amyloid beta 42","Tau protein","Phosphorylated tau<br />protein","Growth-regulated alpha<br />protein","Matrix<br />metalloproteinase 10","Matrix<br />metalloproteinase 7","TNF-related<br />apoptosis-inducing<br />ligand receptor 3"],"tickfont":{"size":12}}},"legend":{"title":{"text":"Diagnosis"},"orientation":"v"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"type":"scatterpolar","mode":"lines+markers","r":[0.28368409270840433,-0.23869554746228758,-0.21267843354168886,-0.20034930721817537,-0.19239788626634546,-0.15054765460683786,-0.1820385545890647,0.28368409270840433],"theta":["Amyloid beta 42","Tau protein","Phosphorylated tau<br />protein","Growth-regulated alpha<br />protein","Matrix<br />metalloproteinase 10","Matrix<br />metalloproteinase 7","TNF-related<br />apoptosis-inducing<br />ligand receptor 3","Amyloid beta 42"],"name":"Control","line":{"color":"#0B1F5E","width":2},"marker":{"color":"#0B1F5E","size":6,"line":{"color":"rgba(31,119,180,1)"}},"fill":"none","fillcolor":"#0B1F5E33","text":["Diagnosis: Control<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value:  0.28","Diagnosis: Control<br>Variable: Tau protein<br>Raw name: tau<br>Value: -0.24","Diagnosis: Control<br>Variable: Phosphorylated tau<br />protein<br>Raw name: p_tau<br>Value: -0.21","Diagnosis: Control<br>Variable: Growth-regulated alpha<br />protein<br>Raw name: GRO_alpha<br>Value: -0.20","Diagnosis: Control<br>Variable: Matrix<br />metalloproteinase 10<br>Raw name: MMP10<br>Value: -0.19","Diagnosis: Control<br>Variable: Matrix<br />metalloproteinase 7<br>Raw name: MMP7<br>Value: -0.15","Diagnosis: Control<br>Variable: TNF-related<br />apoptosis-inducing<br />ligand receptor 3<br>Raw name: TRAIL_R3<br>Value: -0.18","Diagnosis: Control<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value:  0.28"],"hoverinfo":["text","text","text","text","text","text","text","text"],"frame":null},{"type":"scatterpolar","mode":"lines+markers","r":[-0.70285909536709201,0.70394958065149327,0.56558440568229151,0.53279705875603434,0.51165152171929296,0.40035749906433782,0.48410252978630425,-0.70285909536709201],"theta":["Amyloid beta 42","Tau protein","Phosphorylated tau<br />protein","Growth-regulated alpha<br />protein","Matrix<br />metalloproteinase 10","Matrix<br />metalloproteinase 7","TNF-related<br />apoptosis-inducing<br />ligand receptor 3","Amyloid beta 42"],"name":"Impaired","line":{"color":"#E87400","width":2},"marker":{"color":"#E87400","size":6,"line":{"color":"rgba(255,127,14,1)"}},"fill":"none","fillcolor":"#E8740033","text":["Diagnosis: Impaired<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value: -0.70","Diagnosis: Impaired<br>Variable: Tau protein<br>Raw name: tau<br>Value:  0.70","Diagnosis: Impaired<br>Variable: Phosphorylated tau<br />protein<br>Raw name: p_tau<br>Value:  0.57","Diagnosis: Impaired<br>Variable: Growth-regulated alpha<br />protein<br>Raw name: GRO_alpha<br>Value:  0.53","Diagnosis: Impaired<br>Variable: Matrix<br />metalloproteinase 10<br>Raw name: MMP10<br>Value:  0.51","Diagnosis: Impaired<br>Variable: Matrix<br />metalloproteinase 7<br>Raw name: MMP7<br>Value:  0.40","Diagnosis: Impaired<br>Variable: TNF-related<br />apoptosis-inducing<br />ligand receptor 3<br>Raw name: TRAIL_R3<br>Value:  0.48","Diagnosis: Impaired<br>Variable: Amyloid beta 42<br>Raw name: Ab_42<br>Value: -0.70"],"hoverinfo":["text","text","text","text","text","text","text","text"],"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
