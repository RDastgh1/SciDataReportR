# Plot dataset comparison diagnostics

Create diagnostic plots from a
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
result object. This function visualizes dataset-version changes,
including check status, summary metrics, structure changes, and
variable-level value changes.

## Usage

``` r
PlotDatasetComparison(
  CompareObj,
  Plot = c("All", "Checks", "SummaryMetrics", "StructureChanges", "VariableChanges",
    "TopChangedVariables"),
  interactive = TRUE,
  TopN = 10,
  Interactive = lifecycle::deprecated()
)
```

## Arguments

- CompareObj:

  A list returned by
  [`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md).

- Plot:

  Character value specifying which plot to return. Options are `"All"`,
  `"Checks"`, `"SummaryMetrics"`, `"StructureChanges"`,
  `"VariableChanges"`, and `"TopChangedVariables"`. Default is `"All"`.

- interactive:

  Logical; if `TRUE`, plots are converted to interactive `plotly`
  objects using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
  Default is `TRUE`.

- TopN:

  Integer number of variables or records to preview in plots and hover
  text. Default is `10`.

- Interactive:

  **Deprecated** (since 19.15.0). Use `interactive` instead.

## Value

If `Plot = "All"`, a named list of plots. Otherwise, a single plot
object. Plot objects are either `ggplot` objects or `plotly`
htmlwidgets, depending on `Interactive`.

## Details

Use this function after running
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
to quickly inspect what changed between two versions of a dataset.

## Examples

``` r
data(SampleData)

# Build two versions of a keyed dataset to compare
old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data
new_data$age[1:5] <- new_data$age[1:5] + 1

comparison <- CompareDatasets(old_data, new_data, keys = "id")

# Display a single diagnostic plot
PlotDatasetComparison(comparison, Plot = "Checks")

{"x":{"data":[{"orientation":"h","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999947,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999991],"base":[0,0,0,0,0,0,0,0,0,0,0],"x":[1,1,1,1,1,1,1,1,1,1,1],"y":[12,11,10,9,8,7,6,5,4,2,1],"text":["Check: Key Types<br>Status: PASS<br>Count: 0<br>Details: Key storage classes match across dataset versions.<br><br>Examples:<br>None","Check: Duplicate Keys<br>Status: PASS<br>Count: 0<br>Details: No duplicate key combinations detected.<br><br>Examples:<br>None","Check: Records Added<br>Status: PASS<br>Count: 0<br>Details: No added records detected.<br><br>Examples:<br>None","Check: Records Removed<br>Status: PASS<br>Count: 0<br>Details: No removed records detected.<br><br>Examples:<br>None","Check: Variables Added<br>Status: PASS<br>Count: 0<br>Details: No added variables detected after ignoring tibble-style name-repair suffixes.<br><br>Examples:<br>None","Check: Variables Removed<br>Status: PASS<br>Count: 0<br>Details: No removed variables detected after ignoring tibble-style name-repair suffixes.<br><br>Examples:<br>None","Check: Name Repair Differences<br>Status: PASS<br>Count: 0<br>Details: No tibble-style name-repair differences detected among common variables.<br><br>Examples:<br>None","Check: Variables Skipped From Cell Comparison<br>Status: PASS<br>Count: 0<br>Details: No variables were skipped because of ambiguous normalized names.<br><br>Examples:<br>None","Check: Class Changes<br>Status: PASS<br>Count: 0<br>Details: No class changes detected among compared variables.<br><br>Examples:<br>None","Check: Suspicious Changes<br>Status: PASS<br>Count: 0<br>Details: No high change-rate or class-change variables detected.<br><br>Examples:<br>None","Check: Duplicate Keys Excluded From Cell Comparison<br>Status: PASS<br>Count: 0<br>Details: No duplicate matching keys were excluded from cell comparison.<br><br>Examples:<br>None"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(46,125,50,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"PASS","legendgroup":"PASS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"h","width":0.90000000000000036,"base":0,"x":[1],"y":[3],"text":"Check: Values Modified<br>Status: WARNING<br>Count: 5<br>Details: At least one cell-level value changed among compared records.<br><br>Examples:<br>age (5; 1.5%)","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,168,37,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"WARNING","legendgroup":"WARNING","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":40.840182648401829,"r":7.3059360730593621,"b":10.958904109589042,"l":116.16438356164386},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"title":{"text":"Dataset comparison checks","font":{"color":"rgba(0,0,0,1)","family":"","size":17.534246575342465},"x":0,"xref":"paper"},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,12.6],"tickmode":"array","ticktext":["Excluded Keys","Suspicious Changes","Values Modified","Class Changes","Skipped Variables","Name Repair","Variables Removed","Variables Added","Records Removed","Records Added","Duplicate Keys","Key Types"],"tickvals":[1,2,3,4.0000000000000009,5,6,7.0000000000000009,8,9,10,11,12],"categoryorder":"array","categoryarray":["Excluded Keys","Suspicious Changes","Values Modified","Class Changes","Skipped Variables","Name Repair","Variables Removed","Variables Added","Records Removed","Records Added","Duplicate Keys","Key Types"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"title":{"text":"Status","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false,"displayModeBar":false},"source":"A","attrs":{"220d1f6188d3":{"x":{},"y":{},"fill":{},"text":{},"type":"bar"}},"cur_data":"220d1f6188d3","visdat":{"220d1f6188d3":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
