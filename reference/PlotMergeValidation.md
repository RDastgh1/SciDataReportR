# Plot merge validation diagnostics

Create diagnostic plots from a
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
result object. This function visualizes key merge-audit outputs,
including validation check status, key coverage, join-variable auditing,
duplicate-variable agreement, and duplicate-variable conflict counts.
Use this after running
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md)
to quickly inspect whether a merged dataset appears trustworthy.

## Usage

``` r
PlotMergeValidation(
  MergeObj,
  Plot = c("All", "Checks", "Coverage", "JoinAudit", "Agreement", "Conflicts"),
  interactive = TRUE,
  Interactive = lifecycle::deprecated()
)
```

## Arguments

- MergeObj:

  A list returned by
  [`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).

- Plot:

  Character value specifying which plot to return. Options are `"All"`,
  `"Checks"`, `"Coverage"`, `"JoinAudit"`, `"Agreement"`, and
  `"Conflicts"`. Default is `"All"`.

- interactive:

  Logical; if `TRUE`, plots are converted to interactive `plotly`
  objects using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
  Default is `TRUE`.

- Interactive:

  **Deprecated** (since 19.15.0). Use `interactive` instead.

## Value

If `Plot = "All"`, a named list of plots. Otherwise, a single plot
object. Plot objects are either `ggplot` objects or `plotly`
htmlwidgets, depending on `Interactive`.

## Details

The function expects the object returned by
[`ValidateMerge()`](https://rdastgh1.github.io/SciDataReportR/reference/ValidateMerge.md).
It does not re-run any merge validation checks.

## Examples

``` r
set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Display a single diagnostic plot
PlotMergeValidation(validation, Plot = "Checks")

{"x":{"data":[{"orientation":"h","width":[0.89999999999999858,0.89999999999999858,0.89999999999999947,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999991],"base":[0,0,0,0,0,0,0,0,0,0],"x":[1,1,1,1,1,1,1,1,1,1],"y":[10,9,8,7,6,5,4,3,2,1],"text":["Check: Key Types<br>Status: PASS<br>Count: 0<br>Details: Key storage classes match across datasets.","Check: Missing Keys<br>Status: PASS<br>Count: 0<br>Details: No missing key rows detected.","Check: Duplicate Keys<br>Status: PASS<br>Count: 0<br>Details: No duplicate complete key combinations detected.","Check: Coverage<br>Status: PASS<br>Count: 0<br>Details: All complete source key combinations overlap.","Check: Row Inflation<br>Status: PASS<br>Count: 1<br>Details: No meaningful row inflation detected.","Check: Overlapping Variables<br>Status: PASS<br>Count: 0<br>Details: No non-key variables overlap across source datasets.","Check: Unresolved Duplicate Variables<br>Status: PASS<br>Count: 0<br>Details: No unresolved duplicate variable pairs detected.","Check: Variable Conflicts<br>Status: PASS<br>Count: 0<br>Details: No duplicated-variable value conflicts detected.","Check: Suspicious Conflicts<br>Status: PASS<br>Count: 0<br>Details: No low-agreement or class-mismatched duplicated variables detected.","Check: Merge Readiness<br>Status: PASS<br>Count: 0<br>Details: No major merge-integrity blockers detected."],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(46,125,50,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"PASS","legendgroup":"PASS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":40.840182648401829,"r":7.3059360730593621,"b":10.958904109589042,"l":186.30136986301375},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"title":{"text":"Merge validation checks","font":{"color":"rgba(0,0,0,1)","family":"","size":17.534246575342465},"x":0,"xref":"paper"},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,10.6],"tickmode":"array","ticktext":["Merge Readiness","Suspicious Conflicts","Variable Conflicts","Unresolved Duplicate Variables","Overlapping Variables","Row Inflation","Coverage","Duplicate Keys","Missing Keys","Key Types"],"tickvals":[1,2,3,4,5,6.0000000000000009,7,8,9,10],"categoryorder":"array","categoryarray":["Merge Readiness","Suspicious Conflicts","Variable Conflicts","Unresolved Duplicate Variables","Overlapping Variables","Row Inflation","Coverage","Duplicate Keys","Missing Keys","Key Types"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"title":{"text":"Status","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"257b47f27dd6":{"x":{},"y":{},"fill":{},"text":{},"type":"bar"}},"cur_data":"257b47f27dd6","visdat":{"257b47f27dd6":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
