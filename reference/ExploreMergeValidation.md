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

## Examples

``` r
set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Produce the interactive merge-validation dashboard
ExploreMergeValidation(validation)
#> <style>
#>         .sdr-dashboard {
#>           font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
#>           color: #263238;
#>           width: 100%;
#>   max-width: none;
#>   overflow-x: auto;
#>         }
#> 
#>         .sdr-title {
#>           font-size: 28px;
#>           font-weight: 750;
#>           margin: 0 0 6px 0;
#>         }
#> 
#>         .sdr-subtitle {
#>           color: #607D8B;
#>           margin-bottom: 22px;
#>           font-size: 14px;
#>         }
#> 
#>         .sdr-section {
#>           margin-top: 28px;
#>           margin-bottom: 10px;
#>         }
#> 
#>         .sdr-section-title {
#>           font-size: 21px;
#>           font-weight: 750;
#>           margin-bottom: 4px;
#>         }
#> 
#>         .sdr-section-subtitle {
#>           font-size: 13px;
#>           color: #607D8B;
#>           margin-bottom: 12px;
#>         }
#> 
#>         .sdr-accordion {
#>           margin-top: 28px;
#>           margin-bottom: 10px;
#>           border: 1px solid #E0E6ED;
#>           border-radius: 14px;
#>           padding: 10px 16px;
#>           background: #F8FAFC;
#>         }
#> 
#>         .sdr-accordion[open] {
#>           background: #FFFFFF;
#>         }
#> 
#>         .sdr-accordion-summary {
#>           font-size: 21px;
#>           font-weight: 750;
#>           cursor: pointer;
#>           padding: 4px 0;
#>         }
#> 
#>         .sdr-accordion-summary:hover {
#>           color: #1565C0;
#>         }
#> 
#>         .sdr-detail-box {
#>           background: #FAFAFA;
#>           border-left: 4px solid #1565C0;
#>           padding: 12px 14px;
#>           margin: 8px;
#>           border-radius: 8px;
#>         }
#> 
#>         .sdr-detail-title {
#>           font-weight: 750;
#>           margin-bottom: 8px;
#>         }
#> 
#>         .sdr-detail-subtitle {
#>           font-weight: 700;
#>           margin-top: 10px;
#>           margin-bottom: 4px;
#>           color: #455A64;
#>         }
#> 
#>         .sdr-detail-examples {
#>           font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
#>           font-size: 12px;
#>           line-height: 1.45;
#>           color: #37474F;
#>         }
#> 
#>         .sdr-left-value {
#>           color: #C62828;
#>           background: #FFEBEE;
#>           padding: 2px 5px;
#>           border-radius: 5px;
#>           font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
#>         }
#> 
#>         .sdr-right-value {
#>           color: #2E7D32;
#>           background: #E8F5E9;
#>           padding: 2px 5px;
#>           border-radius: 5px;
#>           font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
#>         }
#>         </style>
#> <div class="sdr-dashboard">
#>   <div class="sdr-title">Merge validation explorer</div>
#>   <div class="sdr-subtitle">Interactive review of merge integrity from ValidateMerge().</div>
#>   <div class="sdr-section">
#>     <div class="sdr-section-title">Validation checks</div>
#>     <div class="sdr-section-subtitle">Search, filter, sort, and expand rows to inspect merge-integrity examples. INFO rows summarize rows, columns, and unique keys across the source and merged datasets.</div>
#>     <div id="htmlwidget-36aa3d2a04d42bbc2145" class="reactable html-widget" style="width:100%;height:350px;"></div>
#>     <script type="application/json" data-for="htmlwidget-36aa3d2a04d42bbc2145">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Status":["INFO","INFO","INFO","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS"],"Check":["Rows","Columns","Unique keys","Key Types","Missing Keys","Expected Relationship","Duplicate Keys","Coverage","Row Count","Row Inflation","Overlapping variables","Unresolved duplicates","Variable conflicts","Suspicious conflicts","Merge readiness"],"Count":[50,3,50,0,0,0,0,0,0,1,0,0,0,0,0],"Details":["Left: 50; Right: 50; Merged: 50.","Left: 2; Right: 2; Merged: 3.","Unique key combinations - Left: 50; Right: 50; Merged: 50.","Key storage classes match across datasets.","No missing key rows detected.","Detected relationship matches expected_relationship = 'one-to-one'.","No duplicate key blockers detected for the expected relationship.","All complete source key combinations overlap.","MergedData row count matches LeftData row count.","No meaningful row inflation detected.","No non-key variables overlap across source datasets.","No unresolved duplicate variable pairs detected.","No duplicated-variable value conflicts detected.","No low-agreement or class-mismatched duplicated variables detected.","No major merge-integrity blockers detected."],"Examples":["None","None","None","None","None","None","None","<strong>Matching examples:<\/strong><br>1<br>2<br>3<br>4<br>5<br>6<br>7<br>8<br>9<br>10<br><br><strong>Left only examples:<\/strong><br>None<br><br><strong>Right only examples:<\/strong><br>None","None","Row inflation factor: 1","None","None","None","None","No major merge-integrity blockers detected."]},"columns":[{"id":".details","name":"","type":null,"sortable":false,"resizable":false,"filterable":false,"searchable":false,"width":45,"align":"center","details":[null,null,null,{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Key Types details"]},{"name":"p","attribs":{},"children":["Key storage classes match across datasets."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Missing Keys details"]},{"name":"p","attribs":{},"children":["No missing key rows detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Expected Relationship details"]},{"name":"p","attribs":{},"children":["Detected relationship matches expected_relationship = 'one-to-one'."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Duplicate Keys details"]},{"name":"p","attribs":{},"children":["No duplicate key blockers detected for the expected relationship."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Coverage details"]},{"name":"p","attribs":{},"children":["All complete source key combinations overlap."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["<strong>Matching examples:<\/strong><br>1<br>2<br>3<br>4<br>5<br>6<br>7<br>8<br>9<br>10<br><br><strong>Left only examples:<\/strong><br>None<br><br><strong>Right only examples:<\/strong><br>None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Row Count details"]},{"name":"p","attribs":{},"children":["MergedData row count matches LeftData row count."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Row Inflation details"]},{"name":"p","attribs":{},"children":["No meaningful row inflation detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["Row inflation factor: 1"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Overlapping variables details"]},{"name":"p","attribs":{},"children":["No non-key variables overlap across source datasets."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Unresolved duplicates details"]},{"name":"p","attribs":{},"children":["No unresolved duplicate variable pairs detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Variable conflicts details"]},{"name":"p","attribs":{},"children":["No duplicated-variable value conflicts detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Suspicious conflicts details"]},{"name":"p","attribs":{},"children":["No low-agreement or class-mismatched duplicated variables detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Merge readiness details"]},{"name":"p","attribs":{},"children":["No major merge-integrity blockers detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["No major merge-integrity blockers detected."]}]}]},{"id":"Status","name":"Status","type":"character","cell":[{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#1565C0","background":"#E3F2FD"}},"children":["● INFO"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#1565C0","background":"#E3F2FD"}},"children":["● INFO"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#1565C0","background":"#E3F2FD"}},"children":["● INFO"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]}],"width":130},{"id":"Check","name":"Check","type":"character","minWidth":190},{"id":"Count","name":"Count","type":"numeric","width":90,"align":"center"},{"id":"Details","name":"Details","type":"character","minWidth":350},{"id":"Examples","name":"Examples","type":"character","show":false,"html":true}],"filterable":true,"searchable":true,"pagination":false,"highlight":true,"bordered":true,"striped":true,"compact":true,"width":"100%","height":"350px","dataKey":"0aae257d1db5692111f7c1a0b29b9937"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>
#>   </div>
#> </div>
```
