# Explore dataset comparison results interactively

Create an interactive HTML dashboard from a
[`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)
result object. This function is designed for data review and
quality-control workflows. It displays high-level summary cards, a
traffic-light checks table, and an expandable variable-change explorer
that shows side-by-side old and new values for modified cells.

## Usage

``` r
ExploreDatasetComparison(
  CompareObj,
  Title = "Dataset comparison explorer",
  TopN = 10
)
```

## Arguments

- CompareObj:

  A list returned by
  [`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md).

- Title:

  Character title shown at the top of the dashboard. Default is
  `"Dataset comparison explorer"`.

- TopN:

  Integer number of example variables or records to show in previews and
  expanded sections. Default is `10`.

## Value

An
[`htmltools::tagList()`](https://rstudio.github.io/htmltools/reference/tagList.html)
object containing an interactive dashboard.

## Details

This function is intended for interactive review rather than publication
tables. It returns an HTML object that can be rendered in the RStudio
Viewer, Quarto, R Markdown, or Shiny.

## Examples

``` r
data(SampleData)

old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data
new_data$age[1:5] <- new_data$age[1:5] + 1

comparison <- CompareDatasets(old_data, new_data, keys = "id")

# Produce the interactive comparison dashboard
ExploreDatasetComparison(comparison)
#> <style>
#>         .sdr-dashboard {
#>           font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
#>           color: #263238;
#>            width: 100%;
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
#>         .sdr-card-grid {
#>           display: grid;
#>           grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
#>           gap: 12px;
#>           margin-bottom: 24px;
#>         }
#> 
#>         .sdr-card {
#>           border-radius: 14px;
#>           padding: 14px 16px;
#>           box-shadow: 0 1px 4px rgba(0,0,0,0.08);
#>         }
#> 
#>         .sdr-card-label {
#>           font-size: 13px;
#>           color: #546E7A;
#>           margin-bottom: 6px;
#>         }
#> 
#>         .sdr-card-value {
#>           font-size: 30px;
#>           font-weight: 800;
#>           line-height: 1;
#>         }
#> 
#>         .sdr-card-status {
#>           margin-top: 8px;
#>           font-size: 12px;
#>           font-weight: 700;
#>           text-transform: uppercase;
#>           letter-spacing: 0.04em;
#>           color: #455A64;
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
#>         .sdr-old-value {
#>           color: #C62828;
#>           background: #FFEBEE;
#>           padding: 2px 5px;
#>           border-radius: 5px;
#>           font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
#>         }
#> 
#>         .sdr-new-value {
#>           color: #2E7D32;
#>           background: #E8F5E9;
#>           padding: 2px 5px;
#>           border-radius: 5px;
#>           font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
#>         }
#>         </style>
#> <div class="sdr-dashboard">
#>   <div class="sdr-title">Dataset comparison explorer</div>
#>   <div class="sdr-subtitle">Interactive review of dataset version changes from CompareDatasets().</div>
#>   <div class="sdr-card-grid">
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Added records</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Removed records</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Added variables</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Removed variables</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #F9A825; background: #FFF8E1;">
#>       <div class="sdr-card-label">Modified values</div>
#>       <div class="sdr-card-value">5</div>
#>       <div class="sdr-card-status">WARNING</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Suspicious changes</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Duplicate keys</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>     <div class="sdr-card" style="border-left: 6px solid #2E7D32; background: #E8F5E9;">
#>       <div class="sdr-card-label">Skipped variables</div>
#>       <div class="sdr-card-value">0</div>
#>       <div class="sdr-card-status">PASS</div>
#>     </div>
#>   </div>
#>   <div class="sdr-section">
#>     <div class="sdr-section-title">Validation checks</div>
#>     <div class="sdr-section-subtitle">Search, filter, sort, and expand rows to inspect examples.</div>
#>     <div id="htmlwidget-ac96cb3ee4656e2e9ec3" class="reactable html-widget" style="width:100%;height:auto;"></div>
#>     <script type="application/json" data-for="htmlwidget-ac96cb3ee4656e2e9ec3">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Status":["PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","WARNING","PASS","PASS"],"Check":["Key Types","Duplicate Keys","Records Added","Records Removed","Variables Added","Variables Removed","Name Repair","Skipped Variables","Class Changes","Values Modified","Suspicious Changes","Excluded Keys"],"Count":[0,0,0,0,0,0,0,0,0,5,0,0],"Details":["Key storage classes match across dataset versions.","No duplicate key combinations detected.","No added records detected.","No removed records detected.","No added variables detected after ignoring tibble-style name-repair suffixes.","No removed variables detected after ignoring tibble-style name-repair suffixes.","No tibble-style name-repair differences detected among common variables.","No variables were skipped because of ambiguous normalized names.","No class changes detected among compared variables.","At least one cell-level value changed among compared records.","No high change-rate or class-change variables detected.","No duplicate matching keys were excluded from cell comparison."],"Examples":["None","None","None","None","None","None","None","None","None","age (5 changes; 1.5%)","None","None"]},"columns":[{"id":".details","name":"","type":null,"sortable":false,"resizable":false,"filterable":false,"searchable":false,"width":45,"align":"center","details":[{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Key Types details"]},{"name":"p","attribs":{},"children":["Key storage classes match across dataset versions."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Duplicate Keys details"]},{"name":"p","attribs":{},"children":["No duplicate key combinations detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Records Added details"]},{"name":"p","attribs":{},"children":["No added records detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Records Removed details"]},{"name":"p","attribs":{},"children":["No removed records detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Variables Added details"]},{"name":"p","attribs":{},"children":["No added variables detected after ignoring tibble-style name-repair suffixes."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Variables Removed details"]},{"name":"p","attribs":{},"children":["No removed variables detected after ignoring tibble-style name-repair suffixes."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Name Repair details"]},{"name":"p","attribs":{},"children":["No tibble-style name-repair differences detected among common variables."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Skipped Variables details"]},{"name":"p","attribs":{},"children":["No variables were skipped because of ambiguous normalized names."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Class Changes details"]},{"name":"p","attribs":{},"children":["No class changes detected among compared variables."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Values Modified details"]},{"name":"p","attribs":{},"children":["At least one cell-level value changed among compared records."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["age (5 changes; 1.5%)"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Suspicious Changes details"]},{"name":"p","attribs":{},"children":["No high change-rate or class-change variables detected."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]},{"name":"div","attribs":{"className":"sdr-detail-box"},"children":[{"name":"div","attribs":{"className":"sdr-detail-title"},"children":["Excluded Keys details"]},{"name":"p","attribs":{},"children":["No duplicate matching keys were excluded from cell comparison."]},{"name":"div","attribs":{"className":"sdr-detail-subtitle"},"children":["Examples"]},{"name":"div","attribs":{"className":"sdr-detail-examples"},"children":["None"]}]}]},{"id":"Status","name":"Status","type":"character","cell":[{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#F9A825","background":"#FFF8E1"}},"children":["● WARNING"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]},{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","font-weight":"700","color":"#2E7D32","background":"#E8F5E9"}},"children":["● PASS"]}],"width":130},{"id":"Check","name":"Check","type":"character","minWidth":190},{"id":"Count","name":"Count","type":"numeric","width":90,"align":"center"},{"id":"Details","name":"Details","type":"character","minWidth":350},{"id":"Examples","name":"Examples","type":"character","show":false,"html":true}],"filterable":true,"searchable":true,"defaultPageSize":12,"highlight":true,"bordered":true,"striped":true,"compact":true,"width":"100%","dataKey":"e22a3e7089b7690ec1f414031a636eb3"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>
#>   </div>
#>   <div class="sdr-section">
#>     <div class="sdr-section-title">Variable changes</div>
#>     <div class="sdr-section-subtitle">Expand a variable to review old and new values side by side.</div>
#>     <div id="htmlwidget-e5c8c404fe174e4c81bd" class="reactable html-widget" style="width:100%;height:auto;"></div>
#>     <script type="application/json" data-for="htmlwidget-e5c8c404fe174e4c81bd">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Variable":["age"],"Changes":[5],"OldVariable":["age"],"NewVariable":["age"],"OldClass":["numeric"],"NewClass":["numeric"],"ClassChanged":[false],"NameChanged":[false],"Compared":[true],"SkipReason":[null],"ComparedRecords":[333],"PercentChanged":[1.5],"PercentChangedLabel":["1.5%"],"ClassChangeLabel":["No class change"]},"columns":[{"id":".details","name":"","type":null,"sortable":false,"resizable":false,"filterable":false,"searchable":false,"width":45,"align":"center","details":[{"name":"Reactable","attribs":{"data":{"id":["1","2","3","4","5"],"ChangeType":["Different non-missing values","Different non-missing values","Different non-missing values","Different non-missing values","Different non-missing values"],"Change":["<span class='sdr-old-value'>52<\/span> &rarr; <span class='sdr-new-value'>53<\/span>","<span class='sdr-old-value'>61<\/span> &rarr; <span class='sdr-new-value'>62<\/span>","<span class='sdr-old-value'>77<\/span> &rarr; <span class='sdr-new-value'>78<\/span>","<span class='sdr-old-value'>97<\/span> &rarr; <span class='sdr-new-value'>98<\/span>","<span class='sdr-old-value'>73<\/span> &rarr; <span class='sdr-new-value'>74<\/span>"],"OldClass":["numeric","numeric","numeric","numeric","numeric"],"NewClass":["numeric","numeric","numeric","numeric","numeric"]},"columns":[{"id":"id","name":"id","type":"character"},{"id":"ChangeType","name":"Change type","type":"character","minWidth":190},{"id":"Change","name":"Change","type":"character","html":true,"minWidth":260},{"id":"OldClass","name":"Old class","type":"character","width":120},{"id":"NewClass","name":"New class","type":"character","width":120}],"filterable":true,"searchable":true,"defaultPageSize":10,"highlight":true,"bordered":true,"striped":true,"compact":true,"width":"100%","dataKey":"05a2089aa81b209d680ec30cca0c5f5b","static":false,"nested":true},"children":[]}]},{"id":"Variable","name":"Variable","type":"character","minWidth":220},{"id":"Changes","name":"Changes","type":"numeric","width":100,"align":"center","style":[{"fontWeight":"700","color":"#1565C0"}]},{"id":"OldVariable","name":"Old variable","type":"character","minWidth":180},{"id":"NewVariable","name":"New variable","type":"character","minWidth":180},{"id":"OldClass","name":"Old class","type":"character","width":130},{"id":"NewClass","name":"New class","type":"character","width":130},{"id":"ClassChanged","name":"ClassChanged","type":"logical","show":false},{"id":"NameChanged","name":"NameChanged","type":"logical"},{"id":"Compared","name":"Compared","type":"logical"},{"id":"SkipReason","name":"SkipReason","type":"character"},{"id":"ComparedRecords","name":"Compared records","type":"numeric","width":140,"align":"center"},{"id":"PercentChanged","name":"% Changed","type":"numeric","cell":["1.5%"],"width":120,"align":"center","style":[{"fontWeight":"700","color":"#1565C0"}]},{"id":"PercentChangedLabel","name":"PercentChangedLabel","type":"character","show":false},{"id":"ClassChangeLabel","name":"Class status","type":"character","cell":[{"name":"span","attribs":{"style":{"display":"inline-block","padding":"4px 9px","border-radius":"999px","color":"#2E7D32","background":"#E8F5E9"}},"children":["No class change"]}],"width":140}],"filterable":true,"searchable":true,"defaultSortDesc":true,"defaultSorted":[{"id":"Changes","desc":true}],"defaultPageSize":10,"highlight":true,"bordered":true,"striped":true,"compact":true,"width":"100%","dataKey":"7b56d76ae78775894829503a03973687"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>
#>   </div>
#> </div>
```
