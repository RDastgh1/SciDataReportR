# Create Statistics Table

Generate a table of statistics including means, standard deviations,
counts, and p-values.

## Usage

``` r
CreateStatisticsTable(data, TargetVar, Data = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame containing the variables of interest.

- TargetVar:

  The target variable for which statistics will be calculated.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A formatted HTML table displaying statistics.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
df_Statistics <- df_Labelled[, c("Diagnosis", "age", "AXL", "Ab_42", "p_tau")]

# Group means (SD), missing counts, and a p-value for every variable
statistics_table <- CreateStatisticsTable(
  data = df_Statistics,
  TargetVar = "Diagnosis"
)

# The returned kableExtra table prints directly in Quarto, R Markdown, and
# the reference page. Do not wrap it in `browsable()` a second time.
statistics_table
#> <div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:500px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
#>  <thead>
#>   <tr>
#>    <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Variable </th>
#>    <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Control (N=242) </th>
#>    <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Impaired (N=91) </th>
#>    <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Total (N=333) </th>
#>    <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> p-value </th>
#>   </tr>
#>  </thead>
#> <tbody>
#>   <tr>
#>    <td style="text-align:left;"> **Age** </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;">0.554</span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    N-Miss </td>
#>    <td style="text-align:left;"> 10.0 </td>
#>    <td style="text-align:left;"> 1.0 </td>
#>    <td style="text-align:left;"> 11.0 </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    Mean (SD) </td>
#>    <td style="text-align:left;"> 72.750 (13.261) </td>
#>    <td style="text-align:left;"> 71.778 (13.116) </td>
#>    <td style="text-align:left;"> 72.478 (13.207) </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> **AXL receptor tyrosine kinase** </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;">0.262</span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    N-Miss </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    Mean (SD) </td>
#>    <td style="text-align:left;"> 0.281 (0.461) </td>
#>    <td style="text-align:left;"> 0.343 (0.413) </td>
#>    <td style="text-align:left;"> 0.298 (0.449) </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> **Amyloid beta 42** </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: yellow !important;">&lt; 0.001</span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    N-Miss </td>
#>    <td style="text-align:left;"> 76.0 </td>
#>    <td style="text-align:left;"> 24.0 </td>
#>    <td style="text-align:left;"> 100.0 </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    Mean (SD) </td>
#>    <td style="text-align:left;"> 12.802 (1.409) </td>
#>    <td style="text-align:left;"> 11.274 (1.335) </td>
#>    <td style="text-align:left;"> 12.363 (1.549) </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> **Phosphorylated tau protein** </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;">  </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: yellow !important;">&lt; 0.001</span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    N-Miss </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> 0.0 </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;">    Mean (SD) </td>
#>    <td style="text-align:left;"> 3.942 (0.419) </td>
#>    <td style="text-align:left;"> 4.303 (0.480) </td>
#>    <td style="text-align:left;"> 4.040 (0.465) </td>
#>    <td style="text-align:left;"> <span style="border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color:  !important;"></span> </td>
#>   </tr>
#> </tbody>
#> </table></div>
# }
```
