# Plot cluster centre profiles as a heatmap

Cluster centres across clustering variables. Values are the centres in
the frozen analysis scale, so a centred and scaled model reads directly
as standard deviations from the cohort mean.

`PlotClusterCentreProfile()` shows the same centres as connected lines,
which reads more like the SOM line map when variables have a meaningful
order.

## Usage

``` r
PlotClusterCentreHeatmap(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre"
)

PlotClusterCentreProfile(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre"
)
```

## Arguments

- centers:

  Matrix or data frame of cluster centres, one row per cluster.

- variable_labels:

  Optional display labels for the columns.

- title:

  Plot title.

- value_label:

  Legend title describing the centre scale.

## Value

A `ggplot` object, or `NULL` when centres are unavailable.

## Examples

``` r
# \donttest{
data(SimulatedPhenotypeData)
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
model <- CreateClusterModel_KMeans(
  df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
)
model$ModelInfo$plots$centre_heatmap

# }
```
