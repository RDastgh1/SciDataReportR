# Plot a per-participant silhouette profile

The classic silhouette profile: one bar per participant, sorted within
cluster, with the average silhouette width marked. Bars near zero or
below sit closer to a neighbouring cluster than their own.

## Usage

``` r
PlotClusterSilhouette(silhouette, title = "Silhouette profile")
```

## Arguments

- silhouette:

  A
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html)
  object or a matrix with `cluster`, `neighbor`, and `sil_width`
  columns.

- title:

  Plot title.

## Value

A `ggplot` object, or `NULL` when silhouette widths are unavailable.

## Examples

``` r
# \donttest{
data(SimulatedPhenotypeData)
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
model <- CreateClusterModel_KMeans(
  df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
)
model$ModelInfo$plots$silhouette

# }
```
