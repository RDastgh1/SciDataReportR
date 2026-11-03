# Plot labelled numeric profiles by cluster

Plot labelled numeric profiles by cluster

## Usage

``` r
PlotClusterProfiles(data, ClusterVar, variables, codebook = NULL)
```

## Arguments

- data:

  A data frame.

- ClusterVar:

  Character string naming the cluster/grouping variable.

- variables:

  Character vector of variable names to plot.

- codebook:

  Optional codebook data frame with columns `Variable` and `Label`.

## Value

A `ggplot` object.

## Examples

``` r
# \donttest{
data(SimulatedPhenotypeData)
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
model <- CreateClusterModel_KMeans(
  df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
)
PlotClusterProfiles(model$DataWithClusters, "Cluster", paste0("Var", 1:12))

# }
```
