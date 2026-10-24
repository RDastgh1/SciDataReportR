# Plot a per-cluster diagnostic value

Boxplot of any per-participant diagnostic (posterior probability,
distance to centroid, outlier score) split by cluster.

## Usage

``` r
PlotClusterDiagnostic(individual, value, title = NULL, noise_label = 0L)
```

## Arguments

- individual:

  Per-participant diagnostic table containing `Cluster`.

- value:

  Diagnostic column name.

- title:

  Optional plot title.

- noise_label:

  Cluster value treated as noise, or `NULL` to disable.

## Value

A `ggplot` object.

## Examples

``` r
# \donttest{
data(SimulatedPhenotypeData)
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
model <- CreateClusterModel_MClust(
  df_Training, paste0("Var", 1:12), method = "finalize",
  final_k = 4, final_model = 1
)
#> Warning: PlotClusterBoxplot(): `ClusterVar = "Cluster"` is numeric with 4 observed non-missing values (examples: 1, 2, 3, 4); converting locally to a factor. Few unique values do not establish categorical meaning. Explicitly convert this column to a factor if categorical treatment is intended.
PlotClusterDiagnostic(model$ProbFit$individual, "PosteriorMax")

# }
```
