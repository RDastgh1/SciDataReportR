# Plot cluster fit-review metrics

Plot cluster fit-review metrics

## Usage

``` r
PlotClusterFitReview(
  fit_table,
  x = "Classes",
  metrics = NULL,
  group = NULL,
  title = "Candidate model review"
)
```

## Arguments

- fit_table:

  Candidate-model fit table.

- x:

  Candidate-count column name.

- metrics:

  Numeric metric columns to display. When `NULL`, a concise set of
  decision-relevant raw metrics is selected automatically.

- group:

  Optional candidate-family column, such as `"Model"` for Mclust or
  `"Epsilon"` for HDBSCAN. When `NULL`, it is inferred when possible.

- title:

  Plot title.

## Value

A `ggplot` object.

## Examples

``` r
# \donttest{
data(SimulatedPhenotypeData)
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
model <- CreateClusterModel_KMeans(
  df_Training, paste0("Var", 1:12), k_range = 2:4
)
PlotClusterFitReview(model$ModelInfo$fit_table)

# }
```
