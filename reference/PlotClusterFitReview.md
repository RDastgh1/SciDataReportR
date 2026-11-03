# Plot cluster fit-review metrics

Plot cluster fit-review metrics

## Usage

``` r
PlotClusterFitReview(
  fit_table,
  x = "Classes",
  metrics = NULL,
  group = NULL,
  title = "Candidate model review",
  metric_direction = NULL
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

- metric_direction:

  Optional named character vector giving the preferred direction
  (`"higher"` or `"lower"`) of metrics, for example `c(BIC = "lower")`.
  It overrides the directions recorded by the clustering pipeline in
  `attr(fit_table, "MetricDirection")` and the package defaults.

## Value

A `ggplot` object.

## Details

Each panel is labelled with the metric's preferred direction, such as
"BIC (lower is better)". Directions that do not depend on the backend
are built in (for example AIC lower, ICL and Entropy higher). BIC's
direction does: mclust reports it on a higher-is-better scale, while
tidyLPA (used by
[`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md))
and poLCA report it on a lower-is-better scale. Fit tables returned by
the package's clustering pipelines record the direction they used; for
other tables, supply it via `metric_direction`, otherwise the BIC panel
is shown without a direction. ICL is higher-is-better on both mclust's
and tidyLPA's scale, so in a tidyLPA fit table BIC and ICL panels point
in opposite directions.

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
