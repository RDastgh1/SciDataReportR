# Plot Significant Correlations

Generate scatterplots for significant correlations based on a previously
generated correlation heatmap.

## Usage

``` r
plotSigCorrelations(
  data,
  CorrelationHeatmapObject,
  PVar = "PValue",
  Pthresh = 0.05,
  DataFrame = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset used to generate the scatterplots.

- CorrelationHeatmapObject:

  The output of the PlotCorrelationsHeatmap function.

- PVar:

  P-value column used to filter significance. Defaults to raw
  `"PValue"`; use `"PValueFDR"` for BH-adjusted selection.

- Pthresh:

  The significance threshold (default is 0.05).

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A list of scatterplot objects for significant correlations.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
# \donttest{
# Build a correlation heatmap, then plot the significant pairs. Keep the
# predictor and outcome sets disjoint so no variable is correlated with
# itself.
ch <- PlotCorrelationsHeatmap(
  mtcars,
  predictor_vars = c("wt", "hp", "disp"),
  outcome_vars = c("mpg", "qsec")
)

plots <- plotSigCorrelations(mtcars, ch)

# Display the first significant-correlation scatterplot
plots[[1]]

# }
```
