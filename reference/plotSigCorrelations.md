# Plot Significant Correlations

Generate scatterplots for significant correlations based on a previously
generated correlation heatmap.

## Usage

``` r
plotSigCorrelations(
  DataFrame,
  CorrelationHeatmapObject,
  PVar = "P",
  Pthresh = 0.05
)
```

## Arguments

- DataFrame:

  The dataset used to generate the scatterplots.

- CorrelationHeatmapObject:

  The output of the PlotCorrelationsHeatmap function.

- PVar:

  The name of the column used to filter for significance (default is
  "P").

- Pthresh:

  The significance threshold (default is 0.05).

## Value

A list of scatterplot objects for significant correlations.
