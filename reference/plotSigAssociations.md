# Plot Significant Associations

Generate ggbetweenstats for significant correlations based on a
previously generated anova matrix

## Usage

``` r
plotSigAssociations(DataFrame, AnovaMatrixObject, PVar = "p", Pthresh = 0.05)
```

## Arguments

- DataFrame:

  The dataset used to generate the scatterplots.

- AnovaMatrixObject:

  The output of the PlotAnovaRelationshipsMatrix function.

- PVar:

  The name of the column used to filter for significance (default is
  "P").

- Pthresh:

  The significance threshold (default is 0.05).

## Value

A list of scatterplot objects for significant correlations.
