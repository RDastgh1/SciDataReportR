# Plot Significant Associations

Generate ggbetweenstats for significant correlations based on a
previously generated anova matrix

## Usage

``` r
plotSigAssociations(
  data,
  AnovaMatrixObject,
  PVar = "p",
  Pthresh = 0.05,
  DataFrame = lifecycle::deprecated()
)
```

## Arguments

- data:

  The dataset used to generate the scatterplots.

- AnovaMatrixObject:

  The output of the PlotAnovaRelationshipsMatrix function.

- PVar:

  The name of the column used to filter for significance (default is
  "P").

- Pthresh:

  The significance threshold (default is 0.05).

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A list of scatterplot objects for significant correlations.

## Examples

``` r
# \donttest{
# Build an ANOVA relationships matrix, then plot the significant pairs
av <- PlotAnovaRelationshipsMatrix(
  mtcars,
  CatVars = c("cyl", "gear"),
  ContVars = c("mpg", "wt", "hp")
)

plots <- plotSigAssociations(mtcars, av)

# Display the first significant-association plot
plots[[1]]

# }
```
