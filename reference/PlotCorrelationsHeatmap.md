# Plot correlations heatmap

Computes correlations or partial correlations and plots a heatmap.
Handles:

- continuous + categorical covariates

- labelled data

- non-syntactic names

- sparse real-world datasets

- ordinal variables

- partial correlations via residualization

## Usage

``` r
PlotCorrelationsHeatmap(
  Data,
  xVars = NULL,
  yVars = NULL,
  covars = NULL,
  method = "pearson",
  Relabel = TRUE,
  Ordinal = FALSE,
  min_n = 3,
  eps = 1e-12
)
```

## Arguments

- Data:

  data.frame

- xVars:

  character vector

- yVars:

  character vector

- covars:

  optional covariates

- method:

  pearson/spearman/kendall

- Relabel:

  use labels

- Ordinal:

  include ordinal vars

- min_n:

  minimum complete rows

- eps:

  variance tolerance

## Value

list
