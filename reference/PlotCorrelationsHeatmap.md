# Plot correlations heatmap

Calculates correlations (or partial correlations) and plots a heatmap.
Fully label-aware, robust to non-syntactic names, and safe for
real-world data.

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

  A data frame.

- xVars:

  Character vector of x variables.

- yVars:

  Character vector of y variables.

- covars:

  Optional covariates.

- method:

  Correlation method.

- Relabel:

  Use labels if available.

- Ordinal:

  Convert ordinal variables.

- min_n:

  Minimum N required.

- eps:

  Variance tolerance.

## Value

A list with correlation matrices and plots.
