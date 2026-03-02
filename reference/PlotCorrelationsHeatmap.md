# Plot correlations heatmap

Calculate correlations (or partial correlations if covariates are
supplied) between variables and plot them as a heatmap with significance
stars.

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

  A data frame containing variables to correlate.

- xVars:

  Character vector of x-axis variable names. If NULL, numeric variables
  are auto-selected.

- yVars:

  Character vector of y-axis variable names. If NULL, uses xVars and
  removes diagonal.

- covars:

  Character vector of covariate names. If provided, computes partial
  correlations via ppcor.

- method:

  Correlation method: "pearson", "spearman", or "kendall".

- Relabel:

  Logical. If TRUE, uses sjlabelled variable labels when present.

- Ordinal:

  Logical. If TRUE, ordinal variables are converted to numeric using
  ConvertOrdinalToNumeric().

- min_n:

  Minimum complete cases required for a correlation (or partial
  correlation).

- eps:

  Small tolerance to treat near-constant variables as constant.

## Value

A list with:

- Unadjusted: list(r, p, npairs, plot)

- FDRCorrected: list(r, p, npairs, plot)

- method, Relabel, Covariates, CovariatesMissing

## Details

This version is robust to non-syntactic column names (spaces, hyphens,
Greek letters like β), by preserving names during data.frame coercions
and using tidyselect::all_of() in pivot_longer().
