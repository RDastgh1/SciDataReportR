# Plot differential feature abundance as a volcano plot

Plot differential feature abundance as a volcano plot

## Usage

``` r
PlotFeatureDifferentialAbundance(
  x,
  contrast = NULL,
  significance_metric = c("FDR", "PValue"),
  alpha = NULL,
  min_effect = 0,
  label_top_n = 10,
  label_by = "Variable"
)
```

## Arguments

- x:

  A `feature_differential_abundance` object, as returned by
  [`AnalyzeFeatureDifferentialAbundance()`](https://rdastgh1.github.io/SciDataReportR/reference/AnalyzeFeatureDifferentialAbundance.md)
  or one of its omics wrappers.

- contrast:

  Contrast to plot. Defaults to the first.

- significance_metric:

  Metric used for coloring and the threshold line: `"FDR"` (default) or
  `"PValue"`. The y-axis always shows the P-value.

- alpha:

  Significance threshold. Defaults to the value used in the analysis.

- min_effect:

  Optional minimum absolute effect for a feature to be colored as
  significant.

- label_top_n:

  Number of most significant features to label.

- label_by:

  Column of the results used for labels. Defaults to `"Variable"`.

## Value

A ggplot object.
