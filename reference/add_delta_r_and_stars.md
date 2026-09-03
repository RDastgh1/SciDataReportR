# Add DeltaR values and significance stars to a correlation comparison heatmap

Adds correlation differences and significance stars to a heatmap
returned by
[`PlotCorrelationComparisons()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationComparisons.md).
This is the correlation-comparison counterpart to
[`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md).

## Usage

``` r
add_delta_r_and_stars(
  res,
  star_from = c("fdr", "raw"),
  delta_digits = 2,
  delta_size = 3,
  star_size = 5,
  delta_color = "black",
  star_color = "black",
  delta_nudge_y = -0.18,
  star_nudge_y = 0.2,
  remove_existing_stars = TRUE
)
```

## Arguments

- res:

  An object returned by
  [`PlotCorrelationComparisons()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationComparisons.md).

- star_from:

  Whether stars should use `"fdr"` or `"raw"` comparison p-values.

- delta_digits:

  Number of decimal places used for DeltaR.

- delta_size:

  Text size for DeltaR labels.

- star_size:

  Text size for significance stars.

- delta_color:

  Color for DeltaR labels.

- star_color:

  Color for significance stars.

- delta_nudge_y:

  Vertical position adjustment for DeltaR.

- star_nudge_y:

  Vertical position adjustment for stars.

- remove_existing_stars:

  Logical. Remove the original star-only layer before adding the
  combined annotations.

## Value

A ggplot containing DeltaR values and significance stars.

## Examples

``` r
# res <- PlotCorrelationComparisons(...)
# add_delta_r_and_stars(res)
```
