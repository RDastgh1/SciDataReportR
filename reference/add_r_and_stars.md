# Add r-values and significance stars to a correlations heatmap

Adds correlation coefficients and significance stars to a heatmap object
returned by
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).
This is a downstream annotation helper for SciDataReportR correlation
workflows, not a general-purpose plotting function.

## Usage

``` r
add_r_and_stars(
  res,
  star_from = c("existing", "raw", "fdr", "P", "P_adj", "column"),
  star_col = NULL,
  r_var = "R",
  r_digits = 2,
  r_size = 3,
  r_color = "black",
  r_nudge_y = -0.28,
  star_size = 6,
  star_color = "black",
  star_nudge_y = 0.15,
  remove_existing_stars = TRUE,
  p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  p_labels = c("***", "**", "*", "")
)
```

## Arguments

- res:

  The list returned by
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).
  The object should contain `Unadjusted$plot` and/or `FDRCorrected$plot`
  ggplot objects with plot data columns such as `R`, `P`, `P_adj`,
  `stars`, or `stars_FDR`.

- star_from:

  One of: "existing" (use whatever the chosen plot already mapped as
  stars), "raw" (use `stars` or compute from `P`), "fdr" (use
  `stars_FDR` or compute from `P_adj`), "P","P_adj" (compute from those
  columns), "column" (use `star_col`).

- star_col:

  Column name to use when star_from = "column"

- r_var:

  Column name for correlation values (default "R")

- r_digits:

  Decimal places for r labels

- r_size, star_size:

  Text sizes for r and stars

- r_color, star_color:

  Colors for r and stars

- r_nudge_y, star_nudge_y:

  Vertical nudges (r down, stars up)

- remove_existing_stars:

  If TRUE, remove any pre-existing star text layers

- p_breaks, p_labels:

  Cutpoints/labels for computing stars from p

## Value

A ggplot with r-values and stars added.

## Input requirements

`res` is expected to be the output of
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md).
Use `star_from = "raw"` for unadjusted p-value stars,
`star_from = "fdr"` for FDR-adjusted stars, or `star_from = "existing"`
to use the star mapping already present in the selected heatmap.

## Common downstream use

Add
[`geom_starcaption()`](https://rdastgh1.github.io/SciDataReportR/reference/geom_starcaption.md)
to the returned plot to explain the star annotation thresholds in
report-ready figures.
