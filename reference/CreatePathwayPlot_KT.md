# Create Kynurenine-Tryptophan Pathway Plot

Compatibility alias for
[`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md).
Prefer
[`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md)
in new code. This KT pathway visualization is expected to move to a
future metabolomics-focused package; this compatibility alias will
remain during the transition.

## Usage

``` r
CreatePathwayPlot_KT(...)
```

## Arguments

- ...:

  Arguments passed to
  [`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md).

## Value

A ggplot2 object returned by
[`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md).

## See also

[`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md)
for the canonical function and full examples.

## Examples

``` r
results <- data.frame(
  Metabolite = c("Tryptophan", "Kynurenine", "Quinolinic Acid"),
  correlation = c(0.3, 0.1, 0.45),
  p_value = c(0.01, 0.5, 0.008),
  p_adj = c(0.05, 0.7, 0.03)
)

CreatePathwayPlot_KT(results, "Kynurenine pathway")
```
