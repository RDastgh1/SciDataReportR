# SciDataReportR discrete fill scale

Apply the SciDataReportR qualitative palette to the `fill` aesthetic in
a ggplot. Use this scale for categorical groups represented by bars,
boxes, violins, areas, or other filled geometries.

## Usage

``` r
scale_fill_SciData(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

A ggplot2 discrete fill scale.

## Examples

``` r
ggplot2::ggplot(
  iris,
  ggplot2::aes(
    x = Species,
    y = Sepal.Length,
    fill = Species
  )
) +
  ggplot2::geom_boxplot() +
  scale_fill_SciData()

```
