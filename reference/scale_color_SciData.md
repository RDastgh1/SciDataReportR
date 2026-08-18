# SciDataReportR discrete color scale

Apply the SciDataReportR qualitative palette to the `color` aesthetic in
a ggplot. Use this scale for categorical groups represented by points,
lines, outlines, or other color-mapped geometries.

## Usage

``` r
scale_color_SciData(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`ggplot2::scale_color_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

A ggplot2 discrete color scale.

## Examples

``` r
ggplot2::ggplot(
  iris,
  ggplot2::aes(
    x = Sepal.Length,
    y = Sepal.Width,
    color = Species
  )
) +
  ggplot2::geom_point() +
  scale_color_SciData()

```
