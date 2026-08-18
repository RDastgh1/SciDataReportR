# Apply an evidence-aware p-value color scale

Maps raw p-values to either the Inferno or Viridis color palette using a
threshold-aware transformation. The transformation allocates additional
visual resolution around commonly interpreted p-value thresholds of
0.05, 0.01, and 0.001. Values above 0.05 are progressively desaturated
to reduce their visual emphasis while retaining a continuous
representation of the underlying p-values. The colorbar uses the same
warped coordinates, so its raw p-value labels remain separated around
these thresholds. Use this scale when p-values are mapped to the `color`
aesthetic.

## Usage

``` r
scale_color_pvalue(
  palette = "inferno",
  direction = -1,
  name = "P-value",
  breaks = c(1, 0.1, 0.05, 0.01, 0.001, 1e-05, 1e-08),
  labels = .format_pvalue_labels,
  limits = c(1e-08, 1),
  na.value = "grey80",
  guide = NULL,
  ...
)
```

## Arguments

- palette:

  A character value specifying the color palette. Must be `"inferno"` or
  `"viridis"`. Defaults to `"inferno"`.

- direction:

  A numeric value controlling the direction of the palette. Use `-1`,
  the default, for darker colors to indicate smaller p-values and
  stronger statistical evidence. Use `1` to reverse the palette.

- name:

  The legend title. Defaults to `"P-value"`.

- breaks:

  A numeric vector of raw p-values to display as legend breaks. Defaults
  to commonly interpreted statistical thresholds and reference values.

- labels:

  A function or character vector used to label the legend breaks. By
  default, the legend displays raw p-values rather than transformed
  values.

- limits:

  A numeric vector containing the minimum and maximum p-values
  represented by the scale. Values outside these limits are squished to
  the nearest limit. Defaults to `c(1e-8, 1)`.

- na.value:

  The color assigned to missing p-values. Defaults to `"grey80"`.

- guide:

  A guide function or guide name. The default, `NULL`, uses a taller
  colorbar so the threshold-aware breaks remain legible. Supply
  `"colourbar"` for ggplot2's standard-sized colorbar or another guide
  to override this behavior.

- ...:

  Additional arguments passed to
  [`ggplot2::continuous_scale()`](https://ggplot2.tidyverse.org/reference/continuous_scale.html).

## Value

A continuous ggplot2 color scale.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Four outcomes by six predictors: 24 standardized associations.
screen <- MakeUnivariateRegressionTable(
  data = Labelled,
  outcome_vars = c("tau", "p_tau", "AXL", "Ferritin"),
  predictor_vars = c(
    "age", "sex", "Adiponectin", "Cortisol", "Insulin", "Leptin"
  ),
  Standardize = TRUE
)

# A volcano plot is the natural home for this scale: effect size on x,
# evidence on y, and the colour reinforcing where the thresholds fall so
# the eye is not left to judge distance along a log axis.
ggplot2::ggplot(
  screen$Results,
  ggplot2::aes(
    x = Estimate,
    y = -log10(PValue),
    color = PValue
  )
) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_pvalue() +
  ggplot2::labs(
    title = "Standardized associations across four outcomes",
    subtitle = "Dashed line: p = 0.05. Dotted line: no effect.",
    x = "Estimate per SD", y = expression(-log[10](p))
  ) +
  ggplot2::theme_bw()

# }
```
