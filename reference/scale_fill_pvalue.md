# Apply an evidence-aware p-value fill scale

Maps raw p-values to either the Inferno or Viridis color palette using a
threshold-aware transformation. The transformation allocates additional
visual resolution around commonly interpreted p-value thresholds of
0.05, 0.01, and 0.001. Values above 0.05 are progressively desaturated
to reduce their visual emphasis while retaining a continuous
representation of the underlying p-values. The colorbar uses the same
warped coordinates, so its raw p-value labels remain separated around
these thresholds. Use this scale when p-values are mapped to the `fill`
aesthetic.

## Usage

``` r
scale_fill_pvalue(
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

  The fill color assigned to missing p-values. Defaults to `"grey80"`.

- guide:

  A guide function or guide name. The default, `NULL`, uses a taller
  colorbar so the threshold-aware breaks remain legible. Supply
  `"colourbar"` for ggplot2's standard-sized colorbar or another guide
  to override this behavior.

- ...:

  Additional arguments passed to
  [`ggplot2::continuous_scale()`](https://ggplot2.tidyverse.org/reference/continuous_scale.html).

## Value

A continuous ggplot2 fill scale.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A univariate screen: each biomarker against diagnosis, standardized so the
# effect sizes are comparable.
screen <- MakeUnivariateRegressionTable(
  data = Labelled,
  outcome_vars = "Diagnosis",
  predictor_vars = c(
    "age", "sex", "AXL", "Adiponectin", "Cortisol",
    "Ferritin", "Insulin", "Leptin", "tau", "p_tau"
  ),
  Standardize = TRUE
)

# Bar length is the effect, fill is the evidence for it. Reading the two
# together is the point: a long pale bar is a large estimate nobody should
# rely on, and a short dark bar is a small effect that is real.
ggplot2::ggplot(
  screen$Results,
  ggplot2::aes(
    x = Estimate,
    y = stats::reorder(TermLabel, Estimate),
    fill = PValue
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_pvalue() +
  ggplot2::labs(
    title = "Association with impairment, per standard deviation",
    subtitle = "Dashed line: odds ratio of 1, no association",
    x = "Odds ratio per SD", y = NULL
  ) +
  ggplot2::theme_bw()

# }
```
