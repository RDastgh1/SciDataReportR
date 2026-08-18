# Plot Continuous Distributions

Creates rain-cloud plots (half-violin + box/median + scatter) for one or
more continuous variables, with optional group-wise colouring.

## Usage

``` r
PlotContinuousDistributions(
  data,
  variables = NULL,
  Fill = NULL,
  Relabel = TRUE,
  FacetLabelStyle = c("both", "label_only", "variable_only", "auto"),
  ncol = 3,
  Ordinal = lifecycle::deprecated(),
  TreatOrdinalAs = "Categorical",
  DataFrame = lifecycle::deprecated(),
  Variables = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame containing the variables to be plotted.

- variables:

  Character vector of column names to plot.

- Fill:

  Optional column name for grouping.

- Relabel:

  Logical; use variable labels when available.

- FacetLabelStyle:

  One of "both", "label_only", "variable_only", "auto".

- ncol:

  Number of columns in the facet grid.

- Ordinal:

  Deprecated logical compatibility option; use `TreatOrdinalAs` instead.

- TreatOrdinalAs:

  How ordinal variables are handled. `"Continuous"` includes their
  numeric score; `"Exclude"` omits them. `"Both"` is not meaningful for
  this plot and errors.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

## Value

A ggplot object.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Eight variables across three columns show wrapped labels and multi-row facets.
PlotContinuousDistributions(
  data = Labelled,
  variables = c("AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Ferritin",
                "Gamma_Interferon_induced_Monokin", "MMP7", "tau", "p_tau"),
  ncol = 3
)
#> Warning: Removed 100 rows containing non-finite outside the scale range
#> (`stat_half_ydensity()`).
#> Warning: Removed 100 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 100 rows containing missing values or values outside the scale range
#> (`geom_point_sorted()`).


# Grouped rain-clouds use the Diagnosis fill to compare distributions.
PlotContinuousDistributions(
  data = Labelled,
  variables = c("Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10", "TRAIL_R3"),
  Fill = "Diagnosis",
  ncol = 3
)
#> Warning: Removed 200 rows containing non-finite outside the scale range
#> (`stat_half_ydensity()`).
#> Warning: Removed 200 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 200 rows containing missing values or values outside the scale range
#> (`geom_point_sorted()`).
```
