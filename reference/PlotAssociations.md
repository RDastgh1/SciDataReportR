# Plot Associations

This function generates scatter plots or box plots to visualize the
relationship between two variables.

## Usage

``` r
PlotAssociations(
  data,
  Var1,
  Var2,
  Ordinal = FALSE,
  IncludeTrendSymbol = FALSE,
  DataFrame = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  The data frame containing the variables of interest.

- Var1:

  The name of the first variable.

- Var2:

  The name of the second variable.

- Ordinal:

  Logical, indicating whether ordinal variables should be included.

- IncludeTrendSymbol:

  Logical; append a plus sign to a nominal-trend subtitle.

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Pthresh:

  Significance threshold, default `0.05`. An explicitly supplied value
  takes precedence over a legacy significance threshold argument.

- PVar:

  P-value column used for significance decisions. Defaults to raw
  `"PValue"`; select `"PValueFDR"` to use Benjamini-Hochberg adjusted
  values. Explicit legacy p-column names remain accepted. Existing,
  explicitly named FDR plots and adjusted-value columns are retained.

  Explicit `Pthresh` takes precedence over `Alpha`, `FDRAlpha`,
  `label_q`, or `adjusted_significance_threshold`. Pairwise
  `star_p = "adjusted"` retains its legacy adjustment method unless
  `PVar` is supplied; `PVar = "PValueFDR"` always selects BH.
  `star_p = "none"` suppresses legacy pairwise stars. Stars retain
  conventional p-value tiers within the requested cutoff.

## Value

A ggplot object representing the relationship between the variables.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## References

This function wraps ggstatsplot. Please cite:

Patil, I. (2021). Visualizations with statistical details: The
'ggstatsplot' approach. *Journal of Open Source Software*, 6(61), 3167.
[doi:10.21105/joss.03167](https://doi.org/10.21105/joss.03167)

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Two categorical variables (grouped bar chart)
PlotAssociations(Labelled, "Diagnosis", "Genotype")


# Two continuous variables (scatter plot with correlation)
PlotAssociations(Labelled, "age", "AXL")


# One continuous and one categorical variable (box/violin plot)
PlotAssociations(Labelled, "Diagnosis", "AXL")
```
