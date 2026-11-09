# Plot correlations heatmap

Computes correlations or partial correlations and plots a heatmap.
Handles:

- continuous + categorical covariates

- labelled data

- non-syntactic names

- sparse real-world datasets

- ordinal variables

- partial correlations via residualization

## Usage

``` r
PlotCorrelationsHeatmap(
  data,
  predictor_vars = NULL,
  outcome_vars = NULL,
  covariates = NULL,
  method = "pearson",
  Relabel = TRUE,
  Ordinal = lifecycle::deprecated(),
  TreatOrdinalAs = "Categorical",
  min_n = 3,
  eps = 1e-12,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  triangle = c("full", "upper", "lower"),
  IncludeTrendSymbol = FALSE,
  Data = lifecycle::deprecated(),
  xVars = lifecycle::deprecated(),
  yVars = lifecycle::deprecated(),
  covars = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  data.frame

- predictor_vars:

  character vector

- outcome_vars:

  character vector

- covariates:

  optional covariates

- method:

  pearson/spearman/kendall

- Relabel:

  use labels

- Ordinal:

  Deprecated logical compatibility option; use `TreatOrdinalAs` instead.

- TreatOrdinalAs:

  How ordinal variables are handled. `"Continuous"` includes ordinal
  scores and `"Exclude"` omits them.

- min_n:

  minimum complete rows

- eps:

  variance tolerance

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  With `"matrix"`, FDR correction is applied across the whole p-value
  matrix at once (historical behavior). With `"per_outcome"`, correction
  is applied separately within each outcome: in this function outcomes
  are the columns of the p-value matrix, i.e. `outcome_vars`
  (`outcome_margin = 2`).

- cluster_rows, cluster_columns:

  Logical; cluster predictor rows and/or outcome columns from their
  displayed correlation profiles. Defaults retain the caller-supplied
  variable order.

- triangle:

  Display `"full"` (default), `"upper"`, or `"lower"` half of a
  symmetric correlation matrix. This affects the plot only; returned
  matrices remain complete. Triangle display is applied only when the
  same variables occur on both axes.

- IncludeTrendSymbol:

  Logical; add a plus sign for `0.05 < p <= 0.10`.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- xVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- yVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- covars:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

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

A list. `Unadjusted` and `FDRCorrected` each contain `r`, `p`, `npairs`,
and `plot`. The standardized aliases `p` (same as `Unadjusted`) and
`p_fdr` (same as `FDRCorrected`) are also included.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Square heatmap: the same 10 variables on both axes
vars <- c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin",
          "Alpha_2_Macroglobulin", "Apolipoprotein_A1", "Apolipoprotein_B",
          "C_Reactive_Protein", "Cortisol", "Insulin")

square <- PlotCorrelationsHeatmap(
  Labelled,
  predictor_vars = vars,
  outcome_vars = vars
)
square$Unadjusted$plot

square$FDRCorrected$plot


# Show each symmetric pair once; estimates and p-value matrices stay full.
upper <- PlotCorrelationsHeatmap(
  Labelled,
  predictor_vars = vars,
  outcome_vars = vars,
  triangle = "upper"
)
upper$FDRCorrected$plot


# Rectangular heatmap: different variables on x and y
rectangular <- PlotCorrelationsHeatmap(
  Labelled,
  predictor_vars = c("age", "AXL", "Adiponectin", "Cortisol", "Insulin"),
  outcome_vars = c("Apolipoprotein_A1", "Apolipoprotein_B",
                   "C_Reactive_Protein", "Ferritin", "Leptin")
)
rectangular$Unadjusted$plot
```
