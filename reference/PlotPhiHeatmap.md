# Plot Phi Correlations Between Binary Variables

Computes pairwise phi coefficients between binary categorical variables
with explicit 0/1 coding (1 == PositiveLevel from
[`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md)),
then renders heatmap-style plots with raw and FDR-adjusted significance.

## Usage

``` r
PlotPhiHeatmap(
  data,
  CatVars,
  Relabel = TRUE,
  binary_map = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  triangle = c("full", "upper", "lower"),
  IncludeTrendSymbol = FALSE,
  Data = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  A dataframe.

- CatVars:

  Character vector of binary categorical variable names.

- Relabel:

  Logical; if TRUE, uses sjlabelled variable labels for axes.

- binary_map:

  Optional mapping as returned by
  [`createBinaryMapping()`](https://rdastgh1.github.io/SciDataReportR/reference/createBinaryMapping.md).
  If NULL, a mapping is created internally for `CatVars`.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all p-values at once (historical behavior).
  `"per_outcome"` corrects separately within each y-axis variable
  (`YVar`); the Phi matrix is symmetric, so this treats each variable's
  row of tiles as one family.

- cluster_rows, cluster_columns:

  Logical; cluster y-axis rows and/or x-axis columns using displayed
  Phi-coefficient profiles.

- triangle:

  Display `"full"` (default), `"upper"`, or `"lower"` half of the
  symmetric Phi matrix. This affects the plot only; returned results
  remain complete.

- IncludeTrendSymbol:

  Logical; add a plus sign for `0.05 < p <= 0.10`.

- Data:

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

A list with:

- `Unadjusted`: list(PvalTable, plot)

- `FDRCorrected`: list(PvalTable, plot)

- `method` = "Phi"

- `Relabel`

- `BinaryMapping` (used)

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# CatVars must be binary (exactly two unique non-NA values). Derive a few
# more binary indicators so the matrix has off-diagonal structure to read;
# self-associations on the diagonal are masked out.
Labelled$APOE4 <- ifelse(
  grepl("E4", as.character(Labelled$Genotype)), "Carrier", "Non-carrier")
Labelled$HighTau <- ifelse(
  Labelled$tau > stats::median(Labelled$tau, na.rm = TRUE), "High", "Low")
Labelled$LowAbeta <- ifelse(
  Labelled$Ab_42 < stats::median(Labelled$Ab_42, na.rm = TRUE), "Low", "High")

result <- PlotPhiHeatmap(
  Labelled,
  CatVars = c("Diagnosis", "sex", "APOE4", "HighTau", "LowAbeta")
)

# Raw p-value phi heatmap
result$Unadjusted$plot


# FDR-adjusted phi heatmap
result$FDRCorrected$plot


# Show each binary-variable pair once.
upper <- PlotPhiHeatmap(
  Labelled,
  CatVars = c("Diagnosis", "sex", "APOE4", "HighTau", "LowAbeta"),
  triangle = "upper"
)
upper$FDRCorrected$plot
```
