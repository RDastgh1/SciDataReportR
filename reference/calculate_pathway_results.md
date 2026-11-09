# Calculate Pathway Results for Metabolite Comparisons

Calculate Pathway Results for Metabolite Comparisons

## Usage

``` r
calculate_pathway_results(
  data,
  comparison_var,
  covariates = NULL,
  metabolites,
  comparison_type = "auto",
  use_point_correlation = FALSE,
  df = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  Data frame containing metabolite and comparison data

- comparison_var:

  Character string specifying the comparison variable name

- covariates:

  Character vector of covariate names (optional)

- metabolites:

  Character vector of metabolite names to analyze

- comparison_type:

  Character string: "auto", "binary", or "continuous"

- use_point_correlation:

  Logical, if TRUE uses point correlation for binary comparisons

- df:

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

Data frame with metabolite results including fold change or correlation
values

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- calculate_pathway_results(
  data = my_data,
  comparison_var = "poorSleep",
  covariates = c("Age", "BMI"),
  metabolites = c("Tryptophan", "Kynurenine"),
  comparison_type = "binary"
)
} # }
```
