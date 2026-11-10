# Statistical results and model provenance

Inferential functions expose a numerical `Results` tibble and a
`Metadata` list without replacing their existing plots, matrices, or
reporting tables. Bare plots and tables continue to print normally; use
`x$Results` and `x$Metadata` to inspect their analysis. Positional lists
of plots retain their length and numeric indexing. The compatibility
class stores these two fields as attributes for data frames, HTML
tables, and positional lists.

## Arguments

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

## Results

Standard columns are `Outcome`, `Predictor`, `OutcomeLabel`,
`PredictorLabel`, `Group`, `Term`, `Estimate`, `ConfLow`, `ConfHigh`,
`PValue`, `PValueFDR`, and `N`. Labels come from the supplied codebook
or attached variable labels. Additional columns include `TermLabel`,
`EffectType`, `ReferenceValue`, `Significant`, `Status`, and `Reason`,
plus function-specific statistics. Unsupported quantities are typed
missing values, not fabricated estimates. Diagnostic likelihood ratios
retain legitimate zero and infinite estimates. Failed requested
relationships remain visible with unavailable estimates and a reason.
Empty requests return the same column schema with zero rows.

A legacy unadjusted correlation fallback is unavailable as an adjusted
estimate: `FallbackEstimate` and `FallbackPValue` retain that
alternative, while legacy matrices and plots remain accessible and the
failed adjustment is recorded explicitly.

## Metadata

`Function`, `PackageVersion`, `Call`, and `Arguments` identify the
resolved analysis. `FormulaTemplate`, `Models`, `ReferenceLevels`, and
`Standardization` describe formulas, engines, reference/comparison
levels, and actual scaling values. `FDRFamily` describes correction
membership in words. `ModelCounts` counts primary analyses independently
of coefficient rows; auxiliary fits are recorded separately. `Failed`
identifies unavailable requests and their reasons. Existing metadata
fields are preserved.

## Compatibility

Standardization and statistical methods retain their existing
interpretation. Legacy reporting columns keep their original scales and
names. Canonical estimates and confidence limits use the same scale,
with reference value one for ratios and zero for signed differences.
`PlotForestFromTable(x$Results)` can display supported interval
estimates. Missing intervals do not imply a precise estimate. FDR
calculation remains available independently of the default raw-p
significance decision.

## Examples

``` r
set.seed(42)
df_Example <- data.frame(outcome = rnorm(40), predictor = rnorm(40))
analysis <- MakeUnivariateRegressionTable(df_Example, "outcome", "predictor")
analysis$Results
#> # A tibble: 1 × 21
#>   Outcome OutcomeLabel OutcomeFamily EffectType Predictor PredictorLabel Term   
#>   <chr>   <chr>        <chr>         <chr>      <chr>     <chr>          <chr>  
#> 1 outcome outcome      linear        Estimate   predictor predictor      predic…
#> # ℹ 14 more variables: Level <chr>, TermLabel <chr>, N <int>, Estimate <dbl>,
#> #   StdError <dbl>, ConfLow <dbl>, ConfHigh <dbl>, PValue <dbl>,
#> #   Significant <lgl>, ReferenceValue <dbl>, PValueFDR <dbl>, Group <chr>,
#> #   Status <chr>, Reason <chr>
analysis$Metadata$Models
#> # A tibble: 1 × 9
#>   Outcome Predictor Group Status Reason Partial ModelID Formula           Engine
#>   <chr>   <chr>     <chr> <chr>  <chr>  <lgl>   <chr>   <chr>             <chr> 
#> 1 outcome predictor NA    ok     NA     FALSE   model_1 outcome ~ predic… stats…
```
