# Analyze feature module relationships with sample traits

Analyze feature module relationships with sample traits

## Usage

``` r
AnalyzeFeatureModuleTraits(
  wgcna_obj,
  trait_data,
  sample_id,
  trait_sample_id = sample_id,
  outcome = NULL,
  covariates = NULL,
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- wgcna_obj:

  A `feature_wgcna_obj`.

- trait_data:

  Trait data with a unique sample identifier.

- sample_id:

  Sample-ID column shared with the input data.

- trait_sample_id:

  Sample-ID column in `trait_data`.

- outcome:

  Trait columns to analyze.

- covariates:

  Optional covariate columns.

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

A `feature_module_trait_obj` with association results and join audit.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.
