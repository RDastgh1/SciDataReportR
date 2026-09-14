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
  covariates = NULL
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

## Value

A `feature_module_trait_obj` with association results and join audit.
