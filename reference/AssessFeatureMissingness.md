# Assess missingness in a wide feature matrix

Assess missingness in a wide feature matrix

## Usage

``` r
AssessFeatureMissingness(
  data,
  variables = NULL,
  missing_threshold = 0.5,
  intensity_heuristic = TRUE,
  intensity_cutoff = 0.3
)
```

## Arguments

- data:

  A data frame with samples in rows.

- variables:

  Character vector of numeric feature columns.

- missing_threshold:

  Numeric proportion above which a feature is flagged.

- intensity_heuristic:

  Logical. Whether to calculate the exploratory
  missingness-versus-sample-intensity correlation.

- intensity_cutoff:

  Minimum absolute correlation for an intensity-dependent flag.

## Value

A list with a feature-level `summary` tibble and `plot`.
