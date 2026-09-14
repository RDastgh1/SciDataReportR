# Transform selected features in a wide data frame

Transform selected features in a wide data frame

## Usage

``` r
TransformFeatureMatrix(
  data,
  variables = NULL,
  method = c("none", "log1p", "log2", "log10", "glog"),
  pseudocount = 1,
  glog_lambda = 1
)
```

## Arguments

- data:

  A data frame with samples in rows.

- variables:

  Character vector of numeric feature columns.

- method:

  One of `"none"`, `"log1p"`, `"log2"`, `"log10"`, or `"glog"`.

- pseudocount:

  Numeric value added before plain-log transformations.

- glog_lambda:

  Positive stabilization value for generalized log.

## Value

A list with transformed `data`, a feature-level `summary`, and a
before/after skewness `plot`.
