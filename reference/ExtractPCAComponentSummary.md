# Extract PCA component summaries

Extract variables contributing to each PCA component based on an
absolute loading threshold. Returns both a tidy long-format table and
compact summary tables suitable for reporting.

## Usage

``` r
ExtractPCAComponentSummary(
  PCAObject,
  loading_threshold = 0.4,
  top_n = NULL,
  use_labels = TRUE,
  html_format = TRUE
)
```

## Arguments

- PCAObject:

  Output object from CreatePCAObject().

- loading_threshold:

  Minimum absolute loading required for inclusion. Default is 0.4.

- top_n:

  Optional maximum number of contributors per component. If NULL, all
  contributors above threshold are retained.

- use_labels:

  Logical indicating whether variable labels should be used when
  available. Default TRUE.

- html_format:

  Logical indicating whether negative contributors should be formatted
  using red HTML text. Default TRUE.

## Value

A list containing:

- LongTable:

  A tidy tibble with one row per contributor.

- SummaryTable:

  A compact tibble with one row per component and comma-separated
  contributor summaries.

- SummaryTableLines:

  A compact tibble with one row per component and line-separated
  contributor summaries.

- FormattedSummaryTable:

  A formatted gt table with comma-separated contributors.

- FormattedSummaryTableLines:

  A formatted gt table with line-separated contributors.

## Details

Negative contributors can optionally be formatted in red HTML text for
improved readability in HTML tables.

## Examples

``` r
if (FALSE) { # \dontrun{
pca_obj <- CreatePCAObject(
  Data = mtcars,
  VarsToReduce = colnames(mtcars)
)

summary_obj <- ExtractPCAComponentSummary(pca_obj)

summary_obj$LongTable
summary_obj$FormattedSummaryTable
summary_obj$FormattedSummaryTableLines
} # }
```
