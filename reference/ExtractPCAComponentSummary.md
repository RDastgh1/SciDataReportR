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

Use this after
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md),
at the point where the components exist but do not yet mean anything.
PCA hands back axes named `PC1`, `PC2`, `PC3`

- each a weighted combination of every input variable - and the analysis
  cannot be written up until those names are replaced by something
  interpretable.

The full loadings matrix is the wrong tool for that job: with 40
variables and 8 components it is 320 numbers, most of them near zero and
irrelevant. This function keeps only the loadings whose absolute value
clears `loading_threshold` (0.4 by default), which are the variables
actually driving each component, and reports them per component, sorted
by strength. Reading down the retained variables is what tells you that
PC1 is "an inflammatory axis" or that PC2 separates volume from
thickness.

The sign matters as much as the magnitude. A component with some
variables loading positively and others negatively is a *contrast* - it
scores the balance between two sets of measures rather than their
overall level - and `html_format = TRUE` prints the negative
contributors in red so that structure is visible at a glance instead of
having to be hunted for in a column of numbers.

Once a component has a name, that name is what belongs on the axis of
every downstream plot and in every table of PCA scores.

## Choosing a threshold

`loading_threshold` trades completeness against readability. Raise it
(0.5-0.6) when a component retains so many variables that no theme is
visible; lower it (0.3) when a component comes back nearly empty, which
usually means its variance is spread thinly across many measures rather
than concentrated in a few. `top_n` caps the list per component instead,
which is the better control when what you want is a compact table for a
manuscript rather than a different scientific claim.

## See also

[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md)
to fit the PCA,
[`CreatePCATable()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md)
for the variance-explained table, and
[`ProjectPCA()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectPCA.md)
to score new data on the components once they have been interpreted.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
vars_Biomarkers <- c(
  "AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
  "Apolipoprotein_A1", "Apolipoprotein_B", "C_Reactive_Protein",
  "Cortisol", "Cystatin_C", "Ferritin", "Insulin", "Leptin", "p_tau"
)

# Thirteen correlated biomarkers reduced to a handful of components
pca_obj <- CreatePCAObject(
  data = df_Labelled,
  VarsToReduce = vars_Biomarkers
)

# At this point the components are called RC1, RC2, RC3
summary_obj <- ExtractPCAComponentSummary(pca_obj)

# One row per component, listing only the variables that drive it
summary_obj$FormattedSummaryTableLines


  

Component
```
