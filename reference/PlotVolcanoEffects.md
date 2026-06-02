# Plot volcano-style association effects

Screen a set of predictor variables against one outcome and visualize
the association results as volcano plots. The function supports
continuous outcomes using standardized beta values, and two-group
categorical outcomes using either Cohen-style standardized effects or
log2 fold change effects. Optional covariates are included in the model
for each predictor.

## Usage

``` r
PlotVolcanoEffects(
  Data,
  xVars,
  yVar,
  Covariates = NULL,
  OutcomeType = c("auto", "continuous", "categorical"),
  EffectMetric = c("auto", "cohens_d", "log2fc"),
  AdjustMethod = "fdr",
  Alpha = 0.05,
  Format = c("tiered", "classic", "fdr_only", "directional", "effect_gradient",
    "minimal", "neon"),
  LabelMode = c("none", "top_n", "significant", "fdr", "extreme"),
  TopN = 10,
  Relabel = TRUE,
  Codebook = NULL,
  InteractiveLabels = TRUE
)
```

## Arguments

- Data:

  A data frame.

- xVars:

  Character vector of predictor variable names to screen.

- yVar:

  Character string naming the outcome variable.

- Covariates:

  Optional character vector of covariate variable names.

- OutcomeType:

  Outcome type. One of `"auto"`, `"continuous"`, or `"categorical"`. If
  `"auto"`, numeric outcomes are treated as continuous and nonnumeric
  two-level outcomes are treated as categorical.

- EffectMetric:

  Effect metric. One of `"auto"`, `"cohens_d"`, or `"log2fc"`. Used for
  two-group categorical outcomes. For continuous outcomes, the effect is
  always a standardized beta from a model using scaled predictor and
  scaled outcome.

- AdjustMethod:

  Multiple-comparison correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  is `"fdr"`.

- Alpha:

  Significance threshold for raw and adjusted p-values. Default is
  `0.05`.

- Format:

  Color format. One of `"tiered"`, `"classic"`, `"fdr_only"`,
  `"directional"`, `"effect_gradient"`, `"minimal"`, or `"neon"`.

- LabelMode:

  Labeling mode. One of `"none"`, `"top_n"`, `"significant"`, `"fdr"`,
  or `"extreme"`.

- TopN:

  Number of variables to label when `LabelMode` is `"top_n"` or
  `"extreme"`. Default is `10`.

- Relabel:

  Logical. If `TRUE`, variable labels are used when available. Labels
  are pulled first from `Codebook` if supplied, then from variable label
  attributes. Default is `TRUE`.

- Codebook:

  Optional codebook data frame with columns `Variable` and `Label`.

- InteractiveLabels:

  Logical. If `TRUE`, a `text` aesthetic is added for compatibility with
  `plotly::ggplotly(tooltip = "text")`. Default is `TRUE`.

## Value

A named list with `RawPPlot`, `FDRPlot`, and `ResultsTable`. `RawPPlot`
uses `-log10(PValue)` on the y-axis. `FDRPlot` uses `-log10(FDR)` on the
y-axis. `ResultsTable` is a tibble with one row per analyzed predictor.

## Details

Missingness is handled pairwise by model, so each predictor is analyzed
using all complete observations available for that predictor, the
outcome, and any covariates. Invalid predictor variables are removed
with a warning.

## Examples

``` r
PlotVolcanoEffects(
  Data = mtcars,
  xVars = c("disp", "hp", "drat", "wt", "qsec"),
  yVar = "mpg",
  Covariates = "cyl",
  OutcomeType = "continuous",
  LabelMode = "top_n",
  TopN = 3
)
#> Warning: Ignoring unknown aesthetics: text
#> Warning: Ignoring unknown aesthetics: text
#> $RawPPlot

#> 
#> $FDRPlot

#> 
#> $ResultsTable
#> # A tibble: 5 × 18
#>   Variable Label Outcome OutcomeType Effect EffectType        PValue     N Note 
#>   <chr>    <chr> <chr>   <chr>        <dbl> <chr>              <dbl> <int> <chr>
#> 1 disp     disp  mpg     continuous  -0.423 Standardized be… 5.42e-2    32 NA   
#> 2 hp       hp    mpg     continuous  -0.218 Standardized be… 2.13e-1    32 NA   
#> 3 drat     drat  mpg     continuous   0.166 Standardized be… 2.20e-1    32 NA   
#> 4 wt       wt    mpg     continuous  -0.518 Standardized be… 2.22e-4    32 NA   
#> 5 qsec     qsec  mpg     continuous  -0.131 Standardized be… 2.76e-1    32 NA   
#> # ℹ 9 more variables: FDR <dbl>, NegLog10P <dbl>, NegLog10FDR <dbl>,
#> #   Significant <lgl>, FDRSignificant <lgl>, Direction <chr>,
#> #   SignificanceTier <chr>, Tooltip <chr>, PlotLabel <chr>
#> 

PlotVolcanoEffects(
  Data = mtcars,
  xVars = c("mpg", "disp", "hp", "drat", "wt", "qsec"),
  yVar = "am",
  OutcomeType = "categorical",
  EffectMetric = "cohens_d",
  LabelMode = "top_n",
  TopN = 3
)
#> Warning: Ignoring unknown aesthetics: text
#> Warning: Ignoring unknown aesthetics: text
#> $RawPPlot

#> 
#> $FDRPlot

#> 
#> $ResultsTable
#> # A tibble: 6 × 18
#>   Variable Label Outcome OutcomeType Effect EffectType        PValue     N Note 
#>   <chr>    <chr> <chr>   <chr>        <dbl> <chr>              <dbl> <int> <chr>
#> 1 mpg      mpg   am      categorical  1.20  Cohen-style d 0.000285      32 NA   
#> 2 disp     disp  am      categorical -1.18  Cohen-style d 0.000366      32 NA   
#> 3 hp       hp    am      categorical -0.487 Cohen-style d 0.180         32 NA   
#> 4 drat     drat  am      categorical  1.43  Cohen-style d 0.00000473    32 NA   
#> 5 wt       wt    am      categorical -1.39  Cohen-style d 0.0000113     32 NA   
#> 6 qsec     qsec  am      categorical -0.461 Cohen-style d 0.206         32 NA   
#> # ℹ 9 more variables: FDR <dbl>, NegLog10P <dbl>, NegLog10FDR <dbl>,
#> #   Significant <lgl>, FDRSignificant <lgl>, Direction <chr>,
#> #   SignificanceTier <chr>, Tooltip <chr>, PlotLabel <chr>
#> 
```
