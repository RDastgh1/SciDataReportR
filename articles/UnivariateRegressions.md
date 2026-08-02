# Screening predictors with MakeUnivariateRegressionTable()

## 1 Overview

[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
performs large numbers of regression analyses simultaneously and returns
publication-ready summary tables.

This approach is particularly useful for:

- Biomarker screening
- Variable prioritization
- Covariate discovery
- Exploratory analyses
- Hypothesis generation

This vignette demonstrates:

- Running multiple univariate regressions
- Publication-ready tables
- Detailed regression tables
- Covariate adjustment
- Standardized coefficients
- Forest plot visualization

## 2 Load packages

``` r

library(SciDataReportR)
library(dplyr)
```

## 3 Load example data

``` r

data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData
```

## 4 Define outcomes and predictors

For this example we will evaluate clinical and demographic predictors of
continuous biomarker outcomes.

``` r

OutcomeVars <- c(
  "Calbindin",
  "Ferritin"
)

PredictorVars <- c(
  "Diagnosis",
  "age",
  "sex",
  "Genotype"
)
```

## 5 Create regression tables

``` r

UniObj <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars
)
```

## 6 Publication-ready table

The formatted table is designed for reporting and manuscript
preparation.

``` r

UniObj$FormattedTable
```

| Variable             | Effect   | N   | Estimate (95% CI)          | p-value |
|----------------------|----------|-----|----------------------------|---------|
| Calbindin            |          |     |                            |         |
| Diagnosis : Impaired | Estimate | 333 | 0.936 (-0.075, 1.95)       | 0.069   |
| Age                  | Estimate | 322 | 0.0164 (-0.0185, 0.0513)   | 0.36    |
| Sex : Male           | Estimate | 333 | -1.62 (-2.53, -0.704)      | \<0.001 |
| Genotype : E2E3      | Estimate | 333 | 1.74 (-4.26, 7.75)         | 0.57    |
| Genotype : E2E4      | Estimate | 333 | 3.13 (-3.42, 9.67)         | 0.35    |
| Genotype : E3E3      | Estimate | 333 | 1.81 (-4.08, 7.69)         | 0.55    |
| Genotype : E3E4      | Estimate | 333 | 2.46 (-3.44, 8.37)         | 0.41    |
| Genotype : E4E4      | Estimate | 333 | 1.31 (-4.98, 7.6)          | 0.68    |
| Ferritin             |          |     |                            |         |
| Diagnosis : Impaired | Estimate | 333 | 0.203 (0.0147, 0.392)      | 0.035   |
| Age                  | Estimate | 322 | 0.0018 (-0.00474, 0.00834) | 0.59    |
| Sex : Male           | Estimate | 333 | 0.253 (0.0815, 0.425)      | 0.004   |
| Genotype : E2E3      | Estimate | 333 | 0.218 (-0.907, 1.34)       | 0.7     |
| Genotype : E2E4      | Estimate | 333 | 0.223 ( -1, 1.45)          | 0.72    |
| Genotype : E3E3      | Estimate | 333 | 0.195 (-0.907, 1.3)        | 0.73    |
| Genotype : E3E4      | Estimate | 333 | 0.267 (-0.838, 1.37)       | 0.63    |
| Genotype : E4E4      | Estimate | 333 | -0.059 (-1.24, 1.12)       | 0.92    |

Estimates and standard errors are combined into a compact format and
significance stars are automatically added.

## 7 Detailed regression table

The larger table provides additional model details.

``` r

UniObj$LargeTable
```

| Outcome | Variable | Effect | N | Estimate | SE | 95% CI Low | 95% CI High | p-value |
|----|----|----|----|----|----|----|----|----|
| Calbindin |  |  |  |  |  |  |  |  |
| Calbindin | Diagnosis : Impaired | Estimate | 333 | 0.936 | 0.514 | −0.075 | 1.948 | 0.069 |
| Calbindin | Age | Estimate | 322 | 0.016 | 0.018 | −0.018 | 0.051 | 0.355 |
| Calbindin | Sex : Male | Estimate | 333 | −1.617 | 0.464 | −2.530 | −0.704 | 0.001 |
| Calbindin | Genotype : E2E3 | Estimate | 333 | 1.745 | 3.054 | −4.263 | 7.753 | 0.568 |
| Calbindin | Genotype : E2E4 | Estimate | 333 | 3.127 | 3.326 | −3.415 | 9.670 | 0.348 |
| Calbindin | Genotype : E3E3 | Estimate | 333 | 1.807 | 2.992 | −4.080 | 7.694 | 0.546 |
| Calbindin | Genotype : E3E4 | Estimate | 333 | 2.464 | 3.003 | −3.443 | 8.371 | 0.412 |
| Calbindin | Genotype : E4E4 | Estimate | 333 | 1.310 | 3.195 | −4.976 | 7.596 | 0.682 |
| Ferritin |  |  |  |  |  |  |  |  |
| Ferritin | Diagnosis : Impaired | Estimate | 333 | 0.203 | 0.096 | 0.015 | 0.392 | 0.035 |
| Ferritin | Age | Estimate | 322 | 0.002 | 0.003 | −0.005 | 0.008 | 0.589 |
| Ferritin | Sex : Male | Estimate | 333 | 0.253 | 0.087 | 0.082 | 0.425 | 0.004 |
| Ferritin | Genotype : E2E3 | Estimate | 333 | 0.218 | 0.572 | −0.907 | 1.343 | 0.703 |
| Ferritin | Genotype : E2E4 | Estimate | 333 | 0.223 | 0.623 | −1.001 | 1.448 | 0.720 |
| Ferritin | Genotype : E3E3 | Estimate | 333 | 0.195 | 0.560 | −0.907 | 1.297 | 0.728 |
| Ferritin | Genotype : E3E4 | Estimate | 333 | 0.267 | 0.562 | −0.838 | 1.373 | 0.635 |
| Ferritin | Genotype : E4E4 | Estimate | 333 | −0.059 | 0.598 | −1.236 | 1.118 | 0.921 |

This version is often useful during exploratory analyses and quality
control.

## 8 Including covariates

Covariates can be included in every regression model.

For example, MMP7 can be included when evaluating whether associations
are independent of another biomarker.

``` r

UniObj_Covar <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars,
  covariates = "MMP7"
)
```

``` r

UniObj_Covar$FormattedTable
```

| Variable             | Effect   | N   | Estimate (95% CI)          | p-value |
|----------------------|----------|-----|----------------------------|---------|
| Calbindin            |          |     |                            |         |
| Diagnosis : Impaired | Estimate | 333 | 1.05 (0.00979, 2.1)        | 0.048   |
| Age                  | Estimate | 322 | 0.0161 (-0.0189, 0.051)    | 0.37    |
| Sex : Male           | Estimate | 333 | -1.62 (-2.54, -0.696)      | \<0.001 |
| Genotype : E2E3      | Estimate | 333 | 1.73 (-4.29, 7.74)         | 0.57    |
| Genotype : E2E4      | Estimate | 333 | 3.09 (-3.46, 9.65)         | 0.35    |
| Genotype : E3E3      | Estimate | 333 | 1.79 (-4.11, 7.68)         | 0.55    |
| Genotype : E3E4      | Estimate | 333 | 2.44 (-3.47, 8.36)         | 0.42    |
| Genotype : E4E4      | Estimate | 333 | 1.28 (-5.01, 7.58)         | 0.69    |
| Ferritin             |          |     |                            |         |
| Diagnosis : Impaired | Estimate | 333 | 0.144 (-0.0493, 0.337)     | 0.14    |
| Age                  | Estimate | 322 | 0.00222 (-0.00425, 0.0087) | 0.5     |
| Sex : Male           | Estimate | 333 | 0.223 (0.0513, 0.395)      | 0.011   |
| Genotype : E2E3      | Estimate | 333 | 0.243 (-0.87, 1.36)        | 0.67    |
| Genotype : E2E4      | Estimate | 333 | 0.272 (-0.94, 1.48)        | 0.66    |
| Genotype : E3E3      | Estimate | 333 | 0.221 (-0.869, 1.31)       | 0.69    |
| Genotype : E3E4      | Estimate | 333 | 0.294 (-0.8, 1.39)         | 0.6     |
| Genotype : E4E4      | Estimate | 333 | -0.02 (-1.18, 1.14)        | 0.97    |

This allows investigators to determine whether associations remain
significant after accounting for potential confounding factors.

## 9 Standardized coefficients

When predictors are measured on different scales, standardized
coefficients can simplify interpretation.

``` r

UniObj_Std <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars,
  Standardize = TRUE
)
```

``` r

UniObj_Std$FormattedTable
```

| Variable             | Effect   | N   | Estimate (95% CI)       | p-value |
|----------------------|----------|-----|-------------------------|---------|
| Calbindin            |          |     |                         |         |
| Diagnosis : Impaired | Estimate | 333 | 0.223 (-0.0179, 0.464)  | 0.069   |
| Age                  | Estimate | 322 | 0.0517 (-0.0581, 0.161) | 0.36    |
| Sex : Male           | Estimate | 333 | -0.385 (-0.603, -0.168) | \<0.001 |
| Genotype : E2E3      | Estimate | 333 | 0.416 (-1.02, 1.85)     | 0.57    |
| Genotype : E2E4      | Estimate | 333 | 0.745 (-0.814, 2.31)    | 0.35    |
| Genotype : E3E3      | Estimate | 333 | 0.431 (-0.973, 1.83)    | 0.55    |
| Genotype : E3E4      | Estimate | 333 | 0.587 (-0.821, 2)       | 0.41    |
| Genotype : E4E4      | Estimate | 333 | 0.312 (-1.19, 1.81)     | 0.68    |
| Ferritin             |          |     |                         |         |
| Diagnosis : Impaired | Estimate | 333 | 0.259 (0.0187, 0.5)     | 0.035   |
| Age                  | Estimate | 322 | 0.0303 (-0.0799, 0.14)  | 0.59    |
| Sex : Male           | Estimate | 333 | 0.323 (0.104, 0.542)    | 0.004   |
| Genotype : E2E3      | Estimate | 333 | 0.278 (-1.16, 1.71)     | 0.7     |
| Genotype : E2E4      | Estimate | 333 | 0.285 (-1.28, 1.85)     | 0.72    |
| Genotype : E3E3      | Estimate | 333 | 0.248 (-1.16, 1.65)     | 0.73    |
| Genotype : E3E4      | Estimate | 333 | 0.341 (-1.07, 1.75)     | 0.63    |
| Genotype : E4E4      | Estimate | 333 | -0.0753 (-1.58, 1.43)   | 0.92    |

Standardized coefficients represent changes in standard deviation units
and facilitate comparison across predictors.

## 10 Accessing model objects

The fitted regression models are returned in the output object.

``` r

names(
  UniObj$ModelSummaries
)
```

    NULL

For example:

``` r

summary(
  UniObj$ModelSummaries$Calbindin$Diagnosis
)
```

    Length  Class   Mode
         0   NULL   NULL 

This allows additional diagnostics and custom analyses.

## 11 Creating a forest plot

A forest plot provides a compact visual summary of regression results.

``` r

ForestPlot <- PlotForestFromTable(
  UniObj
)

ForestPlot
```

![](UnivariateRegressions_files/figure-html/unnamed-chunk-13-1.png)

Forest plots make it easy to identify:

- Strong predictors
- Significant predictors
- Direction of effects
- Precision of estimates

## 12 Forest plot with covariate adjustment

The same visualization can be generated using adjusted models.

``` r

PlotForestFromTable(
  UniObj_Covar
)
```

![](UnivariateRegressions_files/figure-html/unnamed-chunk-14-1.png)

Comparing adjusted and unadjusted models can help identify associations
that may be explained by confounding variables.

## 13 Recommended workflow

A common workflow is:

1.  Use
    [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
    to screen large numbers of predictors.
2.  Review the formatted table.
3.  Compare standardized and unstandardized results.
4.  Add important covariates.
5.  Visualize findings using
    [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).
6.  Follow up significant findings using multivariable models.

## 14 Summary

[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
provides a rapid screening framework for evaluating many predictors
across multiple outcomes.

Key features include:

- Multiple outcomes
- Multiple predictors
- Optional covariate adjustment
- Standardized coefficients
- Publication-ready tables
- Detailed regression output
- Forest plot visualization

## 15 Related functions

- [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md)
- [`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
- [`PlotMiningMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotMiningMatrix.md)
- [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
- `CreatePCA()`

## 16 Session information

``` r

sessionInfo()
```

    R version 4.6.1 (2026-06-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.4 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] dplyr_1.2.1            SciDataReportR_20.19.0

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6           xfun_0.60              bayestestR_0.18.1
     [4] ggplot2_4.0.3          insight_1.5.2          rstatix_1.1.0
     [7] lattice_0.22-9         paletteer_1.7.0        vctrs_0.7.3
    [10] tools_4.6.1            generics_0.1.4         datawizard_1.3.1
    [13] tibble_3.3.1           pkgconfig_2.0.3        RColorBrewer_1.1-3
    [16] correlation_0.8.8      S7_0.2.2               RcppParallel_6.2.0
    [19] gt_1.3.0               lifecycle_1.0.5        compiler_4.6.1
    [22] farver_2.1.2           carData_3.0-6          snakecase_0.11.1
    [25] sass_0.4.10            htmltools_0.5.9        yaml_2.3.12
    [28] Formula_1.2-5          pillar_1.11.1          car_3.1-5
    [31] tidyr_1.3.2            statsExpressions_2.0.0 abind_1.4-8
    [34] tidyselect_1.2.1       sjlabelled_1.2.0       digest_0.6.39
    [37] mvtnorm_1.4-2          gtsummary_2.5.1        purrr_1.2.2
    [40] rematch2_2.1.2         labeling_0.4.3         forcats_1.0.1
    [43] ggstatsplot_1.0.0      labelled_2.16.0        fastmap_1.2.0
    [46] grid_4.6.1             cli_3.6.6              magrittr_2.0.5
    [49] patchwork_1.3.2        dichromat_2.0-1        broom_1.0.13
    [52] withr_3.0.3            scales_1.4.0           backports_1.5.1
    [55] estimability_2.0.0     rmarkdown_2.31         emmeans_2.0.4
    [58] otel_0.2.0             hms_1.1.4              coda_0.19-4.1
    [61] evaluate_1.0.5         knitr_1.51             haven_2.5.5
    [64] parameters_0.29.2      rstantools_2.7.0       rlang_1.3.0
    [67] xtable_1.8-8           glue_1.8.1             xml2_1.6.0
    [70] jsonlite_2.0.0         effectsize_1.0.3       R6_2.6.1
    [73] fs_2.1.0              
