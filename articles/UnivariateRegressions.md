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

For this example we will evaluate a panel of biomarkers as predictors of
several clinical outcomes.

``` r

OutcomeVars <- c(
  "Diagnosis",
  "sex"
)

PredictorVars <- c(
  "Genotype",
  "Calbindin",
  "Ferritin",
  "MMP7",
  "Calbindin",
  "Sortilin",
  "Osteopontin"
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

[TABLE]

Estimates and standard errors are combined into a compact format and
significance stars are automatically added.

## 7 Detailed regression table

The larger table provides additional model details.

``` r

UniObj$LargeTable
```

[TABLE]

This version is often useful during exploratory analyses and quality
control.

## 8 Including covariates

Covariates can be included in every regression model.

For example, age is commonly included when evaluating biomarker
associations.

``` r

UniObj_Covar <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars,
  covariates = "age"
)
```

``` r

UniObj_Covar$FormattedTable
```

[TABLE]

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

[TABLE]

Standardized coefficients represent changes in standard deviation units
and facilitate comparison across predictors.

## 10 Accessing model objects

The fitted regression models are returned in the output object.

``` r

names(
  UniObj$ModelSummaries
)
```

    [1] "Diagnosis" "sex"      

For example:

``` r

summary(
  UniObj$ModelSummaries$Diagnosis$Ferritin
)
```


    Call:
    stats::glm(formula = f, family = stats::binomial(), data = ModelData)

    Coefficients:
                Estimate Std. Error z value Pr(>|z|)
    (Intercept)  -1.9184     0.4720  -4.065 4.81e-05 ***
    Ferritin      0.3359     0.1601   2.098   0.0359 *
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    (Dispersion parameter for binomial family taken to be 1)

        Null deviance: 390.60  on 332  degrees of freedom
    Residual deviance: 386.11  on 331  degrees of freedom
    AIC: 390.11

    Number of Fisher Scoring iterations: 4

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
    [1] dplyr_1.2.1            SciDataReportR_20.10.0

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6           xfun_0.59              bayestestR_0.18.1
     [4] ggplot2_4.0.3          insight_1.5.2          rstatix_1.0.0
     [7] lattice_0.22-9         paletteer_1.7.0        vctrs_0.7.3
    [10] tools_4.6.1            generics_0.1.4         datawizard_1.3.1
    [13] tibble_3.3.1           pkgconfig_2.0.3        RColorBrewer_1.1-3
    [16] correlation_0.8.8      S7_0.2.2               gt_1.3.0
    [19] RcppParallel_5.1.11-2  lifecycle_1.0.5        compiler_4.6.1
    [22] farver_2.1.2           stringr_1.6.0          carData_3.0-6
    [25] snakecase_0.11.1       litedown_0.9           sass_0.4.10
    [28] htmltools_0.5.9        yaml_2.3.12            Formula_1.2-5
    [31] pillar_1.11.1          car_3.1-5              tidyr_1.3.2
    [34] broom.helpers_1.22.0   statsExpressions_2.0.0 abind_1.4-8
    [37] commonmark_2.0.0       tidyselect_1.2.1       sjlabelled_1.2.0
    [40] digest_0.6.39          mvtnorm_1.4-1          stringi_1.8.7
    [43] gtsummary_2.5.1        purrr_1.2.2            rematch2_2.1.2
    [46] labeling_0.4.3         forcats_1.0.1          ggstatsplot_1.0.0
    [49] labelled_2.16.0        fastmap_1.2.0          grid_4.6.1
    [52] cli_3.6.6              magrittr_2.0.5         base64enc_0.1-6
    [55] cards_0.8.0            patchwork_1.3.2        dichromat_2.0-0.1
    [58] broom_1.0.13           withr_3.0.3            scales_1.4.0
    [61] backports_1.5.1        estimability_2.0.0     rmarkdown_2.31
    [64] emmeans_2.0.3          otel_0.2.0             hms_1.1.4
    [67] coda_0.19-4.1          evaluate_1.0.5         knitr_1.51
    [70] haven_2.5.5            parameters_0.29.2      markdown_2.0
    [73] rstantools_2.6.0       rlang_1.2.0            xtable_1.8-8
    [76] glue_1.8.1             xml2_1.6.0             jsonlite_2.0.0
    [79] effectsize_1.0.2       R6_2.6.1               fs_2.1.0              
