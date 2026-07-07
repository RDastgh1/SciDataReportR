# Screening Outcomes and Predictors with MakeUnivariateRegressionTable()

## 1 Overview

[`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
screens outcome-predictor relationships and returns a report-facing
table, a detail table, a tidy results dataframe, fitted models, and
metadata.

Use it when you want to quickly ask:

- Which predictors are associated with one or more outcomes?
- Do associations persist after covariate adjustment?
- Are the displayed labels and logistic event levels what I think they
  are?

The function chooses the model family automatically:

- Numeric outcomes use linear regression.
- Two-level factor, character, or logical outcomes use logistic
  regression.
- Logistic estimates are exponentiated by default and reported as odds
  ratios.

> **Note**
>
> [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
> was previously named
> [`UnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md).
> The old name still works as a backwards-compatible synonym, so
> existing scripts continue to run.

## 2 Load packages

``` r

library(SciDataReportR)
library(dplyr)
```

## 3 Load example data

The examples use `SampleData` and `SampleVariableTypes`, then apply
labels and recoding with
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md).

``` r

data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData
```

## 4 Basic linear screening

This example screens several labelled biomarkers against two predictors.

``` r

vars_Outcomes <- c(
  "Calbindin",
  "Ferritin",
  "MMP7",
  "Sortilin"
)

vars_Predictors <- c(
  "Diagnosis",
  "age"
)
```

``` r

Reg_Obj_Un <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors,
  Standardize = TRUE
)
```

## 5 Formatted and detailed tables

`FormattedTable` is the report-facing table.

``` r

Reg_Obj_Un$FormattedTable
```

[TABLE]

`LargeTable` keeps more of the model detail and is useful for QC.

``` r

Reg_Obj_Un$LargeTable
```

[TABLE]

## 6 Tidy results dataframe

`Results` holds the underlying numbers as a plain dataframe: one row per
estimated term, with the estimate, standard error, confidence interval,
p-value, and the labels used in the tables.

``` r

Reg_Obj_Un$Results
```

        Outcome               OutcomeLabel OutcomeFamily EffectType Predictor
    1 Calbindin                  Calbindin        linear   Estimate Diagnosis
    2 Calbindin                  Calbindin        linear   Estimate       age
    3  Ferritin                   Ferritin        linear   Estimate Diagnosis
    4  Ferritin                   Ferritin        linear   Estimate       age
    5      MMP7 Matrix metalloproteinase 7        linear   Estimate Diagnosis
    6      MMP7 Matrix metalloproteinase 7        linear   Estimate       age
    7  Sortilin                   Sortilin        linear   Estimate Diagnosis
    8  Sortilin                   Sortilin        linear   Estimate       age
      PredictorLabel              Term    Level            TermLabel   N
    1      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333
    2            Age               age     <NA>                  Age 322
    3      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333
    4            Age               age     <NA>                  Age 322
    5      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333
    6            Age               age     <NA>                  Age 322
    7      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333
    8            Age               age     <NA>                  Age 322
         Estimate   StdError     ConfLow   ConfHigh       PValue Significant
    1  0.22318256 0.12254159 -0.01787596 0.46424108 6.946687e-02       FALSE
    2  0.05166760 0.05580595 -0.05812530 0.16146050 3.552248e-01       FALSE
    3  0.25934156 0.12232632  0.01870650 0.49997662 3.474297e-02        TRUE
    4  0.03027304 0.05597431 -0.07985109 0.14039718 5.889954e-01       FALSE
    5  0.55090515 0.11937344  0.31607888 0.78573143 5.630656e-06        TRUE
    6 -0.04643922 0.05571349 -0.15605022 0.06317177 4.051640e-01       FALSE
    7  0.41966369 0.12097458  0.18168773 0.65763966 5.914370e-04        TRUE
    8  0.06679376 0.05597311 -0.04332801 0.17691553 2.336285e-01       FALSE
      ReferenceValue
    1              0
    2              0
    3              0
    4              0
    5              0
    6              0
    7              0
    8              0

Because it is an ordinary dataframe, you can filter, sort, or export it
directly, and pass it (or a subset of it) straight to
[`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).

``` r

Reg_Obj_Un$Results %>%
  filter(Significant) %>%
  arrange(PValue)
```

       Outcome               OutcomeLabel OutcomeFamily EffectType Predictor
    1     MMP7 Matrix metalloproteinase 7        linear   Estimate Diagnosis
    2 Sortilin                   Sortilin        linear   Estimate Diagnosis
    3 Ferritin                   Ferritin        linear   Estimate Diagnosis
      PredictorLabel              Term    Level            TermLabel   N  Estimate
    1      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333 0.5509052
    2      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333 0.4196637
    3      Diagnosis DiagnosisImpaired Impaired Diagnosis : Impaired 333 0.2593416
       StdError   ConfLow  ConfHigh       PValue Significant ReferenceValue
    1 0.1193734 0.3160789 0.7857314 5.630656e-06        TRUE              0
    2 0.1209746 0.1816877 0.6576397 5.914370e-04        TRUE              0
    3 0.1223263 0.0187065 0.4999766 3.474297e-02        TRUE              0

## 7 Covariate adjustment

Use `covariates` when each screen should adjust for the same covariates.

``` r

Reg_Obj_Un_Covar <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors,
  covariates = "sex",
  Standardize = TRUE
)

Reg_Obj_Un_Covar$FormattedTable
```

[TABLE]

## 8 Many predictors for one outcome

The same function also works when the question is one outcome against
many candidate predictors.

``` r

Reg_Obj_ManyPredictors <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = "Sortilin",
  predictor_vars = c("Diagnosis", "age", "Calbindin", "Ferritin", "MMP7"),
  Standardize = TRUE
)

Reg_Obj_ManyPredictors$FormattedTable
```

[TABLE]

## 9 Logistic outcomes

Binary outcomes use logistic regression automatically. By default,
estimates are exponentiated and shown as odds ratios.

``` r

Reg_Obj_Logistic <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = "Diagnosis",
  predictor_vars = c("age", "Calbindin", "Ferritin", "MMP7"),
  Standardize = TRUE
)

Reg_Obj_Logistic$FormattedTable
```

[TABLE]

Always check the metadata for logistic models. It records the reference
level and event level.

``` r

Reg_Obj_Logistic$Metadata$Outcomes
```

        Outcome OutcomeLabel OutcomeFamily ReferenceLevel EventLevel
    1 Diagnosis    Diagnosis      logistic        Control   Impaired

## 10 Access fitted models

Fitted models are returned by outcome, then predictor.

``` r

names(Reg_Obj_Un$ModelSummaries)
```

    [1] "Calbindin" "Ferritin"  "MMP7"      "Sortilin" 

``` r

names(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin)
```

    [1] "Diagnosis" "age"       "Calbindin" "Ferritin"  "MMP7"     

``` r

summary(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin$Ferritin)
```


    Call:
    stats::lm(formula = f, data = ModelData)

    Residuals:
         Min       1Q   Median       3Q      Max
    -1.90775 -0.50046 -0.04898  0.50987  2.18216

    Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
    (Intercept) 1.396e-17  4.318e-02    0.00        1
    Ferritin    6.173e-01  4.324e-02   14.28   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 0.7879 on 331 degrees of freedom
    Multiple R-squared:  0.3811,    Adjusted R-squared:  0.3792
    F-statistic: 203.8 on 1 and 331 DF,  p-value: < 2.2e-16

## 11 Reproducibility

``` r

# save.image("MakeUnivariateRegressionTable_workspace.RData")
print(sessionInfo())
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
    [1] dplyr_1.2.1            SciDataReportR_20.11.0

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
    [40] digest_0.6.39          stringi_1.8.7          mvtnorm_1.4-1
    [43] gtsummary_2.5.1        purrr_1.2.2            rematch2_2.1.2
    [46] forcats_1.0.1          ggstatsplot_1.0.0      labelled_2.16.0
    [49] fastmap_1.2.0          grid_4.6.1             cli_3.6.6
    [52] magrittr_2.0.5         base64enc_0.1-6        cards_0.8.1
    [55] patchwork_1.3.2        dichromat_2.0-0.1      broom_1.0.13
    [58] withr_3.0.3            scales_1.4.0           backports_1.5.1
    [61] estimability_2.0.0     rmarkdown_2.31         emmeans_2.0.3
    [64] otel_0.2.0             hms_1.1.4              coda_0.19-4.1
    [67] evaluate_1.0.5         knitr_1.51             haven_2.5.5
    [70] parameters_0.29.2      markdown_2.0           rstantools_2.6.0
    [73] rlang_1.3.0            xtable_1.8-8           glue_1.8.1
    [76] xml2_1.6.0             jsonlite_2.0.0         effectsize_1.0.2
    [79] R6_2.6.1               fs_2.1.0              
