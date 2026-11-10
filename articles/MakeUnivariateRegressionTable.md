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

`FormattedTable` is the report-facing wide table: each outcome has its
own column spanner with an estimate (95% CI) and p-value. Significant
effects carry stars and are bolded.

``` r

Reg_Obj_Un$FormattedTable
```

[TABLE]

`LargeTable` uses the same wide outcome layout but keeps N, estimate,
standard error, confidence limits, and p-value as separate columns for
QC.

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

    # A tibble: 8 × 21
      Outcome   OutcomeLabel OutcomeFamily EffectType Predictor PredictorLabel Term
      <chr>     <chr>        <chr>         <chr>      <chr>     <chr>          <chr>
    1 Calbindin Calbindin    linear        Estimate   Diagnosis Diagnosis      Diag…
    2 Calbindin Calbindin    linear        Estimate   age       Age            age
    3 Ferritin  Ferritin     linear        Estimate   Diagnosis Diagnosis      Diag…
    4 Ferritin  Ferritin     linear        Estimate   age       Age            age
    5 MMP7      Matrix meta… linear        Estimate   Diagnosis Diagnosis      Diag…
    6 MMP7      Matrix meta… linear        Estimate   age       Age            age
    7 Sortilin  Sortilin     linear        Estimate   Diagnosis Diagnosis      Diag…
    8 Sortilin  Sortilin     linear        Estimate   age       Age            age
    # ℹ 14 more variables: Level <chr>, TermLabel <chr>, N <int>, Estimate <dbl>,
    #   StdError <dbl>, ConfLow <dbl>, ConfHigh <dbl>, PValue <dbl>,
    #   Significant <lgl>, ReferenceValue <dbl>, PValueFDR <dbl>, Group <chr>,
    #   Status <chr>, Reason <chr>

Because it is an ordinary dataframe, you can filter, sort, or export it
directly, and pass it (or a subset of it) straight to
[`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md).

``` r

Reg_Obj_Un$Results %>%
  filter(Significant) %>%
  arrange(PValue)
```

    # A tibble: 3 × 21
      Outcome  OutcomeLabel  OutcomeFamily EffectType Predictor PredictorLabel Term
      <chr>    <chr>         <chr>         <chr>      <chr>     <chr>          <chr>
    1 MMP7     Matrix metal… linear        Estimate   Diagnosis Diagnosis      Diag…
    2 Sortilin Sortilin      linear        Estimate   Diagnosis Diagnosis      Diag…
    3 Ferritin Ferritin      linear        Estimate   Diagnosis Diagnosis      Diag…
    # ℹ 14 more variables: Level <chr>, TermLabel <chr>, N <int>, Estimate <dbl>,
    #   StdError <dbl>, ConfLow <dbl>, ConfHigh <dbl>, PValue <dbl>,
    #   Significant <lgl>, ReferenceValue <dbl>, PValueFDR <dbl>, Group <chr>,
    #   Status <chr>, Reason <chr>

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

    NULL

``` r

names(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin)
```

    NULL

``` r

summary(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin$Ferritin)
```

    Length  Class   Mode
         0   NULL   NULL 

## 11 Reproducibility

``` r

# save.image("MakeUnivariateRegressionTable_workspace.RData")
print(sessionInfo())
```

    R version 4.6.2 (2026-10-29)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 26.04.1 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0

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
    [1] dplyr_1.2.1            SciDataReportR_21.21.0

    loaded via a namespace (and not attached):
     [1] gt_1.4.0               utf8_1.2.6             sass_0.4.10
     [4] generics_0.1.4         tidyr_1.3.2            xml2_1.6.0
     [7] rstatix_1.1.0          gtsummary_2.6.1        lattice_0.23-1
    [10] hms_1.1.4              digest_0.6.39          magrittr_2.0.5
    [13] evaluate_1.0.5         grid_4.6.2             estimability_2.0.0
    [16] RColorBrewer_1.1-3     mvtnorm_1.4-2          fastmap_1.2.0
    [19] jsonlite_2.0.0         backports_1.5.1        Formula_1.2-6
    [22] rematch2_2.1.2         purrr_1.2.2            scales_1.4.0
    [25] abind_1.4-8            cli_3.6.6              labelled_2.16.1
    [28] rlang_1.3.0            withr_3.0.3            yaml_2.3.12
    [31] otel_0.3.0             statsExpressions_2.2.1 tools_4.6.2
    [34] datawizard_1.5.0       rstantools_2.7.1       ggstatsplot_1.1.4
    [37] coda_0.19-4.1          ggplot2_4.0.3          sjlabelled_1.2.0
    [40] bayestestR_0.20.0      forcats_1.0.1          broom_1.0.13
    [43] paletteer_1.7.0        vctrs_0.7.3            R6_2.6.1
    [46] lifecycle_1.0.5        emmeans_2.0.4          snakecase_0.11.1
    [49] fs_2.1.0               car_3.1-5              insight_1.7.0
    [52] pkgconfig_2.0.3        RcppParallel_6.2.1     pillar_1.11.1
    [55] gtable_0.3.6           glue_1.8.1             haven_2.5.5
    [58] xfun_0.62              tibble_3.3.1           tidyselect_1.2.1
    [61] parameters_0.30.0      knitr_1.52             dichromat_2.0-1
    [64] farver_2.1.2           xtable_1.8-8           htmltools_0.5.9
    [67] carData_3.0-6          rmarkdown_2.32         compiler_4.6.2
    [70] S7_0.2.2              
