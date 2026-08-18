# Creating comparison tables with MakeComparisonTable()

## 1 Overview

[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
creates publication-ready summary tables for comparing groups across
continuous and categorical variables.

This vignette demonstrates:

- Basic group comparisons
- Adding effect sizes
- Parametric versus nonparametric testing
- Covariate adjustment
- Pairwise comparisons
- Reference-group comparisons
- Automatic use of variable labels

## 2 Load packages

``` r

library(SciDataReportR)
library(dplyr)
```

## 3 Load example data

SciDataReportR includes an example dataset and codebook.

``` r

data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
) 
df_Revalued <- RevaluedObj$RevaluedData
```

## 4 Basic comparison table

The most common use case is comparing groups across a set of variables.

Here we compare diagnostic groups.

``` r

tbl_basic <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  )
)

tbl_basic
```

[TABLE]

Comparison by Diagnosis (values: mean (SD)). {.table .gt_table
.caption-top quarto-bootstrap="false"}

The resulting table contains:

- Descriptive statistics
- Group comparison tests
- P-values
- Automatically applied variable labels

## 5 Adding effect sizes

Effect sizes provide information about the magnitude of group
differences.

``` r

tbl_effectsize <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddEffectSize = TRUE
)

tbl_effectsize
```

[TABLE]

Comparison by Diagnosis (values: mean (SD)). {.table .gt_table
.caption-top style="width:100%;" quarto-bootstrap="false"}

When reporting group differences, effect sizes should generally be
interpreted alongside p-values.

For parametric continuous outcomes, the reported effect size depends on
the number of groups:

- Two groups use absolute Cohen’s (d). With covariates, this is the
  adjusted estimated marginal mean difference divided by the ANCOVA
  residual standard deviation.
- More than two groups use Cohen’s (f). With covariates, partial
  Cohen’s (f) is calculated from the Type II ANCOVA group effect.

Cohen’s (d) and Cohen’s (f) are not numerically equivalent. In a
balanced two-group design only, (d 2f). Interpret each effect-size scale
in its methodological context; conventional magnitude guides are
heuristics, not thresholds for clinical importance. The generated table
caption and methods footnote remain focused on the summaries and
hypothesis tests used.

## 6 Nonparametric analyses

Nonparametric testing uses rank-based methods that are more robust to
non-normal distributions.

``` r

tbl_nonparametric <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  Parametric = FALSE,
  AddEffectSize = TRUE
)

tbl_nonparametric
```

[TABLE]

Comparison by Diagnosis (values: median \[IQR\]). {.table .gt_table
.caption-top style="width:100%;" quarto-bootstrap="false"}

## 7 Including Covariates

Age is frequently included as a covariate in biomedical analyses.

The example below evaluates group differences after accounting for age.

``` r

tbl_covariate <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  Covariates = "age",
  AddEffectSize = TRUE
)

tbl_covariate
```

[TABLE]

Comparison by Diagnosis (values: mean (SD)). p-values adjusted for Age.
{.table .gt_table .caption-top style="width:100%;"
quarto-bootstrap="false"}

Covariate adjustment can help determine whether observed group
differences persist after controlling for potential confounding factors.

## 8 Pairwise comparisons

When comparing more than two groups, pairwise comparisons can identify
which groups differ from one another.

``` r

tbl_pairwise <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Genotype",
  Variables = c(
    "Diagnosis",
    "age",
    "sex",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddPairwise = TRUE,
  AddEffectSize = TRUE
)

tbl_pairwise
```

[TABLE]

Comparison by Genotype (values: mean (SD)). Pairwise p-values are
bonferroni-adjusted within outcome. {.table .gt_table .caption-top
style="width:100%;" quarto-bootstrap="false"}

You can also choose a referent and compare only to that

``` r

tbl_pairwise_referent <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Genotype",
  Variables = c(
    "Diagnosis",
    "age",
    "sex",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddPairwise = TRUE,
  AddEffectSize = TRUE,
  Referent = "E3E3"
)

tbl_pairwise_referent
```

[TABLE]

Comparison by Genotype (values: mean (SD)). Pairwise p-values are
bonferroni-adjusted within outcome. {.table .gt_table .caption-top
quarto-bootstrap="false"}

## 9 Summary

[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
supports:

- Descriptive statistics
- Parametric analyses
- Nonparametric analyses
- Effect sizes
- Covariate adjustment
- Pairwise testing
- Reference-group comparisons
- Automatic variable labeling

## 10 Related functions

- [`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
- `CreateCodebook()`
- `PlotBoxPlot()`
- [`PlotVolcanoEffects()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotVolcanoEffects.md)
- [`CompareDatasets()`](https://rdastgh1.github.io/SciDataReportR/reference/CompareDatasets.md)

## 11 Session information

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
    [1] dplyr_1.2.1            SciDataReportR_20.24.0

    loaded via a namespace (and not attached):
      [1] Exact_3.3              ggstatsplot_1.0.0      sjlabelled_1.2.0
      [4] tidyselect_1.2.1       rootSolve_1.8.2.4      farver_2.1.2
      [7] statsExpressions_2.0.0 S7_0.2.2               fastmap_1.2.0
     [10] bayestestR_0.18.1      labelled_2.16.0        digest_0.6.39
     [13] estimability_2.0.0     lifecycle_1.0.5        lmom_3.3
     [16] magrittr_2.0.5         compiler_4.6.1         rlang_1.3.0
     [19] sass_0.4.10            tools_4.6.1            yaml_2.3.12
     [22] gt_1.3.0               data.table_1.18.4      knitr_1.51
     [25] xml2_1.6.0             RColorBrewer_1.1-3     abind_1.4-8
     [28] expm_1.0-0             withr_3.0.3            purrr_1.2.2
     [31] nnet_7.3-20            grid_4.6.1             datawizard_1.3.1
     [34] xtable_1.8-8           e1071_1.7-17           gtsummary_2.5.1
     [37] paletteer_1.7.0        ggplot2_4.0.3          emmeans_2.0.4
     [40] scales_1.4.0           MASS_7.3-65            dichromat_2.0-1
     [43] insight_1.5.2          cli_3.6.6              mvtnorm_1.4-2
     [46] rmarkdown_2.31         generics_0.1.4         otel_0.2.0
     [49] RcppParallel_6.2.0     rstudioapi_0.19.0      httr_1.4.8
     [52] tzdb_0.5.0             parameters_0.29.2      commonmark_2.0.0
     [55] readxl_1.5.0           gld_2.6.8              proxy_0.4-29
     [58] effectsize_1.0.3       cellranger_1.1.0       base64enc_0.1-6
     [61] vctrs_0.7.3            Matrix_1.7-5           boot_1.3-32
     [64] sandwich_3.1-3         jsonlite_2.0.0         carData_3.0-6
     [67] car_3.1-5              litedown_0.10          hms_1.1.4
     [70] patchwork_1.3.2        rstatix_1.1.0          Formula_1.2-6
     [73] correlation_0.8.8      tidyr_1.3.2            glue_1.8.1
     [76] rematch2_2.1.2         gtable_0.3.6           tibble_3.3.1
     [79] pillar_1.11.1          htmltools_0.5.9        R6_2.6.1
     [82] evaluate_1.0.5         lattice_0.22-9         readr_2.2.0
     [85] markdown_2.0           haven_2.5.5            backports_1.5.1
     [88] cards_0.8.1            broom_1.0.13           snakecase_0.11.1
     [91] rstantools_2.7.0       DescTools_0.99.60      class_7.3-23
     [94] Rcpp_1.1.2             coda_0.19-4.1          xfun_0.60
     [97] fs_2.1.0               zoo_1.9-0              forcats_1.0.1
    [100] pkgconfig_2.0.3       
