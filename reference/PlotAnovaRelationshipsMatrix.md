# Plot ANOVA Relationships Matrix

This function plots the relationship between continuous and categorical
variables using ANOVA or Kruskal-Wallis tests. It generates a "heatmap"
with points colored and shaped based on statistical significance and
effect size. When generalized eta-squared is unavailable, raw p-values
are colored with
[`scale_color_pvalue()`](https://rdastgh1.github.io/SciDataReportR/reference/scale_color_pvalue.md);
the FDR-corrected plot uses adjusted p-values.

## Usage

``` r
PlotAnovaRelationshipsMatrix(
  data,
  CatVars,
  ContVars,
  covariates = NULL,
  Relabel = TRUE,
  Parametric = TRUE,
  Ordinal = FALSE,
  min_n = 4,
  eps = 1e-08,
  Data = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  Covariates = lifecycle::deprecated()
)
```

## Arguments

- data:

  The data frame containing the variables of interest.

- CatVars:

  Character vector of categorical variable names.

- ContVars:

  Character vector of continuous variable names.

- covariates:

  Optional character vector of covariate names for ANCOVA analysis.

- Relabel:

  Logical indicating whether to relabel variables with their labels
  (default is TRUE).

- Parametric:

  Logical indicating whether to use parametric (ANOVA) or non-parametric
  (Kruskal-Wallis) tests (default is TRUE).

- Ordinal:

  Logical, indicating whether ordinal variables should be considered.

- min_n:

  Minimum number of complete observations required for a tested
  relationship.

- eps:

  Small positive value used to avoid zero-size plotting artifacts.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- fdr_scope:

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all p-values at once (historical behavior).
  `"per_outcome"` corrects separately within each outcome: outcomes are
  the continuous variables (`ContVars`).

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

A list containing three ggplot objects: p (scatter plot without multiple
comparison correction), p_FDR (scatter plot with FDR correction), and
pvaltable (data frame of p-values and significance).

## Parametric and non-parametric tests

ANOVA assumes roughly normal residuals and similar variance across
groups. Biomarker concentrations are often right-skewed, and group sizes
are often uneven, both of which make the F-test unreliable.

`Parametric = FALSE` switches to Kruskal-Wallis, which compares ranks
rather than means and assumes neither normality nor equal variance. The
matrix is built and read exactly the same way.

Comparing the two p-value tables shows which conclusions depended on the
choice of test. Relationships that hold under both are the ones to
trust. Where they disagree, the direction is informative: significant
only under ANOVA suggests the result is being carried by skew or a few
extreme values, while significant only under Kruskal-Wallis suggests a
real shift in the bulk of the distribution that outliers were masking
from the mean-based test.

Kruskal-Wallis is a rank test and has no ANCOVA form, so `covariates`
requires `Parametric = TRUE`.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotAnovaRelationshipsMatrix(
  Labelled,
  CatVars = c("Diagnosis", "sex", "Genotype"),
  ContVars = c("age", "ACE_CD143_Angiotensin_Converti",
               "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
               "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
               "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin",
               "Apolipoprotein_A1")
)

# Raw p-value associations
result$Unadjusted$plot


# FDR-adjusted associations
result$FDRCorrected$plot
#> Warning: Using size for a discrete variable is not advised.


# The same matrix using Kruskal-Wallis instead of ANOVA
result_NonParametric <- PlotAnovaRelationshipsMatrix(
  Labelled,
  CatVars = c("Diagnosis", "sex", "Genotype"),
  ContVars = c("age", "ACE_CD143_Angiotensin_Converti",
               "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
               "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
               "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin",
               "Apolipoprotein_A1"),
  Parametric = FALSE
)

result_NonParametric$Unadjusted$plot

result_NonParametric$FDRCorrected$plot
#> Warning: Using size for a discrete variable is not advised.


# Where the two tests disagree
cols_Key <- c("CategoricalVariable", "ContinuousVariable", "p")
df_Compare <- merge(
  result$Unadjusted$PvalTable[, cols_Key],
  result_NonParametric$Unadjusted$PvalTable[, cols_Key],
  by = c("CategoricalVariable", "ContinuousVariable"),
  suffixes = c("_ANOVA", "_KruskalWallis")
)
df_Compare$AgreesAt05 <-
  (df_Compare$p_ANOVA < 0.05) == (df_Compare$p_KruskalWallis < 0.05)

htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(
    dplyr::mutate(
      df_Compare,
      dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 3))
    ),
    height = "320px", full_width = TRUE
  )
)))


 CategoricalVariable 
```
