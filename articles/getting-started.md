# Getting started with SciDataReportR

SciDataReportR helps researchers turn labelled scientific and clinical
data into reproducible exploratory analyses, statistical comparisons,
scientific visualizations, and report-ready outputs.

This article is adapted from the SciDataReportR R/Medicine workflow
narrative. It follows the common path from a raw researcher-native data
frame to metadata-aware reporting and reusable analysis outputs.

The examples use the current workflow-oriented function names. Older
names such as
[`CreatePCATable()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md),
[`CreateZScorePlot()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md),
and
[`Make_DataDictionary()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeDataDictionary.md)
remain available as compatibility aliases.

## Why this workflow matters

Scientific and clinical datasets often contain many variables, labels,
recoding rules, caveats, missing values, and potential relationships.
Before formal statistical modeling, researchers need to understand the
structure of the data, identify quality issues, inspect distributions,
screen associations, and decide whether covariates or dimensionality
reduction are needed.

SciDataReportR organizes that work into repeatable steps:

1.  Import data and create metadata.
2.  Recode and relabel variables.
3.  Profile missingness and distributions.
4.  Create summary and comparison tables.
5.  Screen associations across variable families.
6.  Build dimensionality reduction or projection workflows.
7.  Save reproducible, report-ready outputs.

## Load data

Most SciDataReportR workflows start with a data frame and, when
available, a variable-type table or codebook.

``` r

library(SciDataReportR)
```

``` r

data("SampleData")
data("SampleVariableTypes")

df_Raw <- SampleData
VariableTypes <- SampleVariableTypes
```

> **Note**
>
> The examples use the raw `SampleData` dataset and its shipped
> `SampleVariableTypes` codebook. The raw `age` variable contains `999`
> as a sentinel for missing data;
> [`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
> uses the codebook to convert those values to `NA` before any
> summaries, plots, or models are created.

For a project dataset, the starting point might look like this:

``` r

df_Clinical <- readr::read_csv(here::here("data", "clinical_data.csv"))
VariableTypes_Clinical <- readr::read_csv(
  here::here("data", "variable_types.csv")
)
```

## Create variable metadata

Variable type templates make the analysis plan explicit. They help
distinguish continuous, categorical, binary, ordinal, outcome,
covariate, and feature variables before plotting or testing.

``` r

VariableTypes_Template <- CreateVariableTypesTemplate(df_Raw)
```

[`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md)
is useful when a dataset does not yet have a codebook. For the remaining
examples, use the shipped `SampleVariableTypes` object because it
already contains the labels, recodes, and missing-value rules for
`SampleData`.

Data dictionaries summarize the available variables and make the dataset
easier to review with collaborators.

``` r

df_DataDictionary <- MakeDataDictionary(df_Raw)
FormattedDataDictionary(df_Raw)
```

## Recode and relabel data

When a variable-type table contains recoding or label information, use
it to produce an analysis-ready data frame while preserving scientific
meaning.

``` r

df_Revalued <- RevalueData(df_Raw, VariableTypes)$RevaluedData
```

## Inspect missingness

Missing data patterns are one of the first bottlenecks in scientific
analysis. They affect which variables are interpretable, which models
are feasible, and whether downstream comparisons are biased.

``` r

PlotMissingData(df_Revalued, Relabel = TRUE)
```

## Create report-ready summaries

Table 1 style summaries help reviewers and collaborators understand the
cohort before detailed analyses.

``` r

MakeTable1(df_Revalued, TreatOrdinalAs = "Continuous")
```

Continuous and categorical profiling functions help researchers inspect
variable distributions without hand-building repeated plots.

``` r

vars_Continuous <- getNumVars(df_Revalued, Ordinal = FALSE)
vars_Categorical <- getCatVars(df_Revalued)

CreateSummaryTable(df_Revalued, vars_Continuous, Relabel = TRUE)
PlotContinuousDistributions(df_Revalued, vars_Continuous[1:12], ncol = 3)
PlotCategoricalDistributions(df_Revalued, vars_Categorical)
```

## Screen associations

SciDataReportR supports both focused plots and matrix-style screening.
This helps researchers identify relationships, candidate covariates, and
patterns that may need dimensionality reduction.

``` r

PlotAssociations(df_Revalued, "age", "Adiponectin")
PlotAssociations(df_Revalued, "Diagnosis", "Ab_42")
```

Correlation heatmaps return structured objects that can be reused
downstream.

``` r

correlation_result <- PlotCorrelationsHeatmap(
  df_Revalued,
  xVars = vars_Continuous[1:5],
  yVars = vars_Continuous[20:40],
  method = "pearson",
  covars = NULL,
  Relabel = TRUE,
  Ordinal = FALSE
)

correlation_result$Unadjusted$plot
```

The output of
[`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md)
can be passed to
[`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md),
and the returned plot can use
[`geom_starcaption()`](https://rdastgh1.github.io/SciDataReportR/reference/geom_starcaption.md)
to explain significance stars.

``` r

add_r_and_stars(correlation_result) + geom_starcaption()
```

## Compare groups

Group comparisons are a central part of clinical and life science
reporting.
[`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
provides report-ready summaries with optional covariates, effect sizes,
and pairwise contrasts.

``` r

MakeComparisonTable(
  DataFrame = df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c("age", "tau", "p_tau"),
  AddEffectSize = TRUE
)
```

For high-dimensional group signals, z-score plots provide a compact
visual summary across many variables.

``` r

vars_LabMeasures <- vars_Continuous[10:60]

PlotZScore(
  df_Revalued,
  TargetVar = "Diagnosis",
  Variables = vars_LabMeasures,
  sort = FALSE
)
```

## Build reusable dimensionality reduction workflows

PCA workflows can produce scree plots, loading summaries, and reusable
projection objects. This is useful when a cohort-level transformation
should be applied to future datasets.

``` r

pca_object <- CreatePCAObject(
  df_Revalued,
  vars_LabMeasures,
  minThresh = 0.75,
  scale = TRUE,
  center = TRUE
)

pca_object$p_scree
pca_object$Lollipop
```

Projection functions make the fit/apply split explicit:

``` r

projected_scores <- ProjectPCA(new_data, pca_object)
```

## Reproducibility

SciDataReportR workflows are designed to make repeated analysis steps
visible and reusable. For reports, include session information and
preserve the variable-type table, codebook, and any fit objects used for
projections.

``` r

# save.image(here::here("results", "getting-started.RData"))
print(sessionInfo())
```

## Next steps

Use the reference index by workflow family:

- Data setup, metadata, and codebooks.
- Preprocessing and data quality.
- Statistical comparison workflows.
- Association, regression, and interaction workflows.
- Visualization functions.
- Dimensionality reduction, projection, and clustering.
- Longitudinal and temporal workflows.

Future article additions should expand this starter flow into focused
workflows for codebook harmonization, visualization galleries,
PCA/projection, SOM projection, normative T-scores, and longitudinal
transitions.
