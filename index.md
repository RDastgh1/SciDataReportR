# SciDataReportR

SciDataReportR is an R workflow infrastructure package for turning
labelled life science data frames into reproducible exploratory
analyses, statistical comparisons, dimensionality reduction workflows,
clustering/projection pipelines, scientific visualizations, and
report-ready outputs.

The package is designed for researchers working with practical,
researcher-native data structures: data frames, labelled variables,
codebooks, variable-type templates, cohorts, visits, outcomes,
covariates, and high-dimensional feature sets.

## Why SciDataReportR?

SciDataReportR helps life science researchers operationalize repeated
analysis patterns without rebuilding the same reporting and exploration
code for every study.

- Reduces repeated manual EDA, metadata cleanup, comparison tables,
  visual summaries, and report assembly.
- Keeps labels, variable types, codebooks, and data dictionaries close
  to the analysis.
- Supports high-dimensional scientific screening across outcomes,
  cohorts, covariates, and feature sets.
- Enables reusable transformations and projections for new cohorts.
- Follows practical R package conventions from the tidyverse, usethis,
  roxygen2, and pkgdown ecosystem.

## Installation

You can install the development version of SciDataReportR from
[GitHub](https://github.com/RDastgh1/SciDataReportR) with:

``` r

# install.packages("devtools")
devtools::install_github("RDastgh1/SciDataReportR")
```

## Scientific Workflows

SciDataReportR functions are intended to be used as workflow components
rather than alphabetical utilities.

| Workflow family | Representative functions |
|----|----|
| Data setup and metadata | [`CreateProjectFolders()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateProjectFolders.md), [`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md), [`MakeDataDictionary()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeDataDictionary.md), [`FormattedDataDictionary()`](https://rdastgh1.github.io/SciDataReportR/reference/FormattedDataDictionary.md), [`UpdateDataDictionary()`](https://rdastgh1.github.io/SciDataReportR/reference/UpdateDataDictionary.md) |
| Codebook and harmonization | [`AddToCodebook()`](https://rdastgh1.github.io/SciDataReportR/reference/AddToCodebook.md), [`UpdateCodebook()`](https://rdastgh1.github.io/SciDataReportR/reference/UpdateCodebook.md), [`CombineCodebooks()`](https://rdastgh1.github.io/SciDataReportR/reference/CombineCodebooks.md), [`MergeCodebooks()`](https://rdastgh1.github.io/SciDataReportR/reference/MergeCodebooks.md), [`CodebookMergeApp()`](https://rdastgh1.github.io/SciDataReportR/reference/CodebookMergeApp.md) |
| Cleaning and preprocessing | [`ReplaceMissingCode()`](https://rdastgh1.github.io/SciDataReportR/reference/ReplaceMissingCode.md), [`ReplaceMissingLabels()`](https://rdastgh1.github.io/SciDataReportR/reference/ReplaceMissingLabels.md), [`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md), [`ReValueFactors()`](https://rdastgh1.github.io/SciDataReportR/reference/ReValueFactors.md), [`ConvertOrdinalToNumeric()`](https://rdastgh1.github.io/SciDataReportR/reference/ConvertOrdinalToNumeric.md), [`windsorize()`](https://rdastgh1.github.io/SciDataReportR/reference/windsorize.md), [`IQROutliers()`](https://rdastgh1.github.io/SciDataReportR/reference/IQROutliers.md) |
| Exploratory profiling | [`PlotMissingData()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotMissingData.md), [`PlotContinuousDistributions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotContinuousDistributions.md), [`PlotCategoricalDistributions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCategoricalDistributions.md), [`PlotTimeDistribution()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotTimeDistribution.md), [`PlotMiningMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotMiningMatrix.md) |
| Statistical comparisons | [`MakeTable1()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeTable1.md), [`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md), [`MakeFacetCatComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeFacetCatComparisonTable.md), [`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md), [`Plot2GroupStats()`](https://rdastgh1.github.io/SciDataReportR/reference/Plot2GroupStats.md) |
| Association and regression mining | [`PlotAssociations()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotAssociations.md), [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md), [`PlotDirectionalHeatmaps()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotDirectionalHeatmaps.md), [`UnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/UnivariateRegressionTable.md), [`plotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/plotForestFromTable.md) |
| Dimensionality reduction and projection | [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md), [`plotPCA()`](https://rdastgh1.github.io/SciDataReportR/reference/plotPCA.md), [`ExtractPCAComponentSummary()`](https://rdastgh1.github.io/SciDataReportR/reference/ExtractPCAComponentSummary.md), [`ProjectPCA()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectPCA.md), [`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md), [`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md) |
| Normative modeling | [`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md), [`ApplyNormativeTScores()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyNormativeTScores.md) |
| Clustering and cohort projection | [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md), [`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md) |
| Longitudinal and temporal workflows | [`Merge_ByClosestTime()`](https://rdastgh1.github.io/SciDataReportR/reference/Merge_ByClosestTime.md), [`SummarizeTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/SummarizeTransitions.md), [`PlotSwimmerTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSwimmerTransitions.md) |
| Reporting | [`CreateSummaryTable()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSummaryTable.md), [`CreateStatisticsTable()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateStatisticsTable.md), [`use_EDATemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/use_EDATemplate.md) |

## Visualization Gallery

Visualization is a cross-cutting SciDataReportR workflow layer.
Plot-producing functions also appear in their scientific workflow
families because scientists often discover methods through the figures
they need.

| Visualization need | Functions |
|----|----|
| Data quality | [`PlotMissingData()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotMissingData.md), [`PlotContinuousDistributions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotContinuousDistributions.md), [`PlotCategoricalDistributions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCategoricalDistributions.md), [`IQROutliers()`](https://rdastgh1.github.io/SciDataReportR/reference/IQROutliers.md) |
| Group comparison | [`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md), [`Plot2GroupStats()`](https://rdastgh1.github.io/SciDataReportR/reference/Plot2GroupStats.md), [`PlotPValueComparisons()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPValueComparisons.md), [`PlotSplitViolin()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSplitViolin.md) |
| Association and correlation | [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md), [`PlotPhiHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPhiHeatmap.md), [`PlotPointCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPointCorrelationsHeatmap.md), [`PlotDirectionalHeatmaps()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotDirectionalHeatmaps.md), [`PlotAssociations()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotAssociations.md), [`plotSigCorrelations()`](https://rdastgh1.github.io/SciDataReportR/reference/plotSigCorrelations.md), [`plotSigAssociations()`](https://rdastgh1.github.io/SciDataReportR/reference/plotSigAssociations.md) |
| Regression and interaction | [`PlotPartialRegressionScatter()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPartialRegressionScatter.md), [`plotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/plotForestFromTable.md), [`PlotInteractionEffectsContinuous()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotInteractionEffectsContinuous.md), [`PlotInteractionEffectsMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotInteractionEffectsMatrix.md), [`PlotCatInteractionEffectsMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCatInteractionEffectsMatrix.md), [`PlotNumInteractionEffectsMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotNumInteractionEffectsMatrix.md) |
| Dimensionality reduction and clustering | [`plotPCA()`](https://rdastgh1.github.io/SciDataReportR/reference/plotPCA.md), plots produced from [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md) workflows |
| Longitudinal and temporal | [`PlotTimeDistribution()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotTimeDistribution.md), [`PlotSwimmerTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSwimmerTransitions.md) |
| Specialized scientific plots | [`PlotBlandAltman()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotBlandAltman.md), [`PlotSpiderChart()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSpiderChart.md), [`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md) |

Planned visualization extensions include heatmaps, volcano plots for
two-group categorical contrasts, and volcano-style plots for continuous
predictors with beta values on the x-axis.

## Implemented Methods

SciDataReportR operationalizes practical scientific methods that
commonly sit between raw data cleaning and manuscript-ready reporting:

- Codebook, data dictionary, and labelled-variable workflows.
- Cohort comparison with covariates, effect sizes, pairwise contrasts,
  and report-ready tables.
- Multi-variable z-score comparisons with p-value and FDR signals.
- Correlation, association, regression, and interaction screening.
- PCA creation, summary, visualization, and projection.
- SOM clustering and projection to new data.
- Normative T-score creation and application.
- Longitudinal transition summaries and swimmer-style plotting.
- Scientific visualization outputs tied to reproducible workflows.

## Function Inputs and Workflow Dependencies

Some functions are designed to consume outputs from other SciDataReportR
functions. These relationships are part of the workflow infrastructure
and should be checked when composing analyses.

SciDataReportR now uses workflow-oriented canonical names for new code,
while older public names remain available as compatibility aliases. For
example,
[`CreatePCATable()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCATable.md)
still works, but
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md)
better describes the reusable PCA object returned by that workflow.

| Function | Primary input | Upstream SciDataReportR function | Output | Common downstream use |
|----|----|----|----|----|
| [`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md) | Correlation heatmap result list | [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md) | Annotated `ggplot` | Add [`geom_starcaption()`](https://rdastgh1.github.io/SciDataReportR/reference/geom_starcaption.md) for report-ready star explanations |
| [`geom_starcaption()`](https://rdastgh1.github.io/SciDataReportR/reference/geom_starcaption.md) | No direct object; added with `+` | [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md) or [`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md) | `labs()` caption layer | Explain `*`, `**`, and `***` thresholds on heatmaps |
| [`ApplyNormativeTScores()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyNormativeTScores.md) | New data and normative scoring object | [`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md) | Applied normative T-scores | Score future cohorts with the same model |
| [`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md) | New data and z-score parameters | [`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md) | Projected standardized scores | Apply a prior standardization to new data |
| [`ProjectPCA()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectPCA.md) | New data and PCA object | [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md) | Projected PCA scores | Compare new samples in an existing PCA space |
| [`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md) | New data and SOM cluster solution | [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md) | Projected SOM/cluster assignments | Apply an existing cluster solution to future cohorts |
| [`ProjectRCI()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectRCI.md) | New data and RCI object | [`CreateRCIObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateRCIObject.md) | Projected reliable change metrics | Apply an existing change model |

## Flagship Workflow Examples

``` r

library(SciDataReportR)

# Metadata-aware setup
variable_types <- CreateVariableTypesTemplate(SampleData)
data_dictionary <- MakeDataDictionary(SampleData)

# Cohort comparison workflow
comparison <- MakeComparisonTable(
  DataFrame = SampleData,
  CompVariable = "Diagnosis",
  Variables = c("age", "tau", "p_tau")
)

# Correlation heatmap workflow with downstream annotation
cor_res <- PlotCorrelationsHeatmap(
  Data = SampleData,
  xVars = c("age", "tau", "p_tau"),
  yVars = c("Ab_42", "C_Reactive_Protein")
)

annotated_plot <- add_r_and_stars(cor_res) + geom_starcaption()
```

## Data Structures and Conventions

SciDataReportR is built around practical life science data structures:

- Data frames or tibbles with one row per observation, sample,
  participant, or visit.
- Labelled variables where labels should remain visible in reports and
  plots.
- Variable-type templates that distinguish continuous, categorical,
  binary, ordinal, outcome, covariate, and feature variables.
- Codebooks and data dictionaries that document the meaning and
  structure of datasets.
- Fit/project workflows where a transformation or model learned from one
  cohort can be applied to a future cohort.

## Documentation

- Package website: <https://rdastgh1.github.io/SciDataReportR/>
- Reference index:
  <https://rdastgh1.github.io/SciDataReportR/reference/>
- Issues and feature requests:
  <https://github.com/RDastgh1/SciDataReportR/issues>
- R/Medicine 2024 slides:
  <https://rdastgh1.quarto.pub/rmedicine-2024-scidatareportr>

## Roadmap

Near-term priorities:

- Expand workflow-oriented vignettes for data setup, codebook
  harmonization, EDA/reporting, comparisons, visualizations,
  PCA/projection, SOM projection, normative T-scores, and longitudinal
  transitions.
- Improve function-level documentation with input object requirements,
  return structures, and downstream use.
- Add a richer visualization gallery with representative outputs.
- Strengthen release notes, CI, pkgdown organization, and
  contributor-facing documentation.

Future directions:

- Heatmap workflows for high-dimensional scientific data.
- Volcano plot workflows for categorical two-group differences.
- Volcano-style workflows for continuous predictors with beta values on
  the x-axis.
- Workflow templates for cohort comparison, biomarker discovery,
  longitudinal follow-up, normative scoring, and cluster projection.
- Structured analysis/report cards that summarize inputs, methods,
  parameters, and outputs.
- Optional AI-native report summaries and machine-readable workflow
  metadata.

## Citation, License, and Contributing

SciDataReportR is released under the MIT license. Please cite the
package and link to the project website when using it in scientific
reports, presentations, or manuscripts.

Contributions, bug reports, and workflow ideas are welcome through
[GitHub issues](https://github.com/RDastgh1/SciDataReportR/issues).
