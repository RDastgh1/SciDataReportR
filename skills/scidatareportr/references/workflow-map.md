# SciDataReportR workflow map

Use this map to choose a workflow from the analytical goal. Read `api-reference.md` for exact arguments, return structures, and specialized options.

## Intake, metadata, and codebooks

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Inspect/import a CSV, Excel, SPSS, or common scientific export | `InspectFile()` then `ReadSciData()` | Inspect problematic files before importing; then build metadata with `CreateVariableTypesTemplate()` if no codebook exists. |
| Apply a codebook to values and labels | `RevalueData()` | Use immediately after data and codebook loading, before EDA. Continue with `$RevaluedData`. |
| Create, display, or update a data dictionary | `CreateVariableTypesTemplate()`, `MakeDataDictionary()`, `FormattedDataDictionary()`, `AddToCodebook()`, `UpdateCodebook()` | Keep the codebook authoritative as derived variables are added. |
| Harmonize study codebooks | `CombineCodebooks()`, `MergeCodebooks()`, `CodebookMergeApp()` | Resolve conflicts before combining/revaluing datasets. |
| Repair labelled values or classify variables | `ReplaceMissingCode()`, `ReplaceMissingLabels()`, `ReValueFactors()`, `ConvertOrdinalToNumeric()`, `getNumVars()`, `getCatVars()`, `getBinaryVars()` | Preserve labels; choose an ordinal representation deliberately. |

## Merge and cohort quality control

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Audit a planned or completed merge | `ValidateMerge()`, `safe_merge()`, `PlotMergeValidation()`, `ExploreMergeValidation()` | Check keys, duplicates, unmatched IDs, rows, and variables before using the merged frame. |
| Compare dataset versions or cohorts | `CompareDatasets()`, `ExploreDatasetComparison()`, `PlotDatasetComparison()` | Use a defined key and review discrepancies before pooling. |
| Merge fragmented or time-adjacent records | `MergeFragmentedRecords()`, `Merge_ByClosestTime()` | Document the matching rule and validate results afterward. |

## EDA and descriptive reporting

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Assess missingness | `PlotMissingData()` | Run after codebook-based recoding. |
| Profile distributions | `PlotContinuousDistributions()`, `PlotCategoricalDistributions()`, `PlotTimeDistribution()`, `IQROutliers()` | Select variables using the type helpers where useful. |
| Produce summaries or Table 1 | `CreateSummaryTable()`, `MakeTable1()`, `CreateStatisticsTable()` | Leave labels attached so human-facing output uses them. |
| Assemble a figure panel | `AssemblePlots()` | Supply completed ggplot objects; do not save plots unless requested. |

## Comparisons, associations, and regression

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Compare variables across a group | `MakeComparisonTable()` | Use `group_var`, `variables`, and optional `covariates`; request effect sizes/pairwise contrasts only when needed. |
| Visualize two-group or distribution differences | `Plot2GroupStats()`, `PlotSplitViolin()`, `PlotZScore()`, `PlotPValueComparisons()` | Pair with a defined comparison question, not indiscriminate significance hunting. |
| Screen correlations/associations | `PlotCorrelationsHeatmap()`, `PlotAssociations()`, `PlotMiningMatrix()` | Define the FDR family intentionally; use `add_r_and_stars()` or significant-relationship plotters downstream only after reviewing results. |
| Run univariable screening | `MakeUnivariateRegressionTable()` | Use `ApplyFDRCorrection()` on the appropriate result family, then optionally `PlotForestFromTable()`. |
| Fit mutually adjusted or penalized models | `MultivariableRegressionTable()` | Use only when the modeling question calls for joint adjustment; inspect model diagnostics/results. |
| Explore interaction effects | `PlotInteractionEffectsContinuous()`, `PlotInteractionEffectsMatrix()`, `PlotCatInteractionEffectsMatrix()`, `PlotNumInteractionEffectsMatrix()` | Specify outcome, predictor, and covariates with current argument names. |

## Biomarkers, projection, and clustering

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Evaluate one or many biomarkers | `EvaluateBiomarkerPerformance()`, `ScreenBiomarkerPerformance()` | Choose binary versus continuous outcome handling and validation deliberately. |
| Create and project standardized scores | `CreateZScoreObject()` then `ProjectZScore()` | Reuse the fitted object for future cohorts rather than recalculating parameters. |
| Create and project PCA | `CreatePCAObject()` then `ProjectPCA()`; use `plotPCA()` and `ExtractPCAComponentSummary()` to interpret | Fit in the training cohort and project new data with the same object. `CreatePCATable()` is a compatibility alias. |
| Build/project normative T scores or reliable change | `CreateNormativeTScoreModel()` then `ApplyNormativeTScores()`; `CreateRCIObject()` then `ProjectRCI()` | Keep the trained object and document the reference cohort. |
| Cluster a cohort | `CreateClusterModel_*()` selected by data type; `ProjectCluster()` for future cohorts | Choose continuous, mixed, categorical, density, PCA, SOM, or latent-class workflows based on the data and scientific question. Review stability and fit before finalizing. |
| Use the legacy SOM/LPA workflow | `CreateClusterModel_SOM_MClust()` | Prefer the current constructor over deprecated `Pipeline_SOMClust()`. |

## Longitudinal and specialized figures

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Summarize/reveal transitions | `SummarizeTransitions()`, `PlotSwimmerTransitions()`, `PlotTimeSwimmer()` | Confirm participant IDs, ordering, and visit/time variables first. |
| Specialized plot selection | Search `api-reference.md` for Bland–Altman, pathway, spider, volcano, directional/phi/point heatmaps, swimmer, and cluster figures | Use the dedicated vignette/Rd documentation before selecting a specialized plot. |

## Object flow rules

- `Create*Object()` / `CreateClusterModel_*()` fits an object; use its matching `Project*()` or `Apply*()` function for new data.
- `RevalueData()` produces the labelled analysis frame consumed by EDA, tables, and plots.
- Result objects from correlation and regression workflows feed their documented annotation, FDR, and forest-plot helpers.
- The full generated catalog is the fallback for any export not named above.
