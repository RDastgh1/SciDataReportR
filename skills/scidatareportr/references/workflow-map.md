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
| Compare several groups with one reference | `MakePairwiseHeatmap()` | Supply the scientific referent; its scaling defines all displayed contrasts. Choose the FDR family across groups, variables, or the displayed matrix intentionally. |
| Screen continuous-variable associations | `PlotCorrelationsHeatmap()`, `PlotAssociations()`, `PlotMiningMatrix()` | Define the FDR family intentionally. For symmetric matrices, `triangle = "upper"` or `"lower"` displays each pair once without changing returned results. |
| Screen binary or mixed-variable associations | `PlotPhiHeatmap()`, `PlotPointCorrelationsHeatmap()`, `PlotDirectionalHeatmaps()` | Use Phi for binary-binary and point-biserial for binary-continuous pairs; use the directional heatmap only when one mixed-variable overview is needed. Preserve the binary mapping so the positive level remains interpretable. |
| Compare correlations across two independent groups | `PlotCorrelationComparisons()` | Define reference and comparison groups explicitly when needed; DeltaR is comparison minus reference. Inspect per-cell test availability/status and approximation metadata before interpreting p-values. |
| Run univariable screening | `MakeUnivariateRegressionTable()` | Use `ApplyFDRCorrection()` on the appropriate result family, then optionally `PlotForestFromTable()`. |
| Screen many effects against one outcome | `PlotVolcanoEffects()` | Choose outcome type and effect metric deliberately; inspect FDR results rather than treating labels or color as the analysis. |
| Fit mutually adjusted or penalized models | `MultivariableRegressionTable()` | Use only when the modeling question calls for joint adjustment; inspect model diagnostics/results. |
| Explore interaction effects | `PlotInteractionEffectsContinuous()`, `PlotInteractionEffectsMatrix()`, `PlotCatInteractionEffectsMatrix()`, `PlotNumInteractionEffectsMatrix()` | Specify outcome, predictor, and covariates with current argument names. |
| Evaluate diagnostic evidence from categorical test results | `DiagnosticLikelihoodRatioTable()` then `PlotDiagnosticLRHeatmap()` or `PlotDiagnosticLRForest()` | This is diagnostic accuracy, not a nested-model test. Define outcome-positive and binary-predictor positive levels before calculation; retain zero/infinite likelihood ratios unless a justified continuity correction is requested. |
| Assess agreement of two measurement methods | `PlotBlandAltman()` | Use agreement limits and mean-difference structure; do not substitute a correlation for method agreement. |

## Biomarkers, projection, and clustering

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Evaluate one or many biomarkers | `EvaluateBiomarkerPerformance()`, `ScreenBiomarkerPerformance()` | Choose binary versus continuous outcome handling and validation deliberately. |
| Create and project standardized scores | `CreateZScoreObject()` then `ProjectZScore()` | Reuse the fitted object for future cohorts rather than recalculating parameters. |
| Create and project PCA | `CreatePCAObject()` then `ProjectPCA()`; use `plotPCA()` and `ExtractPCAComponentSummary()` to interpret | Fit in the training cohort and project new data with the same object. `CreatePCATable()` is a compatibility alias. |
| Reduce nominal categorical measures | `CreateMCAObject()` | Use MCA for interpretable categorical dimensions without clustering; review scree, loadings, and category contributions before interpreting dimensions. |
| Build/project normative T scores or reliable change | `CreateNormativeTScoreModel()` then `ApplyNormativeTScores()`; `CreateRCIObject()` then `ProjectRCI()` | Keep the trained object and document the reference cohort. |
| Cluster a cohort | `CreateClusterModel_*()` selected by data type; `ProjectCluster()` for future cohorts | Choose continuous, mixed, categorical, density, PCA, SOM, or latent-class workflows based on the data and scientific question. For many nominal variables requiring both reduction and clustering, use `CreateClusterModel_MCA_MClust()` rather than standalone MCA. Review stability and fit before finalizing. |
| Use the legacy SOM/LPA workflow | `CreateClusterModel_SOM_MClust()` | Prefer the current constructor over deprecated `Pipeline_SOMClust()`. |

## Imaging, pathway, longitudinal, and specialized figures

| Goal | Prefer | Prerequisite and downstream flow |
| --- | --- | --- |
| Derive FreeSurfer bilateral measures or ICV ratios | `DeriveFreesurferVolumes()` | Verify source ASEG/DKT names and the ICV variable; keep the returned `Freesurfer_derivation_log` attribute with derived measures. |
| Compare metabolites and display kynurenine-pathway effects | `calculate_pathway_results()` then `PlotPathway_KT()` | Define whether the comparison is binary or continuous, use the resulting fold-change/correlation table directly, and choose raw versus FDR significance explicitly. |
| Summarize/reveal transitions | `SummarizeTransitions()`, `PlotSwimmerTransitions()`, `PlotTimeSwimmer()` | Confirm participant IDs, ordering, and visit/time variables first. |
| Specialized plot selection | Search `api-reference.md` for spider, swimmer, and cluster figures | Use the dedicated vignette/Rd documentation before selecting a specialized plot. |

## Object flow rules

- `Create*Object()` / `CreateClusterModel_*()` fits an object; use its matching `Project*()` or `Apply*()` function for new data.
- `RevalueData()` produces the labelled analysis frame consumed by EDA, tables, and plots.
- Result objects from correlation, regression, diagnostic-LR, pathway, and clustering workflows feed their documented annotation, visualization, and projection helpers.
- The full generated catalog is the fallback for any export not named above.
