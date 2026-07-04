# SciDataReportR — Exported API Audit

**Date:** 2026-07-03
**Scope:** All 108 exports in `NAMESPACE` (107 functions + the `%!in%` operator), audited against source in `R/` and generated docs in `man/`. Audit only — no code was modified.

**Method:** Every exported function's formals (names, defaults, order) were extracted by parsing the R source; roxygen blocks and generated `.Rd` files were cross-checked for `@param`/`@return`/`@examples` coverage and default agreement; `codetools::findGlobals()` was run against the set of defined + imported symbols to find dangling references.

---

## 0. Function name inventory (casing of the exports themselves)

| Style | Count | Share | Examples |
|---|---|---|---|
| `PascalCase` | 81 | 75% | `PlotMissingData`, `CreateSummaryTable`, `ValidateMerge` |
| `camelCase` | 12 | 11% | `getNumVars`, `plotPCA`, `createBinaryMapping`, `windsorize`, `removeString`, `reverseFactorLevels` |
| `snake_case` | 7 | 6.5% | `safe_merge`, `merge_detail`, `merge_summary_table`, `add_r_and_stars`, `calculate_pathway_results`, `use_EDATemplate`, `geom_starcaption` |
| `Pascal_snake` hybrid | 7 | 6.5% | `Make_DataDictionary`, `Merge_ByClosestTime`, `Pipeline_SOMClust`, `Project_SOMClust`, `Project_ZScore`, `CreatePathwayPlot_KT`, `PlotPathway_KT` |
| operator | 1 | 1% | `%!in%` |

Of the 7 hybrids, 5 are documented compatibility aliases (see §4.1), so the "real" API is essentially PascalCase with a camelCase minority concentrated in the `get*`/`plot*` helpers and a snake_case minority concentrated in the newest additions (the `safe_merge` family) plus the pathway/metabolomics functions.

---

## 1. Parameter naming consistency

Across the 107 exported functions there are **700 parameter slots** using **429 distinct names**.

### 1.1 Casing conventions

Counting every parameter slot:

| Style | Slots | Share |
|---|---|---|
| `PascalCase` (`Data`, `Relabel`, `TopN`) | 430 | 61% |
| `snake_case` (`id_var`, `point_size`, `names_prefix`) | 142 | 20% |
| single lowercase word (`df`, `covars`, `scale`, `center`, `method`) | 89 | 13% |
| `camelCase` (`xVars`, `interVar`, `numComponents`, `minThresh`) | 37 | 5% |
| `Mixed_Case_underscore` (`SigP_YCoord`, `SigFDR_YCoord` in `PlotZScore`) | 2 | <1% |

Per function: only **38 of 97** functions with named parameters use a single style throughout; **52 mix styles within one signature**. Worst offenders (4 styles in one signature): `Plot2GroupStats`, `PlotChiSqCovar`, `PlotCorrelationsHeatmap`, `PlotInteractionEffectsContinuous`.

Notable hazards:

- **`CreatePCAObject` has both `ImputeMethod` and `imputeMethod`** in the same signature (the lowercase one appears to be a legacy/compat argument). Two parameters differing only in case in one function is an active foot-gun for users.
- **`plotPCA` carries three generations of arguments**: deprecated `t = "NULLʺ` (a *string* `"NULL"` as default), legacy `HoverVar`, and current `HoverVars`/`ColorType`.
- Whole-cluster style splits: the Z-score/M-score/normative/SOM family (`CreateZScoreObject`, `CreateMScoreObject`, `CreateNormativeTScoreModel`, `CreateSOMClusterModel`, `MergeFragmentedRecords`, `PlotSwimmerTransitions`, `SummarizeTransitions`, `safe_merge`, `merge_*`) is written in snake_case/lowercase, while the tables/plots family is PascalCase.

### 1.2 Same concept, different names (semantic clusters)

**Cluster A — the primary data frame (first argument in most functions):**

| Name | Count | Functions |
|---|---|---|
| `Data` | ~31 | `CreateMCAObject`, `CreatePCAObject`, `CreateRCIObject`, `CreateSummaryTable`, `MultivariableRegressionTable`, `UnivariateRegressionTable`, most `Plot*` matrix/heatmap functions, `PlotSpiderChart`, `PlotTimeSwimmer`, `ProjectPCA`, `ProjectRCI`, `createBinaryMapping`, … (also `windsorize`, where `Data` is actually a **numeric vector**, not a data frame) |
| `DataFrame` | 18 | `MakeTable1`, `MakeComparisonTable`, `MakeFacetCatComparisonTable`, `MakeDataDictionary`, `FormattedDataDictionary`, `CreateVariableTypesTemplate`, `PlotAssociations`, `PlotBlandAltman`, `PlotCategoricalDistributions`, `PlotContinuousDistributions`, `PlotMissingData`, `PlotPartialRegressionScatter`, `createFacetLabels`, `getBinaryVars`, `getCatVars`, `getNumVars`, `plotSigAssociations`, `plotSigCorrelations` |
| `df` | 12 | `ApplyNormativeTScores`, `CreateMScoreObject`, `CreateNormativeTScoreModel`, `CreateSOMClusterModel`, `CreateZScoreObject`, `IQROutliers`, `MergeFragmentedRecords`, `PrepNumericData`, `ProjectZScore`, `ReplaceMissingLabels`, `calculate_pathway_results`, `reverseFactorLevels` |
| `data` | 3 | `PlotSplitViolin`, `PlotSwimmerTransitions`, `SummarizeTransitions` |
| `Dataframe` | 1 | `UpdateCodebook` (a fourth casing of the same word) |
| `DatatoRevalue` | 2 | `ReValueFactors`, `RevalueData` |
| role-named pairs | — | `OldData`/`NewData` (`CompareDatasets`), `LeftData`/`RightData`/`MergedData` (`ValidateMerge`), `DataFrame1`/`DataFrame2` (`Merge_ByClosestTime`), `df_before`/`df_add` (`safe_merge`), `new_df` (`ProjectSOMCluster`) |

**Cluster B — "which variables to operate on":**

- `Variables` (14 fns: `CreateSummaryTable`, `MakeTable1`, `PlotMissingData`, `PlotSpiderChart`, `PlotZScore`, …) vs lowercase `variables` (6 fns: z/M-score + SOM family, `PrepNumericData`, `reverseFactorLevels`)
- `VarsToReduce` (`CreateMCAObject`, `CreatePCAObject`, `ProjectPCA`)
- `xVars`/`yVars` (6 heatmap/matrix fns) — but `PlotVolcanoEffects` uses `xVars` + singular `yVar`
- `CatVars`/`ContVars` (`PlotAnovaRelationshipsMatrix`, `PlotPointCorrelationsHeatmap`, `PlotPhiHeatmap`, `createBinaryMapping`)
- `OutcomeVars`/`PredictorVars` (`UnivariateRegressionTable`, `MultivariableRegressionTable`, `PlotMiningMatrix`) vs camelCase `outcomeVars`/`predictorVars` (`PlotInteractionEffectsMatrix`) vs singular `outcomeVar`/`predictorVar` (`PlotInteractionEffectsContinuous`)
- Singles: `Variable` (`IQROutliers`), `Var` (`PlotSplitViolin`, `plotPCA`), `Var1`/`Var2` (`PlotAssociations`), `Variable1`/`Variable2` (`PlotBlandAltman`), `IndepVar`/`DepVar` (`PlotPartialRegressionScatter`), `test_var` (`CreateNormativeTScoreModel`), `metabolites` (`calculate_pathway_results`)

**Cluster C — covariates (four spellings of the same idea):**

- `Covariates` — 7 fns (`MakeComparisonTable`, `PlotMiningMatrix`, `PlotVolcanoEffects`, …)
- `covars` — 6 fns (`PlotChiSqCovar`, `PlotCorrelationsHeatmap`, `PlotInteractionEffects*`, `PlotSplitViolin`)
- `Covars` — 2 fns (`UnivariateRegressionTable`, `MultivariableRegressionTable`)
- `covariates` — 2 fns (`CreateNormativeTScoreModel`, `calculate_pathway_results`)

**Cluster D — grouping / comparison variable:**
`GroupVariable` (`PlotPValueComparisons`, `PlotSpiderChart`) · `GroupVar` (`Plot2GroupStats`) · `Group` (`PlotSplitViolin`) · `CompVariable` (`MakeComparisonTable`) · `TargetVar` (`CreateStatisticsTable`, `PlotZScore`) · `ClusterVar` (`PlotClusterBoxplot`) · `comparison_var` (`calculate_pathway_results`) · `interVar` (the 4 interaction functions) · `FacetVariables` (`MakeFacetCatComparisonTable`)

**Cluster E — merge/join keys:**
`Keys` (`CompareDatasets`, `ValidateMerge`) · `by` (`safe_merge`) · `MergeBy` (`Merge_ByClosestTime`) · `key` (`CombineCodebooks`) · `VariableCol` (`MergeCodebooks`, `CodebookMergeApp`)

**Cluster F — participant/record identifier:**
`ID` (`CreateRCIObject`, `ProjectRCI`, `PlotTimeSwimmer`) · `id` (`IQROutliers`) · `id_var` (`MergeFragmentedRecords`, `PlotSwimmerTransitions`, `SummarizeTransitions`) · `id_col` (`CreateSOMClusterModel`)

**Cluster G — codebook / variable-types object:**
`Codebook` (`InspectCategoricalSummary`, `PlotClusterBoxplot`, `PlotTimeSwimmer`, `PlotVolcanoEffects`, `UpdateCodebook`) · `codebook` (`CreateNormativeTScoreModel`) · `CB` (`AddToCodebook`) · `VariableCodebook` (`ReplaceMissingCode`) · `codebooks` (`MergeCodebooks`, `CodebookMergeApp`) · `OldCodebook`/`NewCodebook` (`CombineCodebooks`) · `VarTypes` (`ReValueFactors`, `RevalueData`) — the last two functions treat "codebook" and "variable types template" as different objects, but downstream they serve the same role.

**Cluster H — time/date column:**
`Time` (`PlotTimeSwimmer`) · `time_var` (`PlotSwimmerTransitions`, `SummarizeTransitions`) · `date_var` (same 2 fns + `MergeFragmentedRecords`) · `DateVariable` (`PlotTimeDistribution`) · `TimeVar1`/`TimeVar2` (`Merge_ByClosestTime`) · `time_var_before`/`time_var_add` (`safe_merge`)

**Cluster I — upstream result object (first arg of second-stage functions):**
`PCAObj` (`plotPCA`) vs `PCAObject` (`ExtractPCAComponentSummary`) vs `PCAInput` (`ProjectPCA`, third position) · `CompareObj` (`ExploreDatasetComparison`, `PlotDatasetComparison`) · `MergeObj` (`ExploreMergeValidation`, `PlotMergeValidation`) · `m` (`merge_detail`) · `merge_log` (`merge_summary_table`) · `res` (`add_r_and_stars`) · `results_table` (`PlotPathway_KT`) · `UnivariateRegressionTables` (`plotForestFromTable`) · `AnovaMatrixObject` (`plotSigAssociations`) · `CorrelationHeatmapObject` (`plotSigCorrelations`) · `normative_obj` (`ApplyNormativeTScores`) · `Object` (`ProjectRCI`) · `object` (`ProjectSOMCluster`)

**Cluster J — display/formatting knobs:**

- Digits: `numdecimals` (`CreateSummaryTable`, `MakeDataDictionary`, `FormattedDataDictionary`) vs `ValueDigits`/`pDigits`/`EffectSizeDigits` (comparison tables) vs `r_digits` (`add_r_and_stars`) vs `TooltipDigits` (`PlotSpiderChart`)
- Interactivity: `Interactive` (`PlotDatasetComparison`, `PlotMergeValidation`) vs `MakeInteractive` (`PlotSpiderChart`) vs `make_interactive` (`PlotSwimmerTransitions`)
- Point/font size: `point_size` (`Plot2GroupStats`, `PlotInteractionEffectsContinuous`) vs `PointSize` (`PlotSpiderChart`, `PlotTimeSwimmer`) vs `pSize` (`plotForestFromTable`); `BaseSize` (3 fns) vs `BaseFontSize` (`AssemblePlots`)
- Top-N: `TopN` (5 fns) vs `top_n` (`ExtractPCAComponentSummary`)
- **`Alpha` means two different things**: significance threshold in `PlotVolcanoEffects` (`Alpha = 0.05`) but transparency in `PlotTimeSwimmer` (`Alpha = 0.9`) / `alpha` in `PlotInteractionEffectsContinuous` (`0.6`) / `FillAlpha` in `PlotSpiderChart`. Meanwhile significance thresholds are elsewhere `Pthresh` (`plotSig*`), `FDRAlpha` (`MultivariableRegressionTable`), `label_q` (`Plot2GroupStats`).

**Consistent bright spots** (worth preserving as the de-facto standard): `Relabel` (25 uses, always this spelling, always `TRUE` default), `Ordinal` (16 uses, same spelling — though defaults diverge, see §3.3), `VariableCategories` (6 uses), `Parametric` (5 uses).

### 1.3 Parameter order inconsistency

- **Data-first is the norm** (~70 functions) but the projection family disagrees with itself: `ProjectRCI(Data, Object, ID)`, `ProjectZScore(df, variables, parameters, …)`, `ProjectPCA(Data, VarsToReduce = NULL, PCAInput, …)`, `ApplyNormativeTScores(df, normative_obj, …)` are all data-first, while **`ProjectSOMCluster(object, new_df, …)` is model-first**.
- **Required arguments after optional ones** (breaks positional calling and reads oddly in `args()`):
  - `CreateRCIObject(Data, Variables, DataFormat = c(...), ID, …)` — required `ID` after defaulted `DataFormat`
  - `ProjectPCA(Data, VarsToReduce = NULL, PCAInput, …)` — required `PCAInput` third, after a default
  - `PlotCatInteractionEffectsMatrix(Data, xVars, yVars = NULL, xVarLabels = NULL, yVarLabels = NULL, interVar)` — required `interVar` dead last (sibling `PlotNumInteractionEffectsMatrix` gives `interVar = NULL` a default instead)
  - `Plot2GroupStats(Data, Variables, VariableCategories = NULL, impClust, normalClust, GroupVar, …)` — three required args after a default
  - `MakeComparisonTable(DataFrame, CompVariable = NULL, Variables, ..., Covariates = NULL, …)` — required `Variables` after defaulted `CompVariable`, **and `...` placed mid-signature** (all other functions put `...` last), which forces every later argument to be fully named.
- The two swimmer functions take the same conceptual arguments in different dialects: `PlotTimeSwimmer(Data, ID, Time, State, …)` vs `PlotSwimmerTransitions(data, id_var, time_var, status_var, …)`.

---

## 2. Return value consistency

### 2.1 By family

| Family | Return shape | Consistent? |
|---|---|---|
| `Explore*` (2) | `htmltools::tagList()` dashboard | ✅ Yes |
| `CompareDatasets` / `ValidateMerge` | named list with `SummaryText`, `Summary`, checks | ✅ Structurally parallel |
| `safe_merge` / `merge_detail` / `merge_summary_table` | list(`data`,`validation`,`log`) / invisible input / tibble | ✅ Coherent pipeline |
| `Create*Object` (PCA, MCA, ZScore, MScore, RCI, SOM, NormativeTScore) | list "object" | ⚠️ Mostly — but only `ZScoreObj`, `MScoreObj`, `SciDataReportR_RCI`, and `Pipeline_SOMClust` get an S3 class; `CreatePCAObject`/`CreateMCAObject`/`CreateNormativeTScoreModel` return unclassed lists |
| Table makers | — | ❌ Five different types: gtsummary (`MakeTable1`, `MakeComparisonTable`, `MakeFacetCatComparisonTable`), kableExtra HTML (`CreateSummaryTable`, `CreateStatisticsTable`), gt (`FormattedDataDictionary`), plain data frame (`MakeDataDictionary`), list-with-`FormattedTable` (`UnivariateRegressionTable`, `MultivariableRegressionTable`) |
| `Plot*` (34 fns) | — | ❌ See below |

### 2.2 The `Plot*` split

- **Single ggplot** (~15): `PlotMissingData`, `PlotContinuousDistributions`, `PlotCategoricalDistributions`, `PlotClusterBoxplot`, `PlotTimeSwimmer`, `PlotTimeDistribution`, `PlotSplitViolin`, `PlotAssociations`, `PlotPValueComparisons`, `PlotInteractionEffectsContinuous`, `PlotZScore`, `PlotPathway_KT`, `plotForestFromTable`, `add_r_and_stars`, `geom_starcaption` (labs object)
- **List of plots + tables** (~15): `PlotCorrelationsHeatmap`, `PlotAnovaRelationshipsMatrix`, `PlotChiSqCovar`, `PlotPhiHeatmap`, `PlotPointCorrelationsHeatmap`, `PlotDirectionalHeatmaps`, `PlotMiningMatrix`, `PlotVolcanoEffects`, `Plot2GroupStats`, `PlotBlandAltman`, `PlotPartialRegressionScatter`, `Plot*InteractionEffectsMatrix` (×3), `plotSigAssociations`/`plotSigCorrelations` (lists of plots)
- **htmlwidget**: `plotPCA` (plotly, always)
- **Conditional shape**: `PlotDatasetComparison` and `PlotMergeValidation` (list if `Plot = "All"`, single otherwise; ggplot or plotly depending on `Interactive`), `PlotSpiderChart` (ggplot or plotly via `MakeInteractive`), `PlotSwimmerTransitions` (ggplot, plotly, or list via `return_data`)

Within the list-returning heatmap family, **element names diverge unnecessarily**: `PlotPhiHeatmap`/`PlotPointCorrelationsHeatmap`/`PlotDirectionalHeatmaps` use `Unadjusted`/`FDRCorrected`, but `PlotChiSqCovar` uses `p`/`pvaltable`/`p_FDR`/`pvaltable_FDR` and `PlotAnovaRelationshipsMatrix` uses `p`/`p_FDR`/`pvaltable`. Downstream, `plotSigAssociations` and `plotSigCorrelations` must be told the p-value column name (`PVar`) precisely because these shapes differ.

---

## 3. Documentation gaps

### 3.1 Coverage summary

- **`@param` coverage: complete.** Every exported function documents every argument (verified against `man/*.Rd`; apparent gaps were combined `@param a,b` tags). Watch item: several files (`FormattedDataDictionary.R`, `ReplaceMissingCode.R`, `UpdateCodebook.R`) have a stray blank line between the doc block and `#' @export` — roxygen currently tolerates it, but it's fragile.
- **`@return`: 107/108.** Missing only in `use_EDATemplate` (both copies of it — see §4.2).
- **`@examples`: 29/108 (27%).** The 79 without examples include the entire heatmap/matrix family (`PlotCorrelationsHeatmap`, `PlotMiningMatrix`, `PlotChiSqCovar`, …), all the codebook utilities except `AddToCodebook`/`MergeCodebooks`, both regression-table functions, and all compatibility aliases. Functions **with** examples: `AddToCodebook`, `ApplyNormativeTScores`, `AssemblePlots`, `CreateNormativeTScoreModel`, `CreatePCAObject`, `CreatePCATable`, `CreateVariableTypesTemplate`, `ExtractPCAComponentSummary`, `IQROutliers`, `InspectCategoricalSummary`, `MakeComparisonTable`, `MakeDataDictionary`, `MergeCodebooks`, `PlotClusterBoxplot`, `PlotInteractionEffectsMatrix`, `PlotMissingData`, `PlotPathway_KT`, `PlotSwimmerTransitions`, `PlotTimeSwimmer`, `PlotVolcanoEffects`, `PrepNumericData`, `SummarizeTransitions`, `calculate_pathway_results`, `merge_detail`, `merge_summary_table`, `plotForestFromTable`, `plotPCA`, `safe_merge`, `windsorize`.

### 3.2 Documented default ≠ actual default

1. **`use_EDATemplate`** — doc: *"default: `EDA_Report.qmd`"*. Actual default (from the winning duplicate definition, see §4.2): `"Scripts/script_GetStarted.R"`. Even the original definition's default is `"Reports/EDA_Report.qmd"`, which still doesn't match the doc text. `man/use_EDATemplate.Rd` currently shows the same `\usage` line **twice**.
2. **`plotSigAssociations`** — `@param PVar` says *"default is `"P"`"*; code default is `"p"` ([plotSigAssociations.R:14](R/plotSigAssociations.R)). Sibling `plotSigCorrelations` really does default to `"P"` — so doc and sibling agree, the code diverges (probably deliberately, matching the differing upstream column names — which is itself the §2.2 problem).
3. **`AddToCodebook`** — `@param VariableLabel` says *"default: same as VariableName"* while the formal default is `NA` (the substitution happens inside the body). Technically behavioral, but the stated default is not the literal default.
4. **`plotPCA`** — `t = "NULL"`: the default is the *string* `"NULL"`, not `NULL`. Documented accurately as deprecated, but worth flagging while renaming.

### 3.3 Same parameter, different defaults across functions (not doc errors, but trap-prone)

- `Ordinal` defaults to `TRUE` in 8 functions (`getCatVars`, `getBinaryVars`, `PlotCategoricalDistributions`, `PlotContinuousDistributions`, `PlotChiSqCovar`, `PlotDirectionalHeatmaps`, `PlotPointCorrelationsHeatmap`, `PlotZScore`) and `FALSE` in 8 others (`getNumVars`, `CreateMCAObject`, `CreatePCAObject`, `CreateSummaryTable`, `PlotAssociations`, `PlotAnovaRelationshipsMatrix`, `PlotCorrelationsHeatmap`, `PlotInteractionEffectsMatrix`). Some of that is principled (numeric vs categorical context) but it's undocumented as a convention.
- `minThresh` is **75 (percent)** in `CreateMCAObject` but **0.85 (proportion)** in `CreatePCAObject` — same name, different unit.
- `Standardize` defaults to `TRUE` in `MultivariableRegressionTable` but `FALSE` in `UnivariateRegressionTable`.
- `Interactive = TRUE` (`PlotMergeValidation`) vs `MakeInteractive = FALSE` (`PlotSpiderChart`) vs `make_interactive = FALSE` (`PlotSwimmerTransitions`).

---

## 4. Dead or overlapping functionality

### 4.1 Compatibility aliases (11) — intentional, but they double the API surface

`CalcMScore→CreateMScoreObject`, `CalcZScore→CreateZScoreObject`, `CreateMCATable→CreateMCAObject`, `CreateNormativeTScores→CreateNormativeTScoreModel`, `CreatePCATable→CreatePCAObject`, `CreatePathwayPlot_KT→PlotPathway_KT`, `CreateZScorePlot→PlotZScore`, `Make_DataDictionary→MakeDataDictionary`, `Pipeline_SOMClust→CreateSOMClusterModel`, `Project_SOMClust→ProjectSOMCluster`, `Project_ZScore→ProjectZScore`. All are documented `function(...)` pass-throughs. None emit a deprecation warning, so nothing nudges users toward the canonical name.

### 4.2 Genuine defects found in passing

- **`use_EDATemplate` is defined twice** — once in [use_EDATemplate.R](R/use_EDATemplate.R) (copies the EDA Quarto template) and again in [use_GetStartedScript.R:6](R/use_GetStartedScript.R) (a copy-paste that was never renamed to `use_GetStartedScript`; it copies `script_GetStarted.R` instead). Because the second file loads later alphabetically, **the exported `use_EDATemplate` currently copies the get-started script, not the EDA template**, and there is no `use_GetStartedScript` function at all despite the filename.
- **`CodebookMergeApp`** calls `br()`, `div()`, `h3()`, `h4()`, `observe()`, `observeEvent()` bare (most other shiny calls are `shiny::`-prefixed) while shiny is only in `Suggests` — these will fail at runtime unless the user has attached shiny ([MergeCodebooks.R:242](R/MergeCodebooks.R) ff.).
- **`%notin%`** is defined in [notInOperator.R:18](R/notInOperator.R) and listed as an `@aliases` of `%!in%`, but is **not exported** — the documented alias doesn't work from outside the package.
- **`ReLabelData`** ([RelabelData.R](R/RelabelData.R)) — defined, never exported, never called anywhere in `R/`: dead code. It also contains a misplaced `#' @importFrom` *inside* the function body, and note the casing collision-in-spirit with `RevalueData`/`ReValueFactors`.
- **Empty files**: [newfunction.R](R/newfunction.R) and [SampleFunction.R](R/SampleFunction.R) are 0 lines.
- No exported function references a nonexistent internal helper — `findGlobals` came back clean apart from the shiny symbols above (`everything`/`sym` resolve via the dplyr import).

### 4.3 Overlapping purpose

| Overlap | Notes |
|---|---|
| `PlotTimeSwimmer` vs `PlotSwimmerTransitions` | Two swimmer plots, overlapping purposes, entirely different argument dialects (PascalCase vs snake_case; `ID`/`Time`/`State` vs `id_var`/`time_var`/`status_var`). |
| `MakeDataDictionary` vs `FormattedDataDictionary` vs `CreateVariableTypesTemplate` vs `UpdateDataDictionary` | Same underlying "variable dictionary" concept; data frame vs gt table vs CSV template vs updater. Reasonable as a family, but names don't signal the relationship. |
| `MakeTable1` vs `CreateSummaryTable` | Both descriptive summary tables of a variable list (gtsummary vs kableExtra HTML). |
| Codebook suite: `AddToCodebook`, `UpdateCodebook`, `CombineCodebooks`, `MergeCodebooks`, `CodebookMergeApp` | `CombineCodebooks` (old vs new diff/merge) and `MergeCodebooks` (multi-codebook harmonization) are easy to confuse; `Merge` vs `Combine` semantics are not distinguished anywhere. |
| `windsorize` | Misspelled (standard term is *winsorize*), camelCase-lowercase, operates on a vector but calls its argument `Data`. |
| `RevalueData` vs `ReValueFactors` | Same verb, two casings (`Revalue` vs `ReValue`), closely related purposes. |
| `PlotCorrelationsHeatmap` / `PlotPointCorrelationsHeatmap` / `PlotPhiHeatmap` / `PlotDirectionalHeatmaps` / `PlotChiSqCovar` | Coherent family, but divergent return-element names (§2.2) prevent writing generic downstream code. |

---

## 5. Naming decisions to make

These are the choices that would let a follow-up standardization pass proceed mechanically:

1. **Primary data frame argument**: `Data`, `DataFrame`, or `df`? (Current split ≈ 31 / 18 / 12, plus `data`, `Dataframe`, `DatatoRevalue`.) This is the single highest-leverage decision.
2. **Package-wide parameter casing**: commit to PascalCase (the 61% incumbent) or snake_case (the tidyverse norm and the style of your newest functions)? Decide, too, whether low-level utilities (`safe_merge` family) are allowed a different dialect.
3. **Function-name casing**: keep PascalCase and rename the 12 camelCase (`getNumVars` → `GetNumVars`?) and 2 non-alias snake_case exports, or leave the helpers as-is?
4. **Variable-selection argument**: `Variables` everywhere, or keep specialized names (`VarsToReduce`, `xVars`/`yVars`, `OutcomeVars`/`PredictorVars`)? At minimum: unify `Variables` vs `variables`, and `outcomeVars` vs `OutcomeVars`.
5. **Covariates**: `Covariates`, `Covars`, or `covars`? (Pick one; four spellings today.)
6. **Grouping variable**: `GroupVariable`, `GroupVar`, `Group`, or `CompVariable`?
7. **Merge keys**: `Keys`, `by`, `MergeBy`, or `key`?
8. **ID column**: `ID`, `id`, `id_var`, or `id_col`?
9. **Codebook argument**: `Codebook` everywhere (retiring `CB`, `VariableCodebook`, `VarTypes`)?
10. **Upstream-object argument**: one convention for second-stage functions (e.g. always `Object`, or always the class name like `PCAObj`) — currently 15 different names.
11. **Interactivity flag**: `Interactive`, `MakeInteractive`, or `make_interactive` — and one default (`TRUE` vs `FALSE`).
12. **Digits**: `numdecimals` vs `*Digits` family.
13. **`Alpha`**: reserve for transparency (ggplot convention) and use `Pthresh`/`SigAlpha` for significance, or vice versa?
14. **`Ordinal` default**: document the convention (TRUE in categorical contexts, FALSE in numeric ones) or unify?
15. **`minThresh` units**: percent (75) or proportion (0.85) — MCA and PCA currently disagree.
16. **`Plot*` return contract**: (a) single plot vs list — formalize which functions return which, or add a uniform `ReturnData`-style argument; (b) standardize list element names (`Unadjusted`/`FDRCorrected` vs `p`/`p_FDR`/`pvaltable`).
17. **Aliases**: keep the 11 compatibility aliases silently, add `.Deprecated()` warnings, or drop at next major version?
18. **Housekeeping calls** (not naming, but blocking): fix the `use_EDATemplate` duplicate definition; export or un-document `%notin%`; namespace the bare shiny calls in `CodebookMergeApp`; delete `ReLabelData`, `newfunction.R`, `SampleFunction.R`; decide the fate of `windsorize`'s spelling and `CreatePCAObject`'s `ImputeMethod`/`imputeMethod` pair.
