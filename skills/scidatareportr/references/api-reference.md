# SciDataReportR API reference

Generated from `NAMESPACE` and `man/*.Rd` for SciDataReportR 21.4.0.
Run `Rscript tools/build_scidatareportr_skill_reference.R` after changing exports or Rd documentation.

This catalog has one entry per public export. Usage blocks omit arguments explicitly marked `lifecycle::deprecated()`; consult the compatibility guide for migration help.

## `%!in%`

**Purpose:** Negated "in" Operator

**Canonical usage**
```r
x \%!in\% y

x \%notin\% y
```

**Description:** This custom operator returns TRUE if the element on the left-hand side of the operator is not found in the vector on the right-hand side when using the %in% operator. Otherwise, it returns FALSE.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `%notin%`

**Arguments**
- `x`: The element to be tested for absence in the vector.
- `y`: The vector in which to search for the element.

**Returns:** TRUE if the element is not found in the vector, otherwise FALSE.

**See also:** None documented.

## `%notin%`

**Purpose:** Negated "in" Operator

**Canonical usage**
```r
x \%!in\% y

x \%notin\% y
```

**Description:** This custom operator returns TRUE if the element on the left-hand side of the operator is not found in the vector on the right-hand side when using the %in% operator. Otherwise, it returns FALSE.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `%!in%`

**Arguments**
- `x`: The element to be tested for absence in the vector.
- `y`: The vector in which to search for the element.

**Returns:** TRUE if the element is not found in the vector, otherwise FALSE.

**See also:** None documented.

## `add_biomarker_values`

**Purpose:** Add values to a biomarker performance heatmap

**Canonical usage**
```r
add_biomarker_values(
  plot,
  value_var = "HeatmapValue",
  digits = 2,
  size = 3,
  color = "black"
)
```

**Description:** Adds numeric cell labels to the heatmap returned by codelink[=ScreenBiomarkerPerformance]ScreenBiomarkerPerformance(). This helper is intended for downstream annotation of SciDataReportR biomarker heatmaps while keeping the default heatmap uncluttered and hover-friendly.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `plot`: A biomarker heatmap ggplot, typically codeScreenBiomarkerPerformance(...)$Plots$Heatmap.
- `value_var`: Column in codeplot$data used for labels. Default is code"HeatmapValue".
- `digits`: Number of decimal places. Default is code2.
- `size`: Text size. Default is code3.
- `color`: Text color. Default is code"black".

**Returns:** A ggplot with numeric values added to each non-missing heatmap cell.

**See also:** None documented.

## `add_delta_r_and_stars`

**Purpose:** Add DeltaR values and significance stars to a correlation comparison heatmap

**Canonical usage**
```r
add_delta_r_and_stars(
  res,
  star_from = c("fdr", "raw"),
  delta_digits = 2,
  delta_size = 3,
  star_size = 5,
  delta_color = "black",
  star_color = "black",
  delta_nudge_y = -0.18,
  star_nudge_y = 0.2,
  remove_existing_stars = TRUE
)
```

**Description:** Adds correlation differences and significance stars to a heatmap returned by codelink[=PlotCorrelationComparisons]PlotCorrelationComparisons(). This is the correlation-comparison counterpart to codelink[=add_r_and_stars]add_r_and_stars().

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `res`: An object returned by codelink[=PlotCorrelationComparisons]PlotCorrelationComparisons().
- `star_from`: Whether stars should use code"fdr" or code"raw" comparison p-values.
- `delta_digits`: Number of decimal places used for DeltaR.
- `delta_size`: Text size for DeltaR labels.
- `star_size`: Text size for significance stars.
- `delta_color`: Color for DeltaR labels.
- `star_color`: Color for significance stars.
- `delta_nudge_y`: Vertical position adjustment for DeltaR.
- `star_nudge_y`: Vertical position adjustment for stars.
- `remove_existing_stars`: Logical. Remove the original star-only layer before adding the combined annotations.

**Returns:** A ggplot containing DeltaR values and significance stars.

**See also:** None documented.

## `add_r_and_stars`

**Purpose:** Add r-values and significance stars to a correlations heatmap

**Canonical usage**
```r
add_r_and_stars(
  res,
  star_from = c("existing", "raw", "fdr", "P", "P_adj", "column"),
  star_col = NULL,
  r_var = "R",
  r_digits = 2,
  r_size = 3,
  r_color = "black",
  r_nudge_y = -0.28,
  star_size = 6,
  star_color = "black",
  star_nudge_y = 0.15,
  remove_existing_stars = TRUE,
  p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  p_labels = c("***", "**", "*", "")
)
```

**Description:** Adds correlation coefficients and significance stars to a heatmap object returned by codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap(). This is a downstream annotation helper for SciDataReportR correlation workflows, not a general-purpose plotting function.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `res`: The list returned by codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap(). The object should contain codeUnadjusted$plot and/or codeFDRCorrected$plot ggplot objects with plot data columns such as codeR, codeP, codeP_adj, codestars, or codestars_FDR.
- `star_from`: One of: "existing" (use whatever the chosen plot already mapped as stars), "raw" (use codestars or compute from codeP), "fdr" (use codestars_FDR or compute from codeP_adj), "P","P_adj" (compute from those columns), "column" (use codestar_col).
- `star_col`: Column name to use when star_from = "column"
- `r_var`: Column name for correlation values (default "R")
- `r_digits`: Decimal places for r labels
- `r_size, star_size`: Text sizes for r and stars
- `r_color, star_color`: Colors for r and stars
- `r_nudge_y, star_nudge_y`: Vertical nudges (r down, stars up)
- `remove_existing_stars`: If TRUE, remove any pre-existing star text layers
- `p_breaks, p_labels`: Cutpoints/labels for computing stars from p

**Returns:** A ggplot with r-values and stars added.

**See also:** None documented.

## `AddToCodebook`

**Purpose:** Add a new variable to a codebook

**Canonical usage**
```r
AddToCodebook(
  codebook,
  VariableName,
  VariableLabel = NA,
  VariableType = NA,
  VariableCategory = NA,
  VariableRecode = NA,
  VariableCode = NA,
  VariableExclude = NA,
  VariableNotes = NA,
  ...
)
```

**Description:** Adds one variable entry to a codebook while preserving its existing schema. In addition to the standard codebook fields, named values supplied through code... populate user-defined columns. A new named code... column is added to the codebook (with codeNA for existing rows) and produces a warning.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `codebook`: A data frame representing the codebook. It must contain a codeVariable column.
- `VariableName`: A single, non-missing, non-empty character variable name. It must not already appear in codecodebook$Variable.
- `VariableLabel`: A single label for the variable. Defaults to codeVariableName when codeNA.
- `VariableType`: A single value for the codeType column.
- `VariableCategory`: A single value for the codeCategory column.
- `VariableRecode`: A single value for the codeRecode column.
- `VariableCode`: A single value for the codeCode column.
- `VariableExclude`: A single value for the codeExclude column.
- `VariableNotes`: A single value for the codeNotes column.
- `CB`: strongDeprecated (since 19.15.0). Use codecodebook instead.
- `...`: Named, single atomic values for user-defined codebook columns. Names matching existing columns populate them. New names create a column and warn. Standard fields (codeVariable, codeLabel, codeType, codeCategory, codeRecode, codeCode, codeExclude, and codeNotes) must be supplied through their corresponding formal arguments.

**Returns:** A data frame representing the updated codebook with the new variable added.

**See also:** None documented.

## `ApplyFDRCorrection`

**Purpose:** Apply multiple-comparison correction across a p-value matrix

**Canonical usage**
```r
ApplyFDRCorrection(
  pmat,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  outcome_margin = 2,
  method = "fdr",
  outcome_ids = NULL,
  predictor_ids = NULL,
  symmetric = "auto",
  include_diagonal = FALSE
)
```

**Description:** Central helper used by the SciDataReportR heatmap/matrix family to apply multiple-comparison correction (FDR by default) with a selectable scope. All non-finite p-values (codeNA, codeNaN) are left untouched and excluded from the correction, matching the long-standing behavior of the plotting functions that now delegate to this helper.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `pmat`: A numeric matrix or data frame of p-values, or a plain numeric vector. Matrices and data frames keep their dimensions and dimnames.
- `fdr_scope`: Either code"matrix" (default), code"per_outcome", or code"per_predictor". code"matrix" corrects across all p-values at once (one family). code"per_outcome" corrects separately within each outcome's p-values: for matrix input, groups run along codeoutcome_margin; for vector input, groups are defined by codeoutcome_ids.
- `outcome_margin`: For matrix/data-frame input with codefdr_scope = "per_outcome": code2 (default) if outcomes are columns, code1 if outcomes are rows. Ignored for code"matrix" scope and for vector input.
- `method`: Correction method passed to codelink[stats:p.adjust]stats::p.adjust(). Default code"fdr" (Benjamini-Hochberg).
- `outcome_ids`: Optional vector (same length as codepmat) identifying the outcome each p-value belongs to. Only used - and then required - when codepmat is a vector and codefdr_scope = "per_outcome". This is how the long-format table functions (for example codelink[=PlotPhiHeatmap]PlotPhiHeatmap() or codelink[=PlotChiSqCovar]PlotChiSqCovar()) group their p-values by outcome.
- `predictor_ids`: Optional vector (same length as codepmat) identifying the predictor each p-value belongs to. Only used - and then required - when codepmat is a vector and codefdr_scope = "per_predictor".
- `symmetric`: How to treat a square matrix whose two triangles hold the same p-values. code"auto" (default) detects symmetry and corrects each pair once; codeTRUE requires it (and errors if codepmat is not symmetric); codeFALSE restores the whole-matrix behavior that counts each pair twice. Ignored for vector input. See the Symmetric matrices section.
- `include_diagonal`: Logical; for symmetric input, whether the diagonal (self-comparisons) joins the family being corrected. Default codeFALSE, which excludes it and returns codeNA on the diagonal.

**Returns:** An object of the same shape as codepmat (matrix, data frame, or vector) containing adjusted p-values. Non-finite entries remain codeNA.

**See also:** None documented.

## `ApplyNormativeTScores`

**Purpose:** Apply a normative T-score model to new data

**Canonical usage**
```r
ApplyNormativeTScores(
  data,
  normative_obj,
  score_prefix = "Norm"
)
```

**Description:** Applies a previously fitted normative regression model to new data and computes predicted values, z-scores, and T-scores using the same preprocessing settings used during model development.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the test variable, count variable, and all predictors required by the normative model.
- `normative_obj`: A list returned by codeCreateNormativeTScoreModel().
- `score_prefix`: A character string prefix used when naming output columns. Defaults to code"Norm".
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A tibble containing the original data plus scored columns, each named using the codescore_prefix string as a prefix: describe itemverbscore_prefixRawThe raw input score. itemverbscore_prefixScaledThe transformed analysis-scale score. itemverbscore_prefixPredictedThe predicted score from the normative model. itemverbscore_prefixZThe z-score. itemverbscore_prefixTThe T-score.

**See also:** None documented.

## `AssemblePlots`

**Purpose:** Assemble ggplot objects into a unified multi-panel figure

**Canonical usage**
```r
AssemblePlots(
  Plots,
  ncol = NULL,
  nrow = NULL,
  AutoLayout = TRUE,
  RemoveNULL = TRUE,
  CollectLegend = TRUE,
  LegendPosition = "bottom",
  LegendRelativeSize = 0.08,
  Theme = ggplot2::theme_minimal(),
  BaseFontSize = 12,
  GlobalTheme = NULL,
  GlobalLayers = NULL,
  RemoveTitles = FALSE,
  UseNamesAsTitles = FALSE,
  Align = "hv",
  Axis = "tblr",
  Labels = NULL,
  LabelSize = 14,
  SuggestedBaseWidth = 4,
  SuggestedBaseHeight = 4,
  ReturnMetadata = FALSE
)
```

**Description:** A SciDataReportR wrapper around codelink[cowplot:plot_grid]cowplot::plot_grid() and codelink[cowplot:get_legend]cowplot::get_legend(). Cowplot supplies the underlying grid layout, alignment, panel-label, and shared-legend tools; codeAssemblePlots() packages them into one reporting-oriented workflow.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `Plots`: A ggplot object or list of ggplot objects.
- `ncol`: Optional number of columns.
- `nrow`: Optional number of rows.
- `AutoLayout`: Logical; automatically estimate layout when codencol and codenrow are not supplied.
- `RemoveNULL`: Logical; remove NULL plots before assembly.
- `CollectLegend`: Logical; combine legends into a shared legend.
- `LegendPosition`: Position of shared legend. One of code"top", code"bottom", code"left", code"right", or code"none".
- `LegendRelativeSize`: Relative size allocated to legend area.
- `Theme`: Optional global ggplot theme.
- `BaseFontSize`: Base font size applied globally.
- `GlobalTheme`: Optional additional theme applied globally.
- `GlobalLayers`: Optional list of ggplot layers/scales applied to all plots.
- `RemoveTitles`: Logical; remove plot titles globally.
- `UseNamesAsTitles`: Logical; use plot list names as titles when plot titles are missing.
- `Align`: Plot alignment passed to cowplot.
- `Axis`: Axis alignment passed to cowplot.
- `Labels`: Optional panel labels.
- `LabelSize`: Panel label font size.
- `SuggestedBaseWidth`: Base width per column used for suggested figure dimensions.
- `SuggestedBaseHeight`: Base height per row used for suggested figure dimensions.
- `ReturnMetadata`: Logical; if codeTRUE, returns plot plus metadata.

**Returns:** If codeReturnMetadata = FALSE, returns a ggplot object. If codeReturnMetadata = TRUE, returns a list containing: describe itemPlotCombined plot object itemnrowEstimated number of rows itemncolEstimated number of columns itemSuggestedWidthSuggested figure width itemSuggestedHeightSuggested figure height itemNumPlotsNumber of plots

**See also:** None documented.

## `CalcMScore`

**Purpose:** Calculate robust M-scores for numeric variables

**Canonical usage**
```r
CreateMScoreObject(
  data,
  variables = NULL,
  names_prefix = "M_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE,
  constant = 1.4826
)

CalcMScore(...)
```

**Description:** Calculate median/MAD-based M-scores for selected numeric variables and return both the transformed data and the parameters needed to review or reuse the transformation. codeCalcMScore() has been superseded by codeCreateMScoreObject(). It remains available as a backwards-compatible alias and returns the same reusable M-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateMScoreObject`

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of numeric variables to transform. If codeNULL, numeric variables are detected with codelink[=getNumVars]getNumVars().
- `names_prefix`: Prefix for generated M-score columns.
- `RetainLabels`: Logical; keep variable labels when possible.
- `RenameLabels`: Logical; rename generated labels when labels are retained.
- `center`: Logical; subtract the median before scaling.
- `scale`: Logical; divide by the median absolute deviation.
- `constant`: Scaling constant passed to codelink[stats:mad]stats::mad().
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateMScoreObject]CreateMScoreObject().

**Returns:** An object of class code"MScoreObj", a list with: itemize item codeMScores: data frame of M-score variables only item codeDataWithM: original codedf plus M-score variables item codeParameters: data frame with codeVariable, codeN, codeMedian, and codeMAD item codeCenter: logical flag used item codeScale: logical flag used item codeConstant: MAD scaling constant used

**See also:** None documented.

## `calculate_pathway_results`

**Purpose:** Calculate Pathway Results for Metabolite Comparisons

**Canonical usage**
```r
calculate_pathway_results(
  data,
  comparison_var,
  covariates = NULL,
  metabolites,
  comparison_type = "auto",
  use_point_correlation = FALSE
)
```

**Description:** Calculate Pathway Results for Metabolite Comparisons

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing metabolite and comparison data
- `comparison_var`: Character string specifying the comparison variable name
- `covariates`: Character vector of covariate names (optional)
- `metabolites`: Character vector of metabolite names to analyze
- `comparison_type`: Character string: "auto", "binary", or "continuous"
- `use_point_correlation`: Logical, if TRUE uses point correlation for binary comparisons
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** Data frame with metabolite results including fold change or correlation values

**See also:** None documented.

## `CalcZScore`

**Purpose:** Calculate Z-scores (or standardized scores) and return data + parameters

**Canonical usage**
```r
CreateZScoreObject(
  data,
  variables = NULL,
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)

CalcZScore(...)
```

**Description:** Standardizes each variable to a common scale and, critically, returns the constants used to do it so the identical transformation can be replayed on other data later. codeCalcZScore() has been superseded by codeCreateZScoreObject(). It remains available as a backwards-compatible alias and returns the same reusable Z-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateZScoreObject`

**Arguments**
- `data`: Data frame with variables to standardize.
- `variables`: Character vector of variable names. If NULL, uses SciDataReportR::getNumVars(df).
- `names_prefix`: Prefix to prepend to variable names (default "Z_").
- `RetainLabels`: Logical; if TRUE and Hmisc is available, copy labels.
- `RenameLabels`: Logical; if TRUE, apply the same prefix to labels.
- `center`: Logical; if TRUE, subtract the mean.
- `scale`: Logical; if TRUE, divide by the SD.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateZScoreObject]CreateZScoreObject().

**Returns:** An object of class "ZScoreObj", a list with: itemize item ZScores: data frame of standardized variables only item DataWithZ: original df + standardized variables item Parameters: data frame with Variable, N, Mean, SD item Center: logical flag used item Scale: logical flag used

**See also:** codelink[=ProjectZScore]ProjectZScore() to apply stored parameters to new data, and codelink[=CreateNormativeTScoreModel]CreateNormativeTScoreModel() when the reference values should also be adjusted for covariates such as age or education.

## `CodebookMergeApp`

**Purpose:** Interactive codebook harmonization dashboard

**Canonical usage**
```r
CodebookMergeApp(
  codebooks,
  VariableCol = "Variable",
  auto_type_mapping = TRUE,
  ignore_columns = NULL
)
```

**Description:** Launch a Shiny dashboard for reviewing and harmonizing multiple codebooks before deterministic merging with MergeCodebooks().

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `codebooks`: Named list of codebook data frames.
- `VariableCol`: Name of variable identifier column.
- `auto_type_mapping`: Logical; normalize common type synonyms.
- `ignore_columns`: Optional metadata columns to ignore.

**Returns:** Launches a Shiny app.

**See also:** None documented.

## `CombineCodebooks`

**Purpose:** Combine Two Codebooks with Conflict Detection

**Canonical usage**
```r
CombineCodebooks(
  OldCodebook,
  NewCodebook,
  keys = "Variable"
)
```

**Description:** codeCombineCodebooks compares two codebook data frames (an old version and a new version), identifies added or removed variables and columns, detects cell-by-cell differences, and produces a combined codebook. For records without any differences, it returns a single "Combined" row; for records with differences, it returns both the "Old" and "New" rows, flagged as conflicts.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `OldCodebook`: A data.frame or tibble representing the old codebook. Each row must correspond to a single variable entry.
- `NewCodebook`: A data.frame or tibble representing the new codebook. Must have the same structure (column names) as codeOldCodebook, though extra or missing columns will be handled.
- `keys`: A string giving the name of the key column to identify variables (e.g. "Variable"). Defaults to "Variable".
- `key`: strongDeprecated (since 19.15.0). Use codekeys instead.

**Returns:** A list with elements: itemize item codeadded_variables: character vector of keys present only in codeNewCodebook. item coderemoved_variables: character vector of keys present only in codeOldCodebook. item codecolumns_added: character vector of column names present only in codeNewCodebook. item codecolumns_removed: character vector of column names present only in codeOldCodebook. item codevalue_differences: tibble of cell-level differences (codeRowID, codeField, codeOldValue, codeNewValue). item codecombined_df: tibble containing the merged codebook with versions and conflict flag.

**See also:** None documented.

## `CompareDatasets`

**Purpose:** Compare two versions of a dataset

**Canonical usage**
```r

```

**Description:** Compare an old dataset and a new dataset using one or more key variables. This function identifies record-level, variable-level, and cell-level changes between dataset versions. It is useful when reviewing updated data extracts, revised REDCap exports, cleaned spreadsheet versions, or vendor-delivered dataset updates.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `OldData`: A data frame representing the earlier dataset version.
- `NewData`: A data frame representing the newer dataset version.
- `keys`: Character vector of key variables used to align records across the two datasets. Multiple keys are supported, such as codec("study_id", "TimePoint").
- `Keys`: strongDeprecated (since 19.15.0). Use codekeys instead.

**Returns:** A list with dataset comparison results, including: describe itemSummaryTextA plain-text summary of the dataset comparison. itemSummaryOne-row tibble with core comparison metrics. itemFingerprintTibble comparing rows, columns, and unique key combinations. itemKeyTypesTibble showing key variable classes before coercion. itemChecksTibble summarizing comparison checks and pass/warning/fail status. itemStructureChangesTibble of variables added to or removed from NewData, using normalized variable names. itemAddedVariablesTibble of variables present in NewData but not OldData, using normalized variable names. itemRemovedVariablesTibble of variables present in OldData but not NewData, using normalized variable names. itemAddedRecordsTibble of key combinations present in NewData but not OldData. itemRemovedRecordsTibble of key combinations present in OldData but not NewData. itemDuplicateKeysList containing duplicated key rows from OldData and NewData. itemComparisonKeysList describing matching keys, compared keys, and keys excluded from cell comparison due to duplicate key combinations. itemNameRepairAuditTibble describing variables whose raw names differ after removing tibble-style code...number suffixes. itemComparisonVariableMapTibble mapping normalized variable names to the raw OldData and NewData names used for cell-level comparison. itemClassAuditTibble comparing variable classes for common non-key variables. itemModifiedValuesLong-format tibble of cell-level value changes. itemVariableChangeSummaryTibble summarizing changes by variable. itemTopChangedVariablesTop changed variables by number of modified values. itemSuspiciousChangesTibble of high-change-rate or class-change variables.

**See also:** None documented.

## `ConvertOrdinalToNumeric`

**Purpose:** Prepare ordinal variables for analysis

**Canonical usage**
```r
ConvertOrdinalToNumeric(
  data,
  variables = NULL,
  TreatOrdinalAs = c("Continuous", "Categorical", "Both", "Exclude"),
  Relabel = TRUE,
  ReturnMetadata = FALSE
)
```

**Description:** Applies a consistent ordinal-treatment policy to selected variables. Ordinal score mappings recorded by codelink[=RevalueData]RevalueData() are used when available; otherwise ordered-factor ranks are used.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame containing the variables.
- `variables`: Character vector of variables to consider. If codeNULL, all columns are considered.
- `TreatOrdinalAs`: How ordinal variables are handled: code"Continuous", code"Categorical", code"Both", or code"Exclude".
- `Relabel`: Logical; when codeTreatOrdinalAs = "Both", apply descriptive labels to the derived categorical and continuous variables.
- `ReturnMetadata`: Logical; if codeFALSE (default), return only the transformed data frame. If codeTRUE, return a list containing the data, selected variables, ordinal variables, variable map, and treatment.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A transformed data frame, or a metadata list when codeReturnMetadata = TRUE.

**See also:** None documented.

## `createBinaryMapping`

**Purpose:** Create a Mapping Table for Binary Variables

**Canonical usage**
```r
createBinaryMapping(
  data,
  CatVars,
  prefer = NULL
)
```

**Description:** Identifies binary variables and returns a deterministic mapping for 0/1 coding. itemize item Factors with emphexplicit order (ordered = TRUE) do NOT use heuristics; the highest (last) level is Positive. item Logicals map to Negative = "FALSE", Positive = "TRUE". item Numeric 0/1 (or any 2-value numeric) maps Positive to the numeric maximum. item Characters / unordered factors use minimal heuristics (no race/PWH/sex terms).

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe.
- `CatVars`: Character vector of candidate binary variables.
- `prefer`: Optional named character vector of explicit positive levels, e.g., c(STATUS = "PWH", Smoker = "Yes"). This overrides other rules.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A data.frame with columns: Variable, Label, PositiveLevel, NegativeLevel.

**See also:** codelink[=getBinaryVars]getBinaryVars() to find the candidates.

## `CreateClusterModel_Gower_PAM`

**Purpose:** Fit a projectable Gower-distance PAM model for mixed clinical data

**Canonical usage**
```r
CreateClusterModel_Gower_PAM(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for mixed continuous, binary, ordinal, and nominal clinical measures where medoid exemplars are more interpretable than centroids.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric, logical, factor, or ordered variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory"/code"explore" or code"finalize".
- `k_range`: Candidate medoid counts in exploratory mode.
- `final_k`: Final medoid count in finalized mode.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed.
- `stability_resamples`: Number of 90% participant subsample refits.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large or high-dimensional analyses.

**Returns:** A fitted PAM model with frozen medoids, numeric ranges, categorical levels, silhouette and medoid-distance metrics in codeModelInfo$fit_table, and subsample stability. Mean silhouette is higher-is-better and mean assigned-medoid Gower distance is lower-is-better (a dissimilarity, not a probability). Figures sit beside what they describe: codefit_plot reviews candidates; codeModelInfo$plots holds codesilhouette (the per-participant profile from the selected PAM solution), codesilhouette_by_k, codegower_map, codecategorical_composition, codecategorical_composition_by_cluster, codecategorical_enrichment, and codeprofiles; codeModelInfo$FitDiagnostics$plots holds the medoid-distance histogram; codeProbFit$plots holds assignment-margin figures; and codeStability$plots holds cluster-recovery and complementary stability diagnostics.

**See also:** None documented.

## `CreateClusterModel_HDBSCAN`

**Purpose:** Fit a projectable HDBSCAN model

**Canonical usage**
```r
CreateClusterModel_HDBSCAN(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  minPts_range = 2:10,
  cluster_selection_epsilon_range = c(0, 0.05, 0.1),
  final_minPts = NULL,
  final_cluster_selection_epsilon = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for irregularly shaped numeric clusters, variable density, and data where a meaningful noise/outlier group is expected.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `minPts_range`: Candidate minimum-points settings in exploratory mode.
- `cluster_selection_epsilon_range`: Candidate extraction epsilon values.
- `final_minPts, final_cluster_selection_epsilon`: Finalized density settings.
- `ZScoreType`: Frozen numeric preprocessing. codeScaling is a compatibility alias.
- `Scaling`: Compatibility alias for codeZScoreType.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A fitted HDBSCAN model with cluster/noise assignments, membership probabilities, outlier scores, frozen nearest-core support thresholds, and subsample ARI/Jaccard and noise-recovery metrics. Persistence is higher-is-better, noise proportion is lower-is-better, and the extracted class count is data-derived; membership probability and outlier score are assignment diagnostics rather than candidate fit metrics. Figures sit beside what they describe: codefit_plot reviews the density grid; codeModelInfo$plots holds codedensity_review, codepersistence, codecluster_persistence, and codeprofiles; codeModelInfo$FitDiagnostics$plots holds the nearest-core-distance histogram and codeoutlier_score_by_cluster, both measured against the same frozen reference a projected case is triaged on; codeProbFit$plots holds membership-probability figures; and codeStability$plots holds subsample agreement, per-cluster recovery, and noise recovery.

**See also:** None documented.

## `CreateClusterModel_KMeans`

**Purpose:** Fit a projectable K-means clustering model

**Canonical usage**
```r
CreateClusterModel_KMeans(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  nstart = 50L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for approximately spherical, similarly sized groups in a numeric feature space; use scaling unless all variables are commensurate.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `k_range`: Candidate cluster counts in exploratory mode.
- `final_k`: Number of clusters for a finalized K-means solution.
- `ZScoreType`: Frozen numeric preprocessing. codeScaling is a compatibility alias.
- `Scaling`: Compatibility alias for codeZScoreType.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `nstart`: Number of random K-means starts.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A fitted model with codeModelInfo$fit_table (WSS, between-cluster sum of squares, silhouette, Calinski-Harabasz index, and minimum cluster size) and subsample stability. WSS is lower-is-better; between-cluster sum of squares, silhouette, and Calinski-Harabasz are higher-is-better. Figures sit beside what they describe: codefit_plot reviews candidates; codeModelInfo$plots holds codeelbow, codesilhouette (the per-participant silhouette profile of the selected solution), codesilhouette_by_k, codecalinski_harabasz, codecentre_heatmap, codecentre_profile, and codeprofiles; codeModelInfo$FitDiagnostics$plots holds the distance-to-centroid histogram; codeProbFit$plots holds assignment-margin figures; and codeStability$plots holds cluster-recovery and complementary stability diagnostics. codeModelInfo$ReviewSpace freezes the two-dimensional review space shared by the training and projection maps; it is diagnostic only and does not affect clustering.

**See also:** None documented.

## `CreateClusterModel_LatentClass`

**Purpose:** Fit a projectable latent class model for categorical measures

**Canonical usage**
```r
CreateClusterModel_LatentClass(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ClusterVariableName = "Cluster",
  nrep = 20L,
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Use latent class analysis when the clustering variables are nominal or ordinal questionnaire, symptom, diagnosis, or assay-call items. It estimates class-specific response probabilities and posterior membership.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing categorical variables.
- `variables`: Variables included in the latent class model.
- `method`: Either code"exploratory"/code"explore" or code"finalize".
- `k_range`: Candidate class counts.
- `final_k`: Optional class count; when supplied, fit only this solution.
- `ClusterVariableName`: Output class column name.
- `nrep`: Number of random starts per candidate.
- `seed`: Random seed controlling latent-class random starts.
- `stability_resamples`: Number of 90% participant subsample refits.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large or high-dimensional analyses.

**Returns:** A fitted latent-class model with BIC, AIC, log likelihood, entropy, class-size and subsample stability metrics in codeModelInfo$fit_table. Log likelihood and entropy are higher-is-better; AIC and BIC are lower-is-better. Figures sit beside what they describe: codefit_plot reviews candidates; codeModelInfo$plots holds coderesponse_probabilities, codeitem_profiles, codebic, codeentropy, codeposterior_map, and codecategorical_composition; codeProbFit$plots holds posterior-confidence figures; and codeStability$plots holds cluster-recovery and complementary recovery.

**See also:** None documented.

## `CreateClusterModel_MCA_MClust`

**Purpose:** Fit MCA followed by Mclust for nominal categorical data

**Canonical usage**
```r
CreateClusterModel_MCA_MClust(
  data,
  variables,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  models = c(1L, 2L, 3L, 6L),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  mca_variance_threshold = 75,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Use this pipeline for many nominal categorical variables when a lower-dimensional category space is useful before mixture clustering.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `k_range`: Candidate cluster counts in exploratory mode.
- `models`: Numeric tidyLPA mclust model IDs. The supported models are: code1 (EEI), equal variance and zero covariance; code2 (VVI), varying variance and zero covariance; code3 (EEE), equal variance and equal covariance; and code6 (VVV), varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are estimated separately for each cluster. Models 4 and 5 require OpenMx and are not supported by these pipelines.
- `final_k, final_model`: Finalized cluster count and numeric model ID.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `mca_variance_threshold`: Cumulative MCA inertia percentage retained.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A fitted MCA plus Mclust model containing the frozen MCA object, mixture-model fit and subsample stability metrics, and posterior probabilities. codeModelInfo$plots leads with the reduction layer's codescree and codeloadings, followed by the mixture-model structure figures in MCA score space and codecategorical_composition restated on the original items. BIC/ICL/AIC, entropy, and uncertainty use the same interpretation as Mclust.

**See also:** None documented.

## `CreateClusterModel_MClust`

**Purpose:** Fit a projectable Gaussian-mixture clustering model

**Canonical usage**
```r
CreateClusterModel_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  models = c(1L, 2L, 3L, 6L),
  final_k = NULL,
  final_model = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for continuous measures when clinically meaningful groups may differ in means, variance, or covariance and posterior uncertainty is useful.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `k_range`: Candidate cluster counts in exploratory mode.
- `models`: Numeric tidyLPA mclust model IDs. The supported models are: code1 (EEI), equal variance and zero covariance; code2 (VVI), varying variance and zero covariance; code3 (EEE), equal variance and equal covariance; and code6 (VVV), varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are estimated separately for each cluster. Models 4 and 5 require OpenMx and are not supported by these pipelines.
- `final_k, final_model`: Finalized cluster count and numeric model ID.
- `ZScoreType`: Frozen numeric preprocessing. codeScaling is a compatibility alias.
- `Scaling`: Compatibility alias for codeZScoreType.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A projectable mixture model. codeModelInfo$fit_table contains BIC, ICL, entropy, uncertainty, and subsample stability metrics. codeBIC and codeICL are reported on the scale codemclust returns them on, where higher is better; codeAIC is reported on the conventional scale, where lower is better. Entropy is higher-is-better classification separation and maximum uncertainty is lower-is-better. Shared stability fields are defined in the clustering reference vignette. codeModelInfo$AHP the advisory recommendation. Figures sit beside what they describe: codefit_plot reviews candidates; codeModelInfo$plots holds codebic, codeicl, codeentropy, codecentre_heatmap, codecentre_profile, and codeprofiles; codeModelInfo$FitDiagnostics$plots holds the Mahalanobis distance histogram; codeProbFit$plots holds posterior-confidence figures; and codeStability$plots holds cluster-recovery and complementary stability diagnostics.

**See also:** None documented.

## `CreateClusterModel_PCA_KMeans`

**Purpose:** Fit PCA followed by K-means

**Canonical usage**
```r
CreateClusterModel_PCA_KMeans(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  final_k = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  nstart = 50L,
  pca_variance_threshold = 0.85,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for high-dimensional correlated continuous measures with compact clusters in PCA score space.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `k_range`: Candidate cluster counts in exploratory mode.
- `final_k`: Number of clusters for a finalized PCA + K-means solution.
- `ZScoreType`: Frozen numeric preprocessing. codeScaling is a compatibility alias.
- `Scaling`: Compatibility alias for codeZScoreType.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `nstart`: Number of random K-means starts.
- `pca_variance_threshold`: Cumulative variance retained by the existing PCA workflow.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A frozen PCA and K-means pipeline with full-pipeline subsample stability. codeModelInfo$plots leads with the reduction layer's codescree and codeloadings, followed by the K-means structure figures in score space (including the per-participant codesilhouette profile) and codeprofiles restated in the original measurement scale. WSS/BSS, silhouette, and Calinski-Harabasz use the same interpretation as K-means.

**See also:** None documented.

## `CreateClusterModel_PCA_MClust`

**Purpose:** Fit PCA followed by Mclust

**Canonical usage**
```r
CreateClusterModel_PCA_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  k_range = 2:10,
  models = c(1L, 2L, 3L, 6L),
  final_k = NULL,
  final_model = NULL,
  ZScoreType = NULL,
  Scaling = NULL,
  ClusterVariableName = "Cluster",
  seed = 93421L,
  pca_variance_threshold = 0.85,
  stability_resamples = 0L,
  stability_seed = seed + 1L,
  stability_progress = FALSE,
  stability_cores = NULL
)
```

**Description:** Best for correlated continuous measures when clustering a lower-dimensional, frozen PCA representation is preferable to raw features.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing numeric clustering variables.
- `variables`: Variables used for clustering.
- `method`: Either code"exploratory" or code"finalize".
- `k_range`: Candidate cluster counts in exploratory mode.
- `models`: Numeric tidyLPA mclust model IDs. The supported models are: code1 (EEI), equal variance and zero covariance; code2 (VVI), varying variance and zero covariance; code3 (EEE), equal variance and equal covariance; and code6 (VVV), varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are estimated separately for each cluster. Models 4 and 5 require OpenMx and are not supported by these pipelines.
- `final_k, final_model`: Finalized cluster count and numeric model ID.
- `ZScoreType`: Frozen numeric preprocessing. codeScaling is a compatibility alias.
- `Scaling`: Compatibility alias for codeZScoreType.
- `ClusterVariableName`: Output cluster column name.
- `seed`: Random seed retained for reproducibility.
- `pca_variance_threshold`: Cumulative variance retained by the existing PCA workflow.
- `stability_resamples`: Number of 90% participant subsample refits used to estimate candidate reproducibility. Subsamples are drawn without replacement. Use code0 to disable stability analysis.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.

**Returns:** A frozen PCA and Gaussian-mixture pipeline with full-pipeline subsample stability. codeModelInfo$plots leads with the reduction layer's codescree and codeloadings, followed by the mixture-model structure figures in score space and codeprofiles restated in the original measurement scale. Candidate metrics and the advisory recommendation are in codeModelInfo$fit_table and codeModelInfo$AHP; BIC/ICL/AIC, entropy, and uncertainty use the same interpretation as Mclust.

**See also:** None documented.

## `CreateClusterModel_SOM_HDBSCAN`

**Purpose:** Fit HDBSCAN clusters on a frozen self-organizing map

**Canonical usage**
```r
CreateClusterModel_SOM_HDBSCAN(
  data,
  variables = NULL,
  method = c("exploratory", "finalize"),
  minPts_range = 2:10,
  cluster_selection_epsilon_range = c(0, 0.05, 0.1),
  final_minPts = NULL,
  final_cluster_selection_epsilon = NULL,
  ClusterVariableName = "Cluster",
  seed_som = 934521L,
  seed_hdbscan = 93421L,
  stability_resamples = 0L,
  stability_seed = seed_hdbscan + 1L,
  stability_progress = FALSE,
  stability_cores = NULL,
  ...
)
```

**Description:** Trains the same frozen SOM used by codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust(), then applies HDBSCAN to its node codebook. Participants inherit the cluster (or noise) label of their best-matching node. This is a node-based phenotype model: projected participants are mapped to the original nodes and are not refit with HDBSCAN.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame used to train the SOM and node-level HDBSCAN model.
- `variables`: Numeric variables used for SOM training.
- `method`: Either code"exploratory" or code"finalize".
- `minPts_range`: Candidate HDBSCAN minimum-point settings for SOM nodes.
- `cluster_selection_epsilon_range`: Candidate HDBSCAN extraction epsilon settings.
- `final_minPts`: Optional finalized HDBSCAN minimum-point setting.
- `final_cluster_selection_epsilon`: Optional finalized extraction epsilon.
- `ClusterVariableName`: Name of the appended cluster column.
- `seed_som`: Seed used for SOM training.
- `seed_hdbscan`: Seed retained in the model specification.
- `stability_resamples`: Number of 90% participant subsample refits used for final-model stability. Subsamples are drawn without replacement and reuse the reference model's resolved SOM grid dimensions.
- `stability_seed`: Seed controlling participant subsampling.
- `stability_progress`: Whether to print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.
- `...`: Additional arguments passed to codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Returns:** A codePipeline_SOM_HDBSCAN object containing frozen SOM and HDBSCAN node models, assignments, diagnostics, and a projection specification. Persistence and minimum node-cluster size are higher-is-better, noise proportion is lower-is-better, and extracted class count is data-derived.

**See also:** None documented.

## `CreateClusterModel_SOM_MClust`

**Purpose:** SOM + latent profile clustering pipeline (with AHP and distance baselines)

**Canonical usage**
```r
CreateClusterModel_SOM_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize", "explore"),
  k_range = 2:10,
  models = c(1, 2, 3, 6),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  seed_som = 934521L,
  seed_lpa = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_var = NULL,
  lpa_progress = FALSE,
  lpa_em_itmax = 100L,
  lpa_em_tol = 1e-05,
  lpa_timeout_seconds = 120,
  lpa_drop_zero_sd = TRUE,
  lpa_zero_sd_tol = 1e-08,
  skip_model_after_n_failures = 2L,
  slow_fit_seconds = 120,
  min_nodes_per_cluster = 5,
  high_dist_quantile = 0.95,
  low_prob_threshold = 0.7,
  stability_resamples = 0L,
  stability_seed = 934522L,
  stability_progress = FALSE,
  stability_cores = NULL,
  .NodeClusterFn = NULL
)

Pipeline_SOM_MClust(...)

Pipeline_SOMClust(...)

CreateSOMClusterModel(...)
```

**Description:** End-to-end pipeline to: itemize item Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied Z-score object. item Fit a Self-Organizing Map (SOM; kohonen) on complete cases. item Generate aweSOM visualizations (Circular, Line, Cloud) with optional relabeling using variable labels from the original data frame. item Cluster SOM codebook vectors using latent profile analysis (tidyLPA / mclust backend). item In codemethod = "exploratory", fit a grid of models and select a recommended solution using an Analytic Hierarchy Process (AHP)-style index combining AIC, BIC, and Entropy. item In codemethod = "finalize", fit a user-specified model and number of profiles. item Map node-level clusters and posterior probabilities back to individuals. item Store training variable summaries used later to quantify whether projected cohorts fall outside the original training range. This supports a train once, project many clinical phenotyping workflow: learn phenotype structure in a training cohort, then project new cohorts into the fixed phenotype space without reclustering. Ideal use: correlated continuous clinical or biomarker measures where a topology-preserving map is clinically informative before model-based profiles. Missing data: itemize item SOM and clustering are fit only on rows with complete Z-scores. item The returned codeDataWithClusters has exactly the original rows and columns plus one cluster column; rows not used in SOM/LPA get NA. item The returned codeProbFit$individual is also full length, preserving one row per input row with NA posterior probabilities for rows excluded from SOM/LPA. Z-score behavior: itemize item codeZScoreType = "Center and Scale"/"Center Only"/"Scale Only" computes Z-scores from codedf via codeCreateZScoreObject(). item codeZScoreType = "ZScoreObj" projects Z-scores using an external codeZScoreObj via codeProjectZScore(). item codeZScoreType = "PreZScored" uses existing Z-score columns in codedf as-is and does not re-zscore. Readable SOM + Mclust workflow wrapper for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust(). Compatibility wrapper for codelink[=Pipeline_SOM_MClust]Pipeline_SOM_MClust(). Deprecated alias for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `Pipeline_SOM_MClust`, `Pipeline_SOMClust`, `CreateSOMClusterModel`

**Arguments**
- `data`: Data frame containing the variables to be used in SOM and clustering.
- `variables`: Optional character vector of variable names. If NULL, numeric variables are auto-detected using codeSciDataReportR::getNumVars(df, Ordinal = FALSE). In codeZScoreType = "PreZScored", this can also be NULL if you supply codeZScoreVars or if Z-score columns can be auto-detected by prefix.
- `method`: One of code"exploratory" (default) or code"finalize". In code"exploratory", a grid of models is fit and AHP chooses the recommended solution. In code"finalize", the user must specify codefinal_k and codefinal_model.
- `k_range`: Integer vector of numbers of clusters/profiles to consider in exploratory mode. Default code2:10.
- `models`: Integer vector of model specifications for tidyLPA's mclust backend. Model 1 uses equal variance and zero covariance; model 2 uses varying variance and zero covariance; model 3 uses equal variance and equal covariance; and model 6 uses varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are cluster-specific. Supported values and the default are codec(1, 2, 3, 6). Models 4 and 5 require OpenMx and are intentionally unsupported.
- `final_k`: Integer; number of profiles for codemethod = "finalize".
- `final_model`: Integer; model specification for codemethod = "finalize" (should be one of codemodels).
- `ClusterVariableName`: Name of the cluster column in the output. Defaults to code"Cluster". If this column already exists in codedf, it is overwritten (with a message).
- `ZScoreType`: One of: itemize item code"Center and Scale" (default) item code"Center Only" item code"Scale Only" item code"ZScoreObj" (use an existing ZScore object) item code"PreZScored" (use existing Z-score columns in df as-is)
- `ZScoreObject`: Optional ZScoreObj (from codeCreateZScoreObject() or codeProjectZScore()) to use when codeZScoreType = "ZScoreObj".
- `som_xdim, som_ydim`: Optional integers for SOM grid dimensions. If NULL, a square grid with side length codeceiling(n_complete^(1/3)) is used.
- `som_topo`: SOM topology for codekohonen::somgrid(), default code"hexagonal".
- `som_neigh`: SOM neighbourhood function, default code"gaussian".
- `seed_som, seed_lpa`: Integer seeds for SOM and LPA steps (defaults 934521 and 93421).
- `Relabel`: Logical; if TRUE (default), aweSOM plots are relabeled using variable labels from the emphoriginal codedf (via Hmisc or sjlabelled when available) by stripping the Z-score prefix.
- `ZScorePrefix`: Character prefix used for Z-score columns when codeZScoreType = "PreZScored". Default code"Z_".
- `ZScoreVars`: Optional character vector of Z-score column names to use when codeZScoreType = "PreZScored". If NULL, the function attempts to infer them from codevariables or by detecting columns starting with codeZScorePrefix.
- `id_var`: Optional character scalar. If provided and present in codedf, this column is carried into codeProbFit$individual for convenience.
- `lpa_progress`: Logical; if TRUE, print short progress messages while fitting model/profile combinations.
- `lpa_em_itmax`: Integer; maximum number of EM iterations passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_em_tol`: Numeric; EM convergence tolerance passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_timeout_seconds`: Optional timeout in seconds for individual LPA fits. Use NULL to disable timeouts.
- `lpa_drop_zero_sd`: Logical; if TRUE, remove SOM code dimensions with near-zero standard deviation before LPA.
- `lpa_zero_sd_tol`: Numeric tolerance used when codelpa_drop_zero_sd = TRUE.
- `skip_model_after_n_failures`: Optional integer; skip a model family after this many failures.
- `slow_fit_seconds`: Optional runtime threshold used to flag slow LPA fits in diagnostics.
- `min_nodes_per_cluster`: Optional minimum average SOM nodes per cluster considered before attempting a candidate profile count.
- `high_dist_quantile`: Numeric value between 0 and 1 used to define high SOM-distance flags from the training distance distribution. Default is code0.95.
- `low_prob_threshold`: Numeric posterior probability threshold used to flag uncertain phenotype membership. Default is code0.70.
- `stability_resamples`: Number of 90% participant subsample refits used to assess reproducibility for every successful exploratory candidate. Subsamples are drawn without replacement and reuse the reference model's resolved SOM grid dimensions. Defaults to code0 (disabled); use code50 for an exploratory stability screen.
- `stability_seed`: Integer seed for participant subsampling.
- `stability_progress`: Logical; if TRUE, print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `id_col`: strongDeprecated (since 19.15.0). Use codeid_var instead.
- `.NodeClusterFn`: Internal. A function taking the SOM codebook matrix and returning a list with a codenode_cluster integer vector (one label per SOM node) and, optionally, codefit_table, codeahp_best_row, coderecommendation, codebest_fit_name, and codefit_plot. When supplied, the SOM codebook is clustered by that function and the latent-profile grid is not fitted. Used by codelink[=CreateClusterModel_SOM_HDBSCAN]CreateClusterModel_SOM_HDBSCAN(); not part of the user-facing API.
- `...`: Arguments passed to codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Returns:** A list of class code"Pipeline_SOM_MClust" with components: itemize item codemethod, codevars_used, codeZScoreType, codeZScoreObject, codeZScoreVars, codeClusterVariableName item codeDataWithClusters: original codedf with only the cluster column appended item codefit_plot: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model (plus reproducibility when subsample stability is enabled) item codeModelInfo_SOM: list with codesom_model, codesom_codes, codesom_grid, codetraining_variable_summary, codeSOMFit (distance diagnostics, baselines, and per-cluster flags), codeplots (aweSOM plots) item codeModelInfo_MClust: list with codelpa_models, codefit_table, codeAHP information, and codediagnostics for LPA warnings, failures, runtimes, and preprocessing item codeModelInfo_MClust$Stability: subsample replicate, cluster recovery, and summary tables when codestability_resamples > 0 item codeProbFit: list with codenode (node-level posterior probabilities), codeindividual (full-length per-person mapping and probabilities), and probability plots

**See also:** None documented.

## `createFacetLabels`

**Purpose:** Create facet labels for ggplot2 based on variable labels in a data frame

**Canonical usage**
```r

```

**Description:** This function takes a data frame containing variable labels and creates facet labels suitable for use with ggplot2 facet functions.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing variable labels.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A character vector containing facet labels.

**See also:** codelink[=RevalueData]RevalueData(), which attaches the labels this reads, and codelink[=PlotContinuousDistributions]PlotContinuousDistributions(), which applies the same idea internally through its codeFacetLabelStyle argument.

## `CreateMCAObject`

**Purpose:** Create a reusable MCA object and visualizations

**Canonical usage**
```r
CreateMCAObject(
  data,
  VarsToReduce,
  VariableCategories = NULL,
  minThresh = 75,
  scale = TRUE,
  center = TRUE,
  Relabel = TRUE,
  Ordinal = FALSE,
  numComponents = NULL,
  ImputeMissing = FALSE
)

CreateMCATable(...)
```

**Description:** This function performs Multiple Correspondence Analysis (MCA) on a set of categorical variables, imputes missing data if needed, and generates a set of visualizations and tables to interpret the results. codeCreateMCATable() has been superseded by codeCreateMCAObject(). It remains available as a backwards-compatible alias and returns the same reusable MCA object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateMCATable`

**Arguments**
- `data`: A dataframe containing the data to be analyzed.
- `VarsToReduce`: A vector of column names in codeData to be included in the MCA.
- `VariableCategories`: An optional vector to assign specific categories to the variables in codeVarsToReduce. These will be used to color the loadings plot.
- `minThresh`: A numeric value representing the minimum cumulative variance threshold to determine the number of components. Default is 75%.
- `scale`: Logical, indicating whether to scale the variables. Default is TRUE.
- `center`: Logical, indicating whether to center the variables. Default is TRUE.
- `Relabel`: Logical, if TRUE, the function will replace missing labels in the data using an external helper function codeReplaceMissingLabels. Default is TRUE.
- `Ordinal`: Logical, if TRUE, the function will treat variables as ordinal for MCA. Default is FALSE.
- `numComponents`: An optional integer specifying the number of components to retain. If NULL, the number of components will be determined based on codeminThresh.
- `ImputeMissing`: Logical, if TRUE, missing values will be imputed using codemissRanger. Default is FALSE.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateMCAObject]CreateMCAObject().

**Returns:** A list with the following elements: itemp_screeA codeggplot object representing the scree plot, showing the cumulative and percentage of variance explained by each component. itempcaresultsThe MCA results object, which includes component scores and contributions. itemLoadingTableA data frame with the variable loadings for each component. itemScoresA data frame with the MCA scores for each individual in the data. itemCombinedDataThe original data combined with the MCA scores. itemLollipopA codeggplot object showing a lollipop plot of variable loadings across components.

**See also:** None documented.

## `CreateMCATable`

**Purpose:** Create a reusable MCA object and visualizations

**Canonical usage**
```r
CreateMCAObject(
  data,
  VarsToReduce,
  VariableCategories = NULL,
  minThresh = 75,
  scale = TRUE,
  center = TRUE,
  Relabel = TRUE,
  Ordinal = FALSE,
  numComponents = NULL,
  ImputeMissing = FALSE
)

CreateMCATable(...)
```

**Description:** This function performs Multiple Correspondence Analysis (MCA) on a set of categorical variables, imputes missing data if needed, and generates a set of visualizations and tables to interpret the results. codeCreateMCATable() has been superseded by codeCreateMCAObject(). It remains available as a backwards-compatible alias and returns the same reusable MCA object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateMCAObject`

**Arguments**
- `data`: A dataframe containing the data to be analyzed.
- `VarsToReduce`: A vector of column names in codeData to be included in the MCA.
- `VariableCategories`: An optional vector to assign specific categories to the variables in codeVarsToReduce. These will be used to color the loadings plot.
- `minThresh`: A numeric value representing the minimum cumulative variance threshold to determine the number of components. Default is 75%.
- `scale`: Logical, indicating whether to scale the variables. Default is TRUE.
- `center`: Logical, indicating whether to center the variables. Default is TRUE.
- `Relabel`: Logical, if TRUE, the function will replace missing labels in the data using an external helper function codeReplaceMissingLabels. Default is TRUE.
- `Ordinal`: Logical, if TRUE, the function will treat variables as ordinal for MCA. Default is FALSE.
- `numComponents`: An optional integer specifying the number of components to retain. If NULL, the number of components will be determined based on codeminThresh.
- `ImputeMissing`: Logical, if TRUE, missing values will be imputed using codemissRanger. Default is FALSE.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateMCAObject]CreateMCAObject().

**Returns:** A list with the following elements: itemp_screeA codeggplot object representing the scree plot, showing the cumulative and percentage of variance explained by each component. itempcaresultsThe MCA results object, which includes component scores and contributions. itemLoadingTableA data frame with the variable loadings for each component. itemScoresA data frame with the MCA scores for each individual in the data. itemCombinedDataThe original data combined with the MCA scores. itemLollipopA codeggplot object showing a lollipop plot of variable loadings across components.

**See also:** None documented.

## `CreateMScoreObject`

**Purpose:** Calculate robust M-scores for numeric variables

**Canonical usage**
```r
CreateMScoreObject(
  data,
  variables = NULL,
  names_prefix = "M_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE,
  constant = 1.4826
)

CalcMScore(...)
```

**Description:** Calculate median/MAD-based M-scores for selected numeric variables and return both the transformed data and the parameters needed to review or reuse the transformation. codeCalcMScore() has been superseded by codeCreateMScoreObject(). It remains available as a backwards-compatible alias and returns the same reusable M-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CalcMScore`

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of numeric variables to transform. If codeNULL, numeric variables are detected with codelink[=getNumVars]getNumVars().
- `names_prefix`: Prefix for generated M-score columns.
- `RetainLabels`: Logical; keep variable labels when possible.
- `RenameLabels`: Logical; rename generated labels when labels are retained.
- `center`: Logical; subtract the median before scaling.
- `scale`: Logical; divide by the median absolute deviation.
- `constant`: Scaling constant passed to codelink[stats:mad]stats::mad().
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateMScoreObject]CreateMScoreObject().

**Returns:** An object of class code"MScoreObj", a list with: itemize item codeMScores: data frame of M-score variables only item codeDataWithM: original codedf plus M-score variables item codeParameters: data frame with codeVariable, codeN, codeMedian, and codeMAD item codeCenter: logical flag used item codeScale: logical flag used item codeConstant: MAD scaling constant used

**See also:** None documented.

## `CreateNormativeTScoreModel`

**Purpose:** Create normative T-scores from a regression model

**Canonical usage**
```r
CreateNormativeTScoreModel(
  data,
  test_var,
  count_var,
  covariates,
  reference_var,
  reference_value,
  include_practice_effect = FALSE,
  baseline_count_value = 1,
  reverse_score = FALSE,
  convert_seconds = FALSE,
  seconds_divisor = 1000,
  log_transform = TRUE,
  codebook = NULL,
  return_plots = TRUE
)

CreateNormativeTScores(...)
```

**Description:** Fits a normative regression model in a user-defined reference subgroup and uses the model residual standard deviation to convert observed scores into z-scores and T-scores. This is useful for creating demographically adjusted cognitive norms with optional practice effect adjustment and optional preprocessing such as unit conversion, log transformation, and reverse scoring. codeCreateNormativeTScores() has been superseded by codeCreateNormativeTScoreModel(). It remains available as a backwards-compatible alias and returns the same reusable normative model.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateNormativeTScores`

**Arguments**
- `data`: A data frame containing the test variable, count variable, reference group variable, and covariates.
- `test_var`: A character string naming the raw test score variable.
- `count_var`: A character string naming the visit count or practice count variable.
- `covariates`: A character vector of covariate column names to include in the normative model.
- `reference_var`: A character string naming the variable used to define the normative reference group.
- `reference_value`: The value of codereference_var that defines the normative reference group.
- `include_practice_effect`: Logical. If codeTRUE, codecount_var is included as a predictor and the model is fit using all available visits in the reference group. If codeFALSE, the model is fit only on rows where codecount_var == baseline_count_value.
- `baseline_count_value`: The value of codecount_var used to define the baseline visit when codeinclude_practice_effect = FALSE. Defaults to code1.
- `reverse_score`: Logical. If codeTRUE, the analysis-scale score is multiplied by code-1 so that higher values reflect better performance.
- `convert_seconds`: Logical. If codeTRUE, the raw score is divided by codeseconds_divisor before further processing.
- `seconds_divisor`: Numeric divisor used when codeconvert_seconds = TRUE. Defaults to code1000.
- `log_transform`: Logical. If codeTRUE, applies codelog10() to the analysis score after optional unit conversion and before optional reverse scoring.
- `codebook`: Optional data frame with columns codeVariable and codeLabel. If supplied, plot labels use variable labels when available.
- `return_plots`: Logical. If codeTRUE, returns a list of diagnostic plots.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateNormativeTScoreModel]CreateNormativeTScoreModel().

**Returns:** A list with the following elements: describe itemdataA tibble containing the original data plus codeNormRaw, codeNormScaled, codeNormPredicted, codeNormZ, and codeNormT. itemmodelThe fitted codelm object. itemmodel_summaryA tibble of coefficient estimates. itemmodel_fitA one-row tibble containing model fit statistics. itemtraining_dataThe rows used to fit the normative model. itemplotsA named list of ggplot objects when codereturn_plots = TRUE. itemsettingsA list of preprocessing and modeling settings used.

**See also:** None documented.

## `CreateNormativeTScores`

**Purpose:** Create normative T-scores from a regression model

**Canonical usage**
```r
CreateNormativeTScoreModel(
  data,
  test_var,
  count_var,
  covariates,
  reference_var,
  reference_value,
  include_practice_effect = FALSE,
  baseline_count_value = 1,
  reverse_score = FALSE,
  convert_seconds = FALSE,
  seconds_divisor = 1000,
  log_transform = TRUE,
  codebook = NULL,
  return_plots = TRUE
)

CreateNormativeTScores(...)
```

**Description:** Fits a normative regression model in a user-defined reference subgroup and uses the model residual standard deviation to convert observed scores into z-scores and T-scores. This is useful for creating demographically adjusted cognitive norms with optional practice effect adjustment and optional preprocessing such as unit conversion, log transformation, and reverse scoring. codeCreateNormativeTScores() has been superseded by codeCreateNormativeTScoreModel(). It remains available as a backwards-compatible alias and returns the same reusable normative model.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateNormativeTScoreModel`

**Arguments**
- `data`: A data frame containing the test variable, count variable, reference group variable, and covariates.
- `test_var`: A character string naming the raw test score variable.
- `count_var`: A character string naming the visit count or practice count variable.
- `covariates`: A character vector of covariate column names to include in the normative model.
- `reference_var`: A character string naming the variable used to define the normative reference group.
- `reference_value`: The value of codereference_var that defines the normative reference group.
- `include_practice_effect`: Logical. If codeTRUE, codecount_var is included as a predictor and the model is fit using all available visits in the reference group. If codeFALSE, the model is fit only on rows where codecount_var == baseline_count_value.
- `baseline_count_value`: The value of codecount_var used to define the baseline visit when codeinclude_practice_effect = FALSE. Defaults to code1.
- `reverse_score`: Logical. If codeTRUE, the analysis-scale score is multiplied by code-1 so that higher values reflect better performance.
- `convert_seconds`: Logical. If codeTRUE, the raw score is divided by codeseconds_divisor before further processing.
- `seconds_divisor`: Numeric divisor used when codeconvert_seconds = TRUE. Defaults to code1000.
- `log_transform`: Logical. If codeTRUE, applies codelog10() to the analysis score after optional unit conversion and before optional reverse scoring.
- `codebook`: Optional data frame with columns codeVariable and codeLabel. If supplied, plot labels use variable labels when available.
- `return_plots`: Logical. If codeTRUE, returns a list of diagnostic plots.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateNormativeTScoreModel]CreateNormativeTScoreModel().

**Returns:** A list with the following elements: describe itemdataA tibble containing the original data plus codeNormRaw, codeNormScaled, codeNormPredicted, codeNormZ, and codeNormT. itemmodelThe fitted codelm object. itemmodel_summaryA tibble of coefficient estimates. itemmodel_fitA one-row tibble containing model fit statistics. itemtraining_dataThe rows used to fit the normative model. itemplotsA named list of ggplot objects when codereturn_plots = TRUE. itemsettingsA list of preprocessing and modeling settings used.

**See also:** None documented.

## `CreatePathwayPlot_KT`

**Purpose:** Plot the kynurenine-tryptophan pathway

**Canonical usage**
```r
PlotPathway_KT(
  results_table,
  title = "",
  value_type = "auto",
  metabolite_mapping = NULL,
  use_fdr = FALSE
)

CreatePathwayPlot_KT(...)
```

**Description:** Creates a pathway diagram for the kynurenine-tryptophan metabolic pathway with color-coded fold changes or correlations and significance indicators. codeCreatePathwayPlot_KT() has been superseded by codePlotPathway_KT(). It remains available as a backwards-compatible alias during the pathway plot's planned transition to a metabolomics-focused package.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotPathway_KT`

**Arguments**
- `results_table`: Data frame with columns: Metabolite, p_value, p_adj, and either "% Change" or "correlation"
- `title`: Character string for plot title
- `value_type`: Character string: "auto", "fold_change", or "correlation"
- `metabolite_mapping`: Named character vector mapping results table names to standard names. For example: c("N'-Formylkynurenine" = "N-Formylkynurenine", "Quinolinic Acid(log10)" = "Quinolinic Acid")
- `use_fdr`: Logical: if TRUE uses FDR-adjusted p-values (p_adj) for significance, if FALSE uses raw p-values. Default is FALSE.
- `...`: Arguments passed to codelink[=PlotPathway_KT]PlotPathway_KT().

**Returns:** A ggplot2 object

**See also:** None documented.

## `CreatePCAObject`

**Purpose:** Create a reusable PCA object and visualizations

**Canonical usage**
```r
CreatePCAObject(
  data,
  VarsToReduce,
  VariableCategories = NULL,
  Relabel = TRUE,
  minThresh = 0.85,
  scale = TRUE,
  center = TRUE,
  Ordinal = FALSE,
  numComponents = NULL,
  Mode = c("classic", "omics"),
  backend = c("psych", "prcomp", "irlba"),
  rotate = c("varimax", "none"),
  maxComponents = 20,
  maxScreeComponents = 20,
  VarianceFilter = NULL,
  VarianceFilterMethod = c("top_n", "variance_quantile"),
  MissingnessWarningThreshold = 0.2,
  ParticipantMissingnessWarningThreshold = 0.2,
  MissingDataStrategy = c("complete_cases", "impute", "stop"),
  ImputeMethod = c("missRanger", "median"),
  MaxMissingForImputation = 0.2,
  ImputeRowsWithAllMissing = FALSE,
  SuppressWarnings = FALSE
)

CreatePCATable(...)
```

**Description:** Perform principal component analysis (PCA) on specified variables and return reusable PCA results, scores, loading tables, combined data, and plots. codeCreatePCATable() has been superseded by codeCreatePCAObject(). It remains available as a backwards-compatible alias and returns the same reusable PCA object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreatePCATable`

**Arguments**
- `data`: A data frame containing the variables for PCA.
- `VarsToReduce`: Character vector of raw variable names to include in PCA.
- `VariableCategories`: Optional vector used to color variables in the lollipop loading plot. If supplied, it should be the same length and order as codeVarsToReduce.
- `Relabel`: Logical. If codeTRUE, variable labels are used in output tables and plots when available. If codeFALSE, raw variable names are used. Default is codeTRUE.
- `minThresh`: Numeric threshold for cumulative variance used to select the number of components when codenumComponents = NULL. Default is code0.85.
- `scale`: Logical. If codeTRUE, variables are scaled by their standard deviation before PCA. Default is codeTRUE.
- `center`: Logical. If codeTRUE, variables are centered before PCA. Default is codeTRUE.
- `Ordinal`: Logical. Placeholder retained for backward compatibility. Currently not used inside this function.
- `numComponents`: Optional integer number of components to retain. If codeNULL, the number is selected using codeminThresh.
- `Mode`: Character. Either code"classic" or code"omics". code"classic" preserves the original codepsych::principal() behavior. code"omics" uses capped PCA logic for higher-dimensional data.
- `backend`: Character. PCA backend to use. Options are code"psych", code"prcomp", and code"irlba". Default is code"psych" for compatibility.
- `rotate`: Character. Rotation method. Options are code"varimax" and code"none". Default is code"varimax".
- `maxComponents`: Integer. Maximum number of final components to retain in code"omics" mode when codenumComponents = NULL. Default is code20.
- `maxScreeComponents`: Integer. Maximum number of components used to estimate and plot the scree curve in code"omics" mode. Default is code20.
- `VarianceFilter`: Optional numeric value for variance filtering before PCA in code"omics" mode. If codeVarianceFilterMethod = "top_n", keeps the top codeVarianceFilter most variable variables. If codeVarianceFilterMethod = "variance_quantile", keeps variables with variance at or above the specified quantile.
- `VarianceFilterMethod`: Character. Either code"top_n" or code"variance_quantile". Default is code"top_n".
- `MissingnessWarningThreshold`: Numeric threshold for warning about variable-level missingness. Default is code0.20.
- `ParticipantMissingnessWarningThreshold`: Numeric threshold for warning about participant-level missingness. Default is code0.20.
- `MissingDataStrategy`: Character. Missing-data handling strategy. Options are code"complete_cases", code"impute", and code"stop". Default is code"complete_cases", which fits PCA only on rows complete across the final PCA variables and returns codeNA PCA scores for incomplete rows.
- `ImputeMethod`: Character. Missing-data imputation method used only when codeMissingDataStrategy = "impute". Options are code"missRanger" and code"median". Default is code"missRanger".
- `MaxMissingForImputation`: Numeric value between 0 and 1. When codeMissingDataStrategy = "impute", rows with missingness greater than this value across the final PCA variables are excluded from PCA scoring and receive codeNA component scores. Default is code0.20, meaning rows can be imputed if at least 80 percent of PCA variables are observed.
- `ImputeRowsWithAllMissing`: Logical. If codeFALSE, rows with 100 percent missingness across PCA variables are not imputed even if codeMaxMissingForImputation = 1. Default is codeFALSE.
- `SuppressWarnings`: Logical. If codeTRUE, suppresses PCA-specific warning messages from this function. Default is codeFALSE.
- `imputeMethod`: strongDeprecated (since 19.15.0). Use codeImputeMethod instead. If supplied, it sets codeImputeMethod and uses codeMissingDataStrategy = "impute" unless codeMissingDataStrategy was explicitly set.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreatePCAObject]CreatePCAObject().

**Returns:** A list with the following elements: describe itemp_screeA ggplot scree and cumulative variance plot. itempcaresultsThe PCA result object. In code"classic" mode this is a codepsych::principal() result after codepsych::fa.sort(). In code"omics" mode with code"prcomp" or code"irlba", this is a harmonized list containing loadings, scores, variance information, backend, mode, and rotation. itemLoadingTableA data frame of component loadings with raw variable names, labels, and duplicate-safe plot labels. itemScoresA data frame of component scores aligned to the original row order of codeData. Rows not used for PCA scoring receive codeNA scores. itemCombinedDataThe original input data with component scores appended. itemLollipopA ggplot lollipop loading plot. itemScaleParamsA list containing centering and scaling parameters. itemVarsUsedCharacter vector of variables actually used in PCA after preprocessing. itemVarianceTableA variance table used for optional omics variance filtering, or codeNULL. itemPreprocessingA list documenting variables dropped during preprocessing, missingness summaries, rows used for PCA, rows excluded from PCA, and rows imputed. itemModeThe PCA mode used. itemBackendThe PCA backend used. itemCenterLogical flag indicating whether centering was used. itemScaleLogical flag indicating whether scaling was used.

**See also:** None documented.

## `CreatePCATable`

**Purpose:** Create a reusable PCA object and visualizations

**Canonical usage**
```r
CreatePCAObject(
  data,
  VarsToReduce,
  VariableCategories = NULL,
  Relabel = TRUE,
  minThresh = 0.85,
  scale = TRUE,
  center = TRUE,
  Ordinal = FALSE,
  numComponents = NULL,
  Mode = c("classic", "omics"),
  backend = c("psych", "prcomp", "irlba"),
  rotate = c("varimax", "none"),
  maxComponents = 20,
  maxScreeComponents = 20,
  VarianceFilter = NULL,
  VarianceFilterMethod = c("top_n", "variance_quantile"),
  MissingnessWarningThreshold = 0.2,
  ParticipantMissingnessWarningThreshold = 0.2,
  MissingDataStrategy = c("complete_cases", "impute", "stop"),
  ImputeMethod = c("missRanger", "median"),
  MaxMissingForImputation = 0.2,
  ImputeRowsWithAllMissing = FALSE,
  SuppressWarnings = FALSE
)

CreatePCATable(...)
```

**Description:** Perform principal component analysis (PCA) on specified variables and return reusable PCA results, scores, loading tables, combined data, and plots. codeCreatePCATable() has been superseded by codeCreatePCAObject(). It remains available as a backwards-compatible alias and returns the same reusable PCA object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreatePCAObject`

**Arguments**
- `data`: A data frame containing the variables for PCA.
- `VarsToReduce`: Character vector of raw variable names to include in PCA.
- `VariableCategories`: Optional vector used to color variables in the lollipop loading plot. If supplied, it should be the same length and order as codeVarsToReduce.
- `Relabel`: Logical. If codeTRUE, variable labels are used in output tables and plots when available. If codeFALSE, raw variable names are used. Default is codeTRUE.
- `minThresh`: Numeric threshold for cumulative variance used to select the number of components when codenumComponents = NULL. Default is code0.85.
- `scale`: Logical. If codeTRUE, variables are scaled by their standard deviation before PCA. Default is codeTRUE.
- `center`: Logical. If codeTRUE, variables are centered before PCA. Default is codeTRUE.
- `Ordinal`: Logical. Placeholder retained for backward compatibility. Currently not used inside this function.
- `numComponents`: Optional integer number of components to retain. If codeNULL, the number is selected using codeminThresh.
- `Mode`: Character. Either code"classic" or code"omics". code"classic" preserves the original codepsych::principal() behavior. code"omics" uses capped PCA logic for higher-dimensional data.
- `backend`: Character. PCA backend to use. Options are code"psych", code"prcomp", and code"irlba". Default is code"psych" for compatibility.
- `rotate`: Character. Rotation method. Options are code"varimax" and code"none". Default is code"varimax".
- `maxComponents`: Integer. Maximum number of final components to retain in code"omics" mode when codenumComponents = NULL. Default is code20.
- `maxScreeComponents`: Integer. Maximum number of components used to estimate and plot the scree curve in code"omics" mode. Default is code20.
- `VarianceFilter`: Optional numeric value for variance filtering before PCA in code"omics" mode. If codeVarianceFilterMethod = "top_n", keeps the top codeVarianceFilter most variable variables. If codeVarianceFilterMethod = "variance_quantile", keeps variables with variance at or above the specified quantile.
- `VarianceFilterMethod`: Character. Either code"top_n" or code"variance_quantile". Default is code"top_n".
- `MissingnessWarningThreshold`: Numeric threshold for warning about variable-level missingness. Default is code0.20.
- `ParticipantMissingnessWarningThreshold`: Numeric threshold for warning about participant-level missingness. Default is code0.20.
- `MissingDataStrategy`: Character. Missing-data handling strategy. Options are code"complete_cases", code"impute", and code"stop". Default is code"complete_cases", which fits PCA only on rows complete across the final PCA variables and returns codeNA PCA scores for incomplete rows.
- `ImputeMethod`: Character. Missing-data imputation method used only when codeMissingDataStrategy = "impute". Options are code"missRanger" and code"median". Default is code"missRanger".
- `MaxMissingForImputation`: Numeric value between 0 and 1. When codeMissingDataStrategy = "impute", rows with missingness greater than this value across the final PCA variables are excluded from PCA scoring and receive codeNA component scores. Default is code0.20, meaning rows can be imputed if at least 80 percent of PCA variables are observed.
- `ImputeRowsWithAllMissing`: Logical. If codeFALSE, rows with 100 percent missingness across PCA variables are not imputed even if codeMaxMissingForImputation = 1. Default is codeFALSE.
- `SuppressWarnings`: Logical. If codeTRUE, suppresses PCA-specific warning messages from this function. Default is codeFALSE.
- `imputeMethod`: strongDeprecated (since 19.15.0). Use codeImputeMethod instead. If supplied, it sets codeImputeMethod and uses codeMissingDataStrategy = "impute" unless codeMissingDataStrategy was explicitly set.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreatePCAObject]CreatePCAObject().

**Returns:** A list with the following elements: describe itemp_screeA ggplot scree and cumulative variance plot. itempcaresultsThe PCA result object. In code"classic" mode this is a codepsych::principal() result after codepsych::fa.sort(). In code"omics" mode with code"prcomp" or code"irlba", this is a harmonized list containing loadings, scores, variance information, backend, mode, and rotation. itemLoadingTableA data frame of component loadings with raw variable names, labels, and duplicate-safe plot labels. itemScoresA data frame of component scores aligned to the original row order of codeData. Rows not used for PCA scoring receive codeNA scores. itemCombinedDataThe original input data with component scores appended. itemLollipopA ggplot lollipop loading plot. itemScaleParamsA list containing centering and scaling parameters. itemVarsUsedCharacter vector of variables actually used in PCA after preprocessing. itemVarianceTableA variance table used for optional omics variance filtering, or codeNULL. itemPreprocessingA list documenting variables dropped during preprocessing, missingness summaries, rows used for PCA, rows excluded from PCA, and rows imputed. itemModeThe PCA mode used. itemBackendThe PCA backend used. itemCenterLogical flag indicating whether centering was used. itemScaleLogical flag indicating whether scaling was used.

**See also:** None documented.

## `CreateProjectFolders`

**Purpose:** Create Project Folder Structure

**Canonical usage**
```r
CreateProjectFolders(base_path = ".")
```

**Description:** This function creates a project folder structure with the following directories: itemize item Data/ itemize item Raw/ item Clean/ item Scripts/ item Reports/

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `base_path`: A character string specifying the base directory where the folders should be created. Defaults to the current working directory.

**Returns:** A message indicating that the project structure has been created successfully.

**See also:** None documented.

## `CreateRCIObject`

**Purpose:** Create a Reliable Change Index (RCI) object

**Canonical usage**
```r
CreateRCIObject(
  data,
  variables,
  DataFormat = c("wide", "long"),
  id_var,
  Method = "regression",
  BaselineSpecifier = NULL,
  FollowupSpecifier = NULL,
  SpecifierPosition = c("suffix", "prefix"),
  VisitColumn = NULL,
  VisitOrder = NULL,
  BaselineVisit = NULL,
  Confidence = 0.95,
  Relabel = TRUE
)
```

**Description:** Learn regression-based Reliable Change Index (RCI) models relative to a user-defined reference visit and calculate projected RCI values.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of canonical variable names.
- `DataFormat`: Either "wide" or "long".
- `id_var`: ID column.
- `Method`: Currently only "regression" is supported.
- `BaselineSpecifier`: Baseline visit identifier for wide data.
- `FollowupSpecifier`: Follow-up visit identifier for wide data.
- `SpecifierPosition`: Either "suffix" or "prefix".
- `VisitColumn`: Visit column for long data.
- `VisitOrder`: Optional ordering of visits.
- `BaselineVisit`: Reference visit used for RCI calculations.
- `Confidence`: Confidence interval threshold.
- `Relabel`: Logical; use variable labels when available. #' subsectionInterpretation guide egression-based RCI values are interpreted similarly to z-scores.tabularll RCI cutoff tab Approximate confidence interval cr +/-0.50 tab ~38% cr +/-1.00 tab ~68% cr +/-1.645 tab ~90% cr +/-1.96 tab ~95% cr +/-2.58 tab ~99% cr Traditional Jacobson-Truax RCI thresholds typically use +/-1.96, corresponding to approximately 95% confidence.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `ID`: strongDeprecated (since 19.15.0). Use codeid_var instead.

**Returns:** A SciDataReportR_RCI object.

**See also:** None documented.

## `CreateSOMClusterModel`

**Purpose:** SOM + latent profile clustering pipeline (with AHP and distance baselines)

**Canonical usage**
```r
CreateClusterModel_SOM_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize", "explore"),
  k_range = 2:10,
  models = c(1, 2, 3, 6),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  seed_som = 934521L,
  seed_lpa = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_var = NULL,
  lpa_progress = FALSE,
  lpa_em_itmax = 100L,
  lpa_em_tol = 1e-05,
  lpa_timeout_seconds = 120,
  lpa_drop_zero_sd = TRUE,
  lpa_zero_sd_tol = 1e-08,
  skip_model_after_n_failures = 2L,
  slow_fit_seconds = 120,
  min_nodes_per_cluster = 5,
  high_dist_quantile = 0.95,
  low_prob_threshold = 0.7,
  stability_resamples = 0L,
  stability_seed = 934522L,
  stability_progress = FALSE,
  stability_cores = NULL,
  .NodeClusterFn = NULL
)

Pipeline_SOM_MClust(...)

Pipeline_SOMClust(...)

CreateSOMClusterModel(...)
```

**Description:** End-to-end pipeline to: itemize item Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied Z-score object. item Fit a Self-Organizing Map (SOM; kohonen) on complete cases. item Generate aweSOM visualizations (Circular, Line, Cloud) with optional relabeling using variable labels from the original data frame. item Cluster SOM codebook vectors using latent profile analysis (tidyLPA / mclust backend). item In codemethod = "exploratory", fit a grid of models and select a recommended solution using an Analytic Hierarchy Process (AHP)-style index combining AIC, BIC, and Entropy. item In codemethod = "finalize", fit a user-specified model and number of profiles. item Map node-level clusters and posterior probabilities back to individuals. item Store training variable summaries used later to quantify whether projected cohorts fall outside the original training range. This supports a train once, project many clinical phenotyping workflow: learn phenotype structure in a training cohort, then project new cohorts into the fixed phenotype space without reclustering. Ideal use: correlated continuous clinical or biomarker measures where a topology-preserving map is clinically informative before model-based profiles. Missing data: itemize item SOM and clustering are fit only on rows with complete Z-scores. item The returned codeDataWithClusters has exactly the original rows and columns plus one cluster column; rows not used in SOM/LPA get NA. item The returned codeProbFit$individual is also full length, preserving one row per input row with NA posterior probabilities for rows excluded from SOM/LPA. Z-score behavior: itemize item codeZScoreType = "Center and Scale"/"Center Only"/"Scale Only" computes Z-scores from codedf via codeCreateZScoreObject(). item codeZScoreType = "ZScoreObj" projects Z-scores using an external codeZScoreObj via codeProjectZScore(). item codeZScoreType = "PreZScored" uses existing Z-score columns in codedf as-is and does not re-zscore. Readable SOM + Mclust workflow wrapper for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust(). Compatibility wrapper for codelink[=Pipeline_SOM_MClust]Pipeline_SOM_MClust(). Deprecated alias for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateClusterModel_SOM_MClust`, `Pipeline_SOM_MClust`, `Pipeline_SOMClust`

**Arguments**
- `data`: Data frame containing the variables to be used in SOM and clustering.
- `variables`: Optional character vector of variable names. If NULL, numeric variables are auto-detected using codeSciDataReportR::getNumVars(df, Ordinal = FALSE). In codeZScoreType = "PreZScored", this can also be NULL if you supply codeZScoreVars or if Z-score columns can be auto-detected by prefix.
- `method`: One of code"exploratory" (default) or code"finalize". In code"exploratory", a grid of models is fit and AHP chooses the recommended solution. In code"finalize", the user must specify codefinal_k and codefinal_model.
- `k_range`: Integer vector of numbers of clusters/profiles to consider in exploratory mode. Default code2:10.
- `models`: Integer vector of model specifications for tidyLPA's mclust backend. Model 1 uses equal variance and zero covariance; model 2 uses varying variance and zero covariance; model 3 uses equal variance and equal covariance; and model 6 uses varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are cluster-specific. Supported values and the default are codec(1, 2, 3, 6). Models 4 and 5 require OpenMx and are intentionally unsupported.
- `final_k`: Integer; number of profiles for codemethod = "finalize".
- `final_model`: Integer; model specification for codemethod = "finalize" (should be one of codemodels).
- `ClusterVariableName`: Name of the cluster column in the output. Defaults to code"Cluster". If this column already exists in codedf, it is overwritten (with a message).
- `ZScoreType`: One of: itemize item code"Center and Scale" (default) item code"Center Only" item code"Scale Only" item code"ZScoreObj" (use an existing ZScore object) item code"PreZScored" (use existing Z-score columns in df as-is)
- `ZScoreObject`: Optional ZScoreObj (from codeCreateZScoreObject() or codeProjectZScore()) to use when codeZScoreType = "ZScoreObj".
- `som_xdim, som_ydim`: Optional integers for SOM grid dimensions. If NULL, a square grid with side length codeceiling(n_complete^(1/3)) is used.
- `som_topo`: SOM topology for codekohonen::somgrid(), default code"hexagonal".
- `som_neigh`: SOM neighbourhood function, default code"gaussian".
- `seed_som, seed_lpa`: Integer seeds for SOM and LPA steps (defaults 934521 and 93421).
- `Relabel`: Logical; if TRUE (default), aweSOM plots are relabeled using variable labels from the emphoriginal codedf (via Hmisc or sjlabelled when available) by stripping the Z-score prefix.
- `ZScorePrefix`: Character prefix used for Z-score columns when codeZScoreType = "PreZScored". Default code"Z_".
- `ZScoreVars`: Optional character vector of Z-score column names to use when codeZScoreType = "PreZScored". If NULL, the function attempts to infer them from codevariables or by detecting columns starting with codeZScorePrefix.
- `id_var`: Optional character scalar. If provided and present in codedf, this column is carried into codeProbFit$individual for convenience.
- `lpa_progress`: Logical; if TRUE, print short progress messages while fitting model/profile combinations.
- `lpa_em_itmax`: Integer; maximum number of EM iterations passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_em_tol`: Numeric; EM convergence tolerance passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_timeout_seconds`: Optional timeout in seconds for individual LPA fits. Use NULL to disable timeouts.
- `lpa_drop_zero_sd`: Logical; if TRUE, remove SOM code dimensions with near-zero standard deviation before LPA.
- `lpa_zero_sd_tol`: Numeric tolerance used when codelpa_drop_zero_sd = TRUE.
- `skip_model_after_n_failures`: Optional integer; skip a model family after this many failures.
- `slow_fit_seconds`: Optional runtime threshold used to flag slow LPA fits in diagnostics.
- `min_nodes_per_cluster`: Optional minimum average SOM nodes per cluster considered before attempting a candidate profile count.
- `high_dist_quantile`: Numeric value between 0 and 1 used to define high SOM-distance flags from the training distance distribution. Default is code0.95.
- `low_prob_threshold`: Numeric posterior probability threshold used to flag uncertain phenotype membership. Default is code0.70.
- `stability_resamples`: Number of 90% participant subsample refits used to assess reproducibility for every successful exploratory candidate. Subsamples are drawn without replacement and reuse the reference model's resolved SOM grid dimensions. Defaults to code0 (disabled); use code50 for an exploratory stability screen.
- `stability_seed`: Integer seed for participant subsampling.
- `stability_progress`: Logical; if TRUE, print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `id_col`: strongDeprecated (since 19.15.0). Use codeid_var instead.
- `.NodeClusterFn`: Internal. A function taking the SOM codebook matrix and returning a list with a codenode_cluster integer vector (one label per SOM node) and, optionally, codefit_table, codeahp_best_row, coderecommendation, codebest_fit_name, and codefit_plot. When supplied, the SOM codebook is clustered by that function and the latent-profile grid is not fitted. Used by codelink[=CreateClusterModel_SOM_HDBSCAN]CreateClusterModel_SOM_HDBSCAN(); not part of the user-facing API.
- `...`: Arguments passed to codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Returns:** A list of class code"Pipeline_SOM_MClust" with components: itemize item codemethod, codevars_used, codeZScoreType, codeZScoreObject, codeZScoreVars, codeClusterVariableName item codeDataWithClusters: original codedf with only the cluster column appended item codefit_plot: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model (plus reproducibility when subsample stability is enabled) item codeModelInfo_SOM: list with codesom_model, codesom_codes, codesom_grid, codetraining_variable_summary, codeSOMFit (distance diagnostics, baselines, and per-cluster flags), codeplots (aweSOM plots) item codeModelInfo_MClust: list with codelpa_models, codefit_table, codeAHP information, and codediagnostics for LPA warnings, failures, runtimes, and preprocessing item codeModelInfo_MClust$Stability: subsample replicate, cluster recovery, and summary tables when codestability_resamples > 0 item codeProbFit: list with codenode (node-level posterior probabilities), codeindividual (full-length per-person mapping and probabilities), and probability plots

**See also:** None documented.

## `CreateStatisticsTable`

**Purpose:** Create Statistics Table

**Canonical usage**
```r

```

**Description:** Generate a table of statistics including means, standard deviations, counts, and p-values.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame containing the variables of interest.
- `TargetVar`: The target variable for which statistics will be calculated.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A formatted HTML table displaying statistics.

**See also:** None documented.

## `CreateSummaryTable`

**Purpose:** Create Summary Table

**Canonical usage**
```r
CreateSummaryTable(
  data,
  variables = NULL,
  digits = 2,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  ScrollBoxHeight = "700px"
)
```

**Description:** Generate a descriptive summary table for specified variables in a dataset.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset containing the variables of interest.
- `variables`: A character vector specifying the variables for which summary statistics will be calculated.
- `digits`: Number of decimal places to round the summary statistics.
- `Relabel`: Logical, indicating whether to use variable labels as column headers.
- `Ordinal`: Deprecated logical compatibility option; use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. This numeric descriptive table accepts code"Continuous" or code"Exclude".
- `ScrollBoxHeight`: Height of the scroll box for displaying the table.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `numdecimals`: strongDeprecated (since 19.15.0). Use codedigits instead.

**Returns:** A formatted HTML table displaying summary statistics.

**See also:** None documented.

## `CreateVariableTypesTemplate`

**Purpose:** Create a Template for Variable Types

**Canonical usage**
```r
CreateVariableTypesTemplate(
  data,
  CSVFileName = NULL,
  GuessCategorical = TRUE
)
```

**Description:** Generates the variable-types table for a data frame, already in the format codelink[=RevalueData]RevalueData() expects: one row per column, with the variable's name, its label, a guessed type, and empty columns waiting to be filled in. Optionally writes it straight to CSV.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the variables to be summarized.
- `CSVFileName`: A string specifying the path and name of the CSV file to save the summary. If NULL (the default), the CSV file will not be created.
- `GuessCategorical`: A logical variable specifying if the function should guess what variables are categorical based on having <= 5 unique values
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A data frame with the following columns: describe itemVariableThe names of the variables in the input data frame. itemLabelThe labels of the variables, if available; otherwise, the variable names. itemTypeThe data types of the variables, converted to more user-friendly descriptions. itemCategoryA placeholder column for categorizing variables (default is NA). itemRecodeA placeholder column for recoding information (default is NA). itemCodeA placeholder column for code information (default is NA). itemNotesA placeholder column for any additional notes (default is an empty string). itemExcludeA placeholder column for exclusion flags (default is NA).

**See also:** codelink[=RevalueData]RevalueData() to apply an edited template, codelink[=UpdateDataDictionary]UpdateDataDictionary() to add rows for new variables without losing existing edits, and codelink[=FormattedDataDictionary]FormattedDataDictionary() to render the finished codebook.

## `CreateZScoreObject`

**Purpose:** Calculate Z-scores (or standardized scores) and return data + parameters

**Canonical usage**
```r
CreateZScoreObject(
  data,
  variables = NULL,
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)

CalcZScore(...)
```

**Description:** Standardizes each variable to a common scale and, critically, returns the constants used to do it so the identical transformation can be replayed on other data later. codeCalcZScore() has been superseded by codeCreateZScoreObject(). It remains available as a backwards-compatible alias and returns the same reusable Z-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CalcZScore`

**Arguments**
- `data`: Data frame with variables to standardize.
- `variables`: Character vector of variable names. If NULL, uses SciDataReportR::getNumVars(df).
- `names_prefix`: Prefix to prepend to variable names (default "Z_").
- `RetainLabels`: Logical; if TRUE and Hmisc is available, copy labels.
- `RenameLabels`: Logical; if TRUE, apply the same prefix to labels.
- `center`: Logical; if TRUE, subtract the mean.
- `scale`: Logical; if TRUE, divide by the SD.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=CreateZScoreObject]CreateZScoreObject().

**Returns:** An object of class "ZScoreObj", a list with: itemize item ZScores: data frame of standardized variables only item DataWithZ: original df + standardized variables item Parameters: data frame with Variable, N, Mean, SD item Center: logical flag used item Scale: logical flag used

**See also:** codelink[=ProjectZScore]ProjectZScore() to apply stored parameters to new data, and codelink[=CreateNormativeTScoreModel]CreateNormativeTScoreModel() when the reference values should also be adjusted for covariates such as age or education.

## `CreateZScorePlot`

**Purpose:** Plot Z-score group differences with statistical significance

**Canonical usage**
```r
PlotZScore(
  data,
  TargetVar,
  variables,
  VariableCategories = NULL,
  Relabel = TRUE,
  sort = TRUE,
  RemoveXAxisLabels = TRUE,
  TreatOrdinalAs = "Continuous",
  Parametric = TRUE,
  SigP_YCoord = 1.5,
  SigFDR_YCoord = 1.6
)

CreateZScorePlot(...)
```

**Description:** This function generates a Z-score plot to compare multiple variables across different groups. It offers options for parametric or non-parametric tests, ordinal treatment, and custom labeling. Significant p-values and FDR-adjusted p-values are highlighted on the plot. codeCreateZScorePlot() has been superseded by codePlotZScore(). It remains available as a backwards-compatible alias and returns the same scientific visualization.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `PlotZScore`

**Arguments**
- `data`: A dataframe containing the data to be analyzed.
- `TargetVar`: A string specifying the column name of the grouping variable.
- `variables`: A vector of strings specifying the column names of the variables to be analyzed.
- `VariableCategories`: An optional vector categorizing the variables.
- `Relabel`: Logical; if TRUE, variables will be relabeled using their labels from the dataframe.
- `sort`: Logical; if TRUE, variables will be sorted by category and p-value.
- `RemoveXAxisLabels`: Logical; if TRUE, X-axis labels will be removed.
- `Ordinal`: strongDeprecated (since 20.20.0). Use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. This numeric plot accepts code"Continuous" or code"Exclude".
- `Parametric`: Logical; if TRUE, parametric tests (t-test/ANOVA) will be used; otherwise, non-parametric tests (Wilcoxon/Kruskal-Wallis) will be used.
- `SigP_YCoord`: Numeric; the y-coordinate for marking significant p-values.
- `SigFDR_YCoord`: Numeric; the y-coordinate for marking significant FDR-adjusted p-values.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `...`: Arguments passed to codelink[=PlotZScore]PlotZScore().

**Returns:** A ggplot object representing the Z-score plot.

**See also:** None documented.

## `DeriveFreesurferVolumes`

**Purpose:** Derive Freesurfer bilateral measures and optional ICV-adjusted ratios

**Canonical usage**
```r
DeriveFreesurferVolumes(
  data,
  icv_var = NULL,
  derive_icv_ratios = TRUE,
  bilateral_method = c("sum", "mean"),
  verbose = TRUE
)
```

**Description:** Automatically derives bilateral Freesurfer measures from ASEG and DKT outputs, using either sums or means, and can create intracranial-volume-adjusted ratios.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing Freesurfer ASEG and/or DKT variables.
- `icv_var`: Optional single character string naming the intracranial volume column to use for ratios. When codeNULL and ratios are requested, codeEstimatedTotalIntraCranialVol or codeeTIV is detected automatically.
- `derive_icv_ratios`: Logical. If codeTRUE, derive ICV-adjusted ratios. Default is codeTRUE.
- `bilateral_method`: Character string specifying whether matched left/right measures are combined using a code"sum" or code"mean". Default is code"sum". Output names retain the verb_total suffix for compatibility.
- `verbose`: Logical. If codeTRUE, prints a short summary of derived variables. Default is codeTRUE.

**Returns:** A data frame containing only newly derived variables, with the same number of rows as codedata. A derivation log is stored in the attribute code"Freesurfer_derivation_log".

**See also:** None documented.

## `EvaluateBiomarkerPerformance`

**Purpose:** Evaluate biomarker performance

**Canonical usage**
```r
EvaluateBiomarkerPerformance(
  data,
  outcome_var,
  biomarker_var,
  covariates = NULL,
  PositiveLevel = NULL,
  OutcomeType = c("auto", "binary", "continuous"),
  ThresholdMethod = c("youden", "sensitivity", "specificity", "custom"),
  ThresholdValue = NULL,
  RawThresholdValue = NULL,
  ProbabilityThresholdValue = NULL,
  Validation = c("none", "bootstrap", "cross_validation"),
  BootstrapR = 500,
  CVFolds = 10,
  CIBootstrapR = 500,
  CILevel = 0.95,
  CalibrationGroups = 10,
  Seed = 123,
  Relabel = TRUE,
  codebook = NULL,
  Verbose = TRUE
)
```

**Description:** Evaluates a continuous or categorical biomarker against a binary or continuous outcome. Binary analyses use ordinary logistic regression. Models with separation, non-convergence, or aliased coefficients return stable unavailable metrics.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `outcome_var`: One outcome variable name.
- `biomarker_var`: One biomarker variable name.
- `covariates`: Optional covariate variable names.
- `PositiveLevel`: Positive binary outcome level, or codeNULL to use the second observed level.
- `OutcomeType`: One of code"auto", code"binary", or code"continuous".
- `ThresholdMethod`: One of code"youden", code"sensitivity", code"specificity", or code"custom".
- `ThresholdValue`: Target sensitivity or specificity.
- `RawThresholdValue`: Custom raw-biomarker threshold.
- `ProbabilityThresholdValue`: Custom predicted-probability threshold.
- `Validation`: One of code"none", code"bootstrap", or code"cross_validation".
- `BootstrapR`: Number of bootstrap optimism-correction resamples.
- `CVFolds`: Number of cross-validation folds.
- `CIBootstrapR`: Number of bootstrap confidence-interval resamples.
- `CILevel`: Confidence level.
- `CalibrationGroups`: Maximum grouped-calibration bins.
- `Seed`: Random seed.
- `Relabel`: Use codebook labels, then label attributes, for presentation.
- `codebook`: Optional data frame with codeVariable and codeLabel.
- `Verbose`: Print positive-level information.

**Returns:** A stable named list containing models, performance, thresholds, predictions, calibration, validation, plots, and metadata.

**See also:** None documented.

## `ExploreDatasetComparison`

**Purpose:** Explore dataset comparison results interactively

**Canonical usage**
```r
ExploreDatasetComparison(
  CompareObj,
  Title = "Dataset comparison explorer",
  TopN = 10
)
```

**Description:** Create an interactive HTML dashboard from a codeCompareDatasets() result object. This function is designed for data review and quality-control workflows. It displays high-level summary cards, a traffic-light checks table, and an expandable variable-change explorer that shows side-by-side old and new values for modified cells.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `CompareObj`: A list returned by codeCompareDatasets().
- `Title`: Character title shown at the top of the dashboard. Default is code"Dataset comparison explorer".
- `TopN`: Integer number of example variables or records to show in previews and expanded sections. Default is code10.

**Returns:** An codehtmltools::tagList() object containing an interactive dashboard.

**See also:** None documented.

## `ExploreMergeValidation`

**Purpose:** Explore merge validation results interactively

**Canonical usage**
```r
ExploreMergeValidation(
  MergeObj,
  Title = "Merge validation explorer",
  TopN = 10,
  TableHeight = 350,
  Detail = c("Compact", "Full")
)
```

**Description:** Create an interactive HTML dashboard from a codeValidateMerge() result object. This function is designed for merge quality-control workflows. It displays a traffic-light checks table (with rows/columns/unique-key context folded in as informational rows), a coverage explorer, an expandable duplicate-variable conflict explorer, and suggested actions.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `MergeObj`: A list returned by codeValidateMerge().
- `Title`: Character title shown at the top of the dashboard. Default is code"Merge validation explorer".
- `TopN`: Integer number of example variables or records to show in previews and expanded sections. Default is code10.
- `TableHeight`: Height in pixels for scrollable reactable tables. Default is code350.
- `Detail`: Either code"Compact" (default) or code"Full". In code"Compact" mode, the coverage explorer and conflicts explorer render as collapsed click-to-expand accordion sections labeled with their item counts. In code"Full" mode, the same sections are expanded by default. In both modes, sections with nothing to show (no unmatched keys, no duplicated variables, no suggested actions) are omitted entirely.

**Returns:** An codehtmltools::tagList() object containing an interactive dashboard.

**See also:** None documented.

## `ExtractPCAComponentSummary`

**Purpose:** Extract PCA component summaries

**Canonical usage**
```r
ExtractPCAComponentSummary(
  PCAObject,
  loading_threshold = 0.4,
  top_n = NULL,
  use_labels = TRUE,
  html_format = TRUE
)
```

**Description:** Extract variables contributing to each PCA component based on an absolute loading threshold. Returns both a tidy long-format table and compact summary tables suitable for reporting.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `PCAObject`: Output object from CreatePCAObject().
- `loading_threshold`: Minimum absolute loading required for inclusion. Default is 0.4.
- `top_n`: Optional maximum number of contributors per component. If NULL, all contributors above threshold are retained.
- `use_labels`: Logical indicating whether variable labels should be used when available. Default TRUE.
- `html_format`: Logical indicating whether negative contributors should be formatted using red HTML text. Default TRUE.

**Returns:** A list containing: itemLongTable A tidy tibble with one row per contributor. itemSummaryTable A compact tibble with one row per component and comma-separated contributor summaries. itemSummaryTableLines A compact tibble with one row per component and line-separated contributor summaries. itemFormattedSummaryTable A formatted gt table with comma-separated contributors. itemFormattedSummaryTableLines A formatted gt table with line-separated contributors.

**See also:** codelink[=CreatePCAObject]CreatePCAObject() to fit the PCA, codelink[=CreatePCATable]CreatePCATable() for the variance-explained table, and codelink[=ProjectPCA]ProjectPCA() to score new data on the components once they have been interpreted.

## `FormattedDataDictionary`

**Purpose:** Create a formatted data dictionary table

**Canonical usage**
```r
FormattedDataDictionary(
  data,
  digits = 2
)
```

**Description:** This function generates a formatted data dictionary table using the specified data frame. The table includes variable names, labels, types, and additional formatting based on variable types.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame for which the data dictionary is to be created.
- `digits`: Number of decimals to display for numeric variables (default: 2).
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `numdecimals`: strongDeprecated (since 19.15.0). Use codedigits instead.

**Returns:** A formatted data dictionary table (gt object).

**See also:** None documented.

## `FreezeTableHeader`

**Purpose:** Freeze the header row of a long table when scrolling

**Canonical usage**
```r
FreezeTableHeader(
  x,
  height = NULL,
  width = NULL,
  header_background = "white",
  bootstrap_options = c("striped", "hover", "condensed"),
  full_width = FALSE,
  font_size = NULL,
  ...
)
```

**Description:** Render a table with a "sticky" header row that stays visible while the reader scrolls, so the column meanings are never lost in a long table. This is handy for the tall tables produced by codelink[=MakeComparisonTable]MakeComparisonTable() and other codegtsummary helpers when they are knit into HTML Quarto or R Markdown documents.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `x`: A codegtsummary object (such as the output of codelink[=MakeComparisonTable]MakeComparisonTable()), a data frame or tibble, or an existing pkgkableExtra / codeknitr_kable object.
- `height`: Optional CSS height for a scroll box, for example code"400px" or code"60vh". When supplied, the table scrolls within a box of this height with the header frozen. When codeNULL (default), the header sticks during page-level scrolling instead.
- `width`: Optional CSS width for the scroll box, for example code"100%". Only used when codeheight is supplied.
- `header_background`: Background color for the frozen header row. Defaults to code"white" so the header cleanly covers rows scrolling underneath it.
- `bootstrap_options`: Character vector of pkgkableExtra bootstrap styling options passed to codekableExtra::kable_styling(). Defaults to codec("striped", "hover", "condensed").
- `full_width`: Logical; passed to codekableExtra::kable_styling(). Default codeFALSE.
- `font_size`: Optional numeric font size passed to codekableExtra::kable_styling().
- `...`: Additional arguments passed to codekableExtra::kable_styling().

**Returns:** A pkgkableExtra HTML table object with a frozen header, suitable for printing in a Quarto or R Markdown chunk.

**See also:** codelink[=MakeComparisonTable]MakeComparisonTable()

## `geom_starcaption`

**Purpose:** Add a Caption Explaining Star Annotations

**Canonical usage**
```r
geom_starcaption()
```

**Description:** This function adds a caption to a ggplot explaining the meaning of star annotations (*, **, ***). It is most commonly added to correlation heatmaps produced by codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap() or to downstream plots derived from that function with codelink[=add_r_and_stars]add_r_and_stars().

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- None documented.

**Returns:** A codelabs() object that can be added to a ggplot, especially SciDataReportR heatmaps that use star annotations.

**See also:** None documented.

## `getBinaryVars`

**Purpose:** Identify Binary Variables

**Canonical usage**
```r
getBinaryVars(
  data,
  Ordinal = TRUE,
  Revalued = TRUE
)
```

**Description:** This function identifies and returns a list of binary variables in a dataframe. Binary variables are defined as having exactly two unique values or levels. The function supports options for handling ordinal factors and revalued data.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe to analyze for binary variables.
- `Ordinal`: Logical. If TRUE, ordinal factors are included in the search for binary variables. Default is TRUE.
- `Revalued`: Logical. If TRUE, the function checks factors and their levels; otherwise, it checks for variables with two unique values.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A character vector containing the names of binary variables.

**See also:** codelink[=createBinaryMapping]createBinaryMapping() to fix which level counts as positive, and codelink[=getCatVars]getCatVars() / codelink[=getNumVars]getNumVars() for the other partitions.

## `getCatVars`

**Purpose:** Get Categorical Variables

**Canonical usage**
```r

```

**Description:** Extracts categorical variables from a data frame.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame from which to extract categorical variables.
- `Ordinal`: Logical, indicating whether to include ordinal variables.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A character vector containing the names of categorical variables.

**See also:** codelink[=getNumVars]getNumVars() and codelink[=getBinaryVars]getBinaryVars() for the other partitions.

## `getNumVars`

**Purpose:** Get Numeric Variables

**Canonical usage**
```r

```

**Description:** Extracts numeric variables from a data frame.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame from which to extract numeric variables.
- `Ordinal`: Logical, indicating whether to include ordinal variables.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A character vector containing the names of numeric variables.

**See also:** codelink[=getCatVars]getCatVars(), codelink[=getBinaryVars]getBinaryVars(), and codelink[=ConvertOrdinalToNumeric]ConvertOrdinalToNumeric() for the ordinal policy these share.

## `InsertValues`

**Purpose:** Insert Values into a String Array

**Canonical usage**
```r
InsertValues(vec, values, after_value, location = "after")
```

**Description:** This function inserts one or more strings into an array of strings at specified locations.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `vec`: A character vector in which values will be inserted.
- `values`: A character vector or a single string to insert into the original vector.
- `after_value`: A character string indicating the element in codevec before or after which the values will be inserted. Should have a length of one
- `location`: A character string specifying whether to insert codevalues "before" or "after" the codeafter_value. Defaults to "after".

**Returns:** A character vector with the values inserted.

**See also:** None documented.

## `InspectCategoricalSummary`

**Purpose:** Inspect categorical variables

**Canonical usage**
```r
InspectCategoricalSummary(
  data,
  variables = NULL,
  codebook = NULL,
  IncludeMissing = TRUE,
  MissingLabel = "(Missing)",
  RetainLabels = TRUE,
  SortLevelsBy = c("Frequency", "Value", "None"),
  SortVariablesBy = c("Input", "Label", "MissingPercent", "UniqueLevels", "TotalN"),
  Descending = TRUE,
  MaxLevels = 30,
  Plot = TRUE,
  PlotType = c("bar", "lollipop"),
  FacetScales = c("free_y", "fixed"),
  UsePercent = TRUE,
  LabelBars = TRUE,
  WrapLabels = 35,
  BaseSize = 11
)
```

**Description:** Provides a SciDataReportR-native categorical inspection summary and plot inspired by codeinspectdf::inspect_cat(), created because inspectdf is archived and no longer available on CRAN.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `variables`: Optional character vector of categorical variables to summarize. If codeNULL, character, factor, logical, labelled, and haven-labelled columns are detected automatically.
- `codebook`: Optional data frame with codeVariable and codeLabel columns.
- `IncludeMissing`: Logical. If codeTRUE, missing values are included as a level and in percentages.
- `MissingLabel`: Character label to display for missing values.
- `RetainLabels`: Logical. If codeTRUE, use variable labels from codeCodebook and value labels from labelled variables when available.
- `SortLevelsBy`: One of code"Frequency", code"Value", or code"None".
- `SortVariablesBy`: One of code"Input", code"Label", code"MissingPercent", code"UniqueLevels", or code"TotalN".
- `Descending`: Logical. If codeTRUE, sort selected summaries in descending order.
- `MaxLevels`: Positive integer giving the maximum number of levels to display per variable in the plot before collapsing lower-ranked levels to code"(Other)".
- `Plot`: Logical. If codeTRUE, return a ggplot object.
- `PlotType`: One of code"bar" or code"lollipop".
- `FacetScales`: One of code"free_y" or code"fixed".
- `UsePercent`: Logical. If codeTRUE, plot percentages; otherwise plot counts.
- `LabelBars`: Logical. If codeTRUE, add readable labels to bars or points.
- `WrapLabels`: Integer number of characters used to wrap displayed level and facet labels.
- `BaseSize`: Base font size passed to codetheme_minimal().
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `Codebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A named list with codeSummary, a tibble containing categorical counts and percentages, and codePlot, a ggplot object when codePlot = TRUE or codeNULL when codePlot = FALSE.

**See also:** None documented.

## `InspectFile`

**Purpose:** Inspect a scientific data file before import

**Canonical usage**
```r
InspectFile(
  path,
  sheet = NULL,
  preview_rows = 20,
  check_styles = TRUE,
  check_sheets = TRUE,
  check_header = TRUE,
  quiet = FALSE
)
```

**Description:** codeInspectFile() checks common import risks before a file is read into an analysis workflow. It is especially useful for Excel workbooks, where multiple sheets, metadata rows, unnamed columns, duplicate column names, and workbook formatting can change how the data should be interpreted.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `path`: Path to the file.
- `sheet`: Sheet name or index for Excel files. If codeNULL, the first sheet is inspected.
- `preview_rows`: Number of rows to preview when checking header and column name issues.
- `check_styles`: Logical. For code.xlsx files, check whether workbook styles or formatting exist. This can be slower for large workbooks.
- `check_sheets`: Logical. For Excel files, check whether the workbook has multiple sheets.
- `check_header`: Logical. For Excel files, attempt to detect whether the header row is not row 1.
- `quiet`: Logical. If codeFALSE, print a compact inspection summary.

**Returns:** Invisibly returns a list containing file inspection metadata.

**See also:** None documented.

## `IQROutliers`

**Purpose:** Detect outliers using the Tukey IQR rule and visualize results

**Canonical usage**
```r
IQROutliers(
  data,
  Variable,
  id_var = NULL,
  group = NULL
)
```

**Description:** This function identifies potential outliers in a numeric variable using the Tukey interquartile range (IQR) rule (Tukey, 1977). It returns a tibble of the detected outlier rows and a ggplot visualization showing the variable across groups with outlier points highlighted.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame or tibble containing the variable to evaluate.
- `Variable`: A string specifying the name of the numeric variable to test.
- `id_var`: A string specifying the identifier column to include in the returned outlier table. If codeNULL, no ID column is included in the returned table. Defaults to codeNULL.
- `group`: A string specifying the grouping or batch column to use on the x-axis of the diagnostic plot. If codeNULL, the function will produce a single combined boxplot across all rows. Defaults to codeNULL.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `id`: strongDeprecated (since 19.15.0). Use codeid_var instead.

**Returns:** A list with two elements: itemize item codeoutlierdf: a tibble containing only rows flagged as outliers, including the ID (when requested), variable, group (when requested), and outlier flag. item codep: a ggplot2 object showing a boxplot and jittered points colored by outlier status.

**See also:** None documented.

## `KeepEnv`

**Purpose:** Keep selected objects in an environment and remove everything else

**Canonical usage**
```r
KeepEnv(
  Keep,
  Env = parent.frame(),
  DryRun = FALSE,
  Invert = FALSE,
  Quiet = FALSE
)
```

**Description:** Keeps only the objects specified in codeKeep within codeEnv and removes all other objects in that environment. Optionally returns a summary of what was removed.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `Keep`: Character vector of object names to keep.
- `Env`: Environment to clean. Defaults to the calling environment.
- `DryRun`: If TRUE, does not remove anything and only reports what would be removed.
- `Invert`: If TRUE, removes only codeKeep and keeps everything else.
- `Quiet`: If TRUE, suppresses messages.

**Returns:** Invisible list with codekept and coderemoved vectors.

**See also:** None documented.

## `Make_DataDictionary`

**Purpose:** Create a data dictionary for a data frame

**Canonical usage**
```r
MakeDataDictionary(
  data,
  digits = 2
)

Make_DataDictionary(...)
```

**Description:** This function generates a stable data dictionary for a data frame using codecodebook::skim_codebook(). It is designed to work even when skim output omits type-specific summary columns, such as numeric summaries for data frames without numeric variables. codeMake_DataDictionary() has been superseded by codeMakeDataDictionary(). It remains available as a backwards-compatible alias and returns the same data dictionary.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `MakeDataDictionary`

**Arguments**
- `data`: A data frame.
- `digits`: Number of decimals to display for numeric variables.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `numdecimals`: strongDeprecated (since 19.15.0). Use codedigits instead.
- `...`: Arguments passed to codelink[=MakeDataDictionary]MakeDataDictionary().

**Returns:** A data frame with one row per variable and stable summary columns.

**See also:** None documented.

## `MakeComparisonTable`

**Purpose:** Make comparison table with covariate adjustment, effect sizes, and pairwise contrasts

**Canonical usage**
```r
MakeComparisonTable(
  data,
  group_var = NULL,
  variables,
  ...,
  covariates = NULL,
  value_digits = 2,
  p_digits = 3,
  AddEffectSize = FALSE,
  effect_size_digits = 2,
  AddPairwise = FALSE,
  PairwiseMethod = "bonferroni",
  Parametric = TRUE,
  ParametricDisplay = NULL,
  IncludeOverallN = FALSE,
  IncludeMissing = FALSE,
  suppress_warnings = FALSE,
  Referent = NULL,
  IncludeOverallStats = FALSE,
  ShowPositiveBinaryOnLabel = TRUE,
  CatMethod = c("auto", "chisq", "fisher"),
  MultiCatAdjusted = c("multinomial_LR", "none"),
  ShowNotes = c("auto", "always", "never"),
  NotesPosition = c("last", "after_test", "before_pairwise"),
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical"
)
```

**Description:** Create a label-aware comparison table using codegtsummary::tbl_summary() with optional global hypothesis tests, covariate-adjusted tests, effect sizes, and pairwise comparisons.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `group_var`: Character scalar naming the grouping variable.
- `variables`: Character vector of variables to summarize.
- `...`: Optional additional variable names supplied individually.
- `covariates`: Optional character vector of covariates for adjusted models.
- `value_digits`: Number of digits for descriptive statistics.
- `p_digits`: Number of digits for p-values.
- `AddEffectSize`: Logical; add effect-size columns.
- `effect_size_digits`: Number of digits for effect sizes.
- `AddPairwise`: Logical; add pairwise comparison columns.
- `PairwiseMethod`: P-value adjustment method. Use code"none" for no adjustment.
- `Parametric`: Logical; use parametric tests for continuous outcomes.
- `ParametricDisplay`: Logical; display continuous summaries as mean (SD). If codeFALSE, display median linkIQR. Defaults to codeParametric.
- `IncludeOverallN`: Logical; add N column.
- `IncludeMissing`: Logical; include missing rows in summaries.
- `suppress_warnings`: Logical; suppress selected gtsummary warnings.
- `Referent`: Optional reference group for pairwise comparisons.
- `IncludeOverallStats`: Logical; add overall summary column.
- `ShowPositiveBinaryOnLabel`: Logical; for binary variables, show only the positive level where identifiable.
- `CatMethod`: Categorical test method. One of code"auto", code"chisq", code"fisher".
- `MultiCatAdjusted`: Adjusted multicategory method. Currently code"multinomial_LR" or code"none".
- `ShowNotes`: Whether to show the Analysis notes column. One of code"auto", code"always", or code"never".
- `NotesPosition`: Analysis notes column position. One of code"last", code"after_test", or code"before_pairwise".
- `Relabel`: Logical; if TRUE (default), use attached variable labels.
- `TreatOrdinalAs`: How ordinal variables are treated: code"Categorical", code"Continuous", code"Both", or code"Exclude".
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `CompVariable`: strongDeprecated (since 19.15.0). Use codegroup_var instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.
- `ValueDigits`: strongDeprecated (since 19.15.0). Use codevalue_digits instead.
- `pDigits`: strongDeprecated (since 19.15.0). Use codep_digits instead.
- `EffectSizeDigits`: strongDeprecated (since 19.15.0). Use codeeffect_size_digits instead.

**Returns:** A codegtsummary object.

**See also:** None documented.

## `MakeDataDictionary`

**Purpose:** Create a data dictionary for a data frame

**Canonical usage**
```r
MakeDataDictionary(
  data,
  digits = 2
)

Make_DataDictionary(...)
```

**Description:** This function generates a stable data dictionary for a data frame using codecodebook::skim_codebook(). It is designed to work even when skim output omits type-specific summary columns, such as numeric summaries for data frames without numeric variables. codeMake_DataDictionary() has been superseded by codeMakeDataDictionary(). It remains available as a backwards-compatible alias and returns the same data dictionary.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `Make_DataDictionary`

**Arguments**
- `data`: A data frame.
- `digits`: Number of decimals to display for numeric variables.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `numdecimals`: strongDeprecated (since 19.15.0). Use codedigits instead.
- `...`: Arguments passed to codelink[=MakeDataDictionary]MakeDataDictionary().

**Returns:** A data frame with one row per variable and stable summary columns.

**See also:** None documented.

## `MakeFacetCatComparisonTable`

**Purpose:** Create a merged gtsummary table by faceting comparisons across multiple categorical variables

**Canonical usage**
```r
MakeFacetCatComparisonTable(
  data,
  FacetVariables,
  variables,
  covariates = NULL,
  value_digits = 2,
  p_digits = 3,
  AddEffectSize = FALSE,
  effect_size_digits = 2,
  AddPairwise = FALSE,
  PairwiseMethod = "bonferroni",
  Parametric = TRUE,
  ParametricDisplay = NULL,
  IncludeOverallN = FALSE,
  IncludeMissing = FALSE,
  suppress_warnings = FALSE,
  Referent = NULL,
  IncludeOverallStats = FALSE,
  ShowPositiveBinaryOnLabel = TRUE,
  CompFun = MakeComparisonTable,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  ...
)
```

**Description:** Generates a series of comparison tables using codeMakeComparisonTable() for each categorical variable (facet) in the provided list and merges them side-by-side using codegtsummary::tbl_merge(). This function extends the functionality of codeMakeComparisonTable() by automatically detecting which facet variables are categorical (factor or character) and producing a faceted summary of how the main comparison variable (e.g., Cluster, TreatmentArm) differs across multiple categorical dimensions such as Race, Sex, or HIV status.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing all variables to be analyzed.
- `FacetVariables`: A character vector of variable names to facet by. The function automatically selects those that are categorical (codefactor or codecharacter).
- `variables`: A character string naming the variable(s) being compared (e.g., "Cluster").
- `covariates`: Optional character vector of covariate names to adjust for.
- `value_digits`: Number of decimal digits to display for numeric values (default = 2).
- `p_digits`: Number of decimal digits to display for p-values (default = 3).
- `AddEffectSize`: Logical; if TRUE, include effect sizes (default = FALSE).
- `effect_size_digits`: Decimal digits for effect size values (default = 2).
- `AddPairwise`: Logical; if TRUE, include pairwise comparisons (default = FALSE).
- `PairwiseMethod`: Method for pairwise comparison p-value adjustment (default = "bonferroni").
- `Parametric`: Logical; if TRUE, use parametric tests (default = TRUE).
- `ParametricDisplay`: Optional vector specifying which statistics to display for parametric tests.
- `IncludeOverallN`: Logical; if TRUE, adds overall N to the table (default = FALSE).
- `IncludeMissing`: Logical; if TRUE, includes missing categories (default = FALSE).
- `suppress_warnings`: Logical; suppress internal warnings (default = FALSE).
- `Referent`: Optional string specifying the referent category for binary or categorical comparisons.
- `IncludeOverallStats`: Logical; if TRUE, adds overall descriptive statistics (default = FALSE).
- `ShowPositiveBinaryOnLabel`: Logical; if TRUE, labels binary variables with positive outcome (default = TRUE).
- `CompFun`: Comparison function to apply; defaults to codeMakeComparisonTable.
- `Relabel`: Logical; if TRUE (default), use attached variable labels.
- `TreatOrdinalAs`: How ordinal variables are treated in each table.
- `...`: Additional arguments passed to the comparison function.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.
- `ValueDigits`: strongDeprecated (since 19.15.0). Use codevalue_digits instead.
- `pDigits`: strongDeprecated (since 19.15.0). Use codep_digits instead.
- `EffectSizeDigits`: strongDeprecated (since 19.15.0). Use codeeffect_size_digits instead.

**Returns:** A codegtsummary table created by merging each facet's codeMakeComparisonTable() output side-by-side using codegtsummary::tbl_merge(). Each facet variable is labeled with its own tab spanner header for clarity.

**See also:** codelink[=MakeComparisonTable]MakeComparisonTable() for a single grouping, and codelink[=MakeTable1]MakeTable1() for a plain descriptive table with no grouping at all.

## `MakePairwiseHeatmap`

**Purpose:** Make a pairwise referent heatmap

**Canonical usage**
```r
MakePairwiseHeatmap(
  data,
  group_var,
  variables,
  Referent,
  covariates = NULL,
  Parametric = TRUE,
  adjust_scope = c("per_group", "per_variable", "matrix", "none"),
  p_adjust_method = c("fdr", "bonferroni", "holm", "none"),
  star_p = c("raw", "adjusted", "none"),
  adjusted_outline = TRUE,
  adjusted_significance_threshold = 0.05,
  adjusted_outline_color = "black",
  adjusted_outline_linewidth = 1,
  low_color = "#52BCA3FF",
  mid_color = "white",
  high_color = "#E58606FF",
  fill_midpoint = 0,
  fill_limits = NULL,
  fill_oob = scales::squish,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_caption = FALSE,
  x_axis_text_angle = 0,
  return_models = FALSE,
  star_color = "black",
  star_size = 4
)
```

**Description:** Build a heatmap of pairwise group contrasts against a required referent group. Continuous outcomes are always transformed using the referent group before modeling: Z-scores when codeParametric = TRUE, and M-scores when codeParametric = FALSE.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `group_var`: Character scalar naming the grouping variable.
- `variables`: Character vector of continuous outcome variables.
- `Referent`: Character scalar naming the referent level of codegroup_var.
- `covariates`: Optional character vector of covariates.
- `Parametric`: Logical. If codeTRUE, outcomes are Z-scored before modeling. If codeFALSE, outcomes are M-scored and HC3 robust covariance is used for estimated marginal mean contrasts.
- `adjust_scope`: Multiple-comparison correction scope. code"per_group" adjusts across variables within each group-vs-referent contrast; code"per_variable" adjusts across group contrasts within each variable; code"matrix" adjusts across all displayed cells; code"none" applies no correction.
- `p_adjust_method`: Method passed to codelink[stats:p.adjust]stats::p.adjust(). Use code"none" for no correction.
- `star_p`: Which p-values should drive cell stars: raw, adjusted, or none.
- `adjusted_outline`: Logical; outline cells significant after adjustment.
- `adjusted_significance_threshold`: Threshold for adjusted-significant outlines.
- `adjusted_outline_color, adjusted_outline_linewidth`: Appearance of the adjusted-significant outline.
- `low_color, mid_color, high_color`: Diverging heatmap colors.
- `fill_midpoint`: Numeric midpoint for the fill scale.
- `fill_limits`: Optional numeric vector of length 2. If codeNULL, symmetric limits are computed from the observed estimated mean differences.
- `fill_oob`: Out-of-bounds handler for the fill scale.
- `cluster_rows, cluster_columns`: Logical; optionally cluster rows or columns based on estimated mean differences.
- `show_caption`: Logical; add an explanatory caption to the plot.
- `x_axis_text_angle`: Numeric angle for x-axis labels. Defaults to code0.
- `return_models`: Logical; include fitted model objects in the return.
- `star_color, star_size`: Appearance of p-value stars.

**Returns:** An object of class code"SciDataReportRPairwiseHeatmap" with codePlot, codeResults, codeModels, codeSettings, codeScalingParameters, and codeWarnings. codeResults includes readable audit columns such as codeTest, codeContrast, codeAdjustment, and codeModelFormula.

**See also:** None documented.

## `MakeTable1`

**Purpose:** Create Summary Table using gtsummary

**Canonical usage**
```r
MakeTable1(
  data,
  variables = NULL,
  TreatOrdinalAs = "Categorical",
  Relabel = TRUE,
  AutoDetectDistribution = FALSE,
  IncludeMissing = "ifany"
)
```

**Description:** This function is a wrapper around codegtsummary::tbl_summary that ensures continuous variables are treated as continuous.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataframe to create the summary table from.
- `variables`: Optional. A character vector specifying the names of variables to include in the summary table. If NULL, all variables are included.
- `TreatOrdinalAs`: Character. Specifies how ordinal variables should be treated. Can be "Continuous", "Categorical", or "Both".
- `Relabel`: Logical; if TRUE (default), use attached variable labels.
- `AutoDetectDistribution`: Logical. If TRUE, the function will attempt to automatically detect the distribution of variables. Default is FALSE.
- `IncludeMissing`: Character matching gtsummary criteria. Can be "no", "ifany", or "always". Default is "ifany"
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A summary table created using gtsummary.

**See also:** None documented.

## `MakeUnivariateRegressionTable`

**Purpose:** Univariate Regression Table

**Canonical usage**
```r
MakeUnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical"
)

UnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE
)
```

**Description:** Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified). codeUnivariateRegressionTable() was renamed to codeMakeUnivariateRegressionTable() in SciDataReportR 20.5.0 to match the package's verbMake* naming convention. It remains available as a backwards-compatible synonym.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `UnivariateRegressionTable`

**Arguments**
- `data`: Dataframe containing the variables
- `outcome_vars`: Character vector of outcome variable names
- `predictor_vars`: Character vector of predictor variable names
- `covariates`: Character vector of covariate variable names (default: NULL)
- `Standardize`: Logical indicating whether to standardize numeric variables (default: FALSE)
- `Method`: Character. Regression method to use. code"auto" detects linear regression for numeric outcomes and logistic regression for two-level outcomes. code"lm" and code"logistic" force one model family for all outcomes.
- `LogisticExponentiate`: Logical. If codeTRUE, logistic regression estimates are exponentiated and reported as odds ratios.
- `ReturnModels`: Logical. If codeTRUE, return fitted model objects in codeModelSummaries. Default is codeFALSE to keep large screening runs lighter.
- `Relabel`: Logical; if TRUE (default), display attached variable labels.
- `TreatOrdinalAs`: How ordinal outcomes and predictors are handled.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `OutcomeVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `PredictorVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `Covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing: itemize item FormattedTable: A codegt table with formatted regression results item LargeTable: A codegt table with unformatted regression results item Results: A tidy dataframe with one row per estimated term. Columns: codeOutcome, codeOutcomeLabel, codeOutcomeFamily, codeEffectType, codePredictor, codePredictorLabel, codeTerm, codeLevel, codeTermLabel, codeN, codeEstimate, codeStdError, codeConfLow, codeConfHigh, codePValue, codeSignificant, and codeReferenceValue. This dataframe can be filtered and passed directly to codelink[=PlotForestFromTable]PlotForestFromTable(). item ModelSummaries: A list of fitted model objects when codeReturnModels = TRUE, otherwise codeNULL item Metadata: Outcome families and analysis settings

**See also:** codelink[=PlotForestFromTable]PlotForestFromTable() to visualize codeResults, codelink[=MultivariableRegressionTable]MultivariableRegressionTable() for mutually adjusted models, and codelink[=ApplyFDRCorrection]ApplyFDRCorrection() for multiple-comparison correction.

## `Merge_ByClosestTime`

**Purpose:** Merge Two Data Frames by Closest Time

**Canonical usage**
```r
Merge_ByClosestTime(
  DataFrame1,
  DataFrame2,
  TimeVar1,
  TimeVar2,
  keys = NULL,
  is_date = FALSE
)
```

**Description:** This function merges two data frames based on the closest time in the specified time columns. It optionally merges using one or more matching variables (e.g., IDs). The resulting merged data frame contains the closest time matches and time differences.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `DataFrame1`: A data frame containing the first set of data.
- `DataFrame2`: A data frame containing the second set of data.
- `TimeVar1`: The name of the time variable in DataFrame1 (as a string).
- `TimeVar2`: The name of the time variable in DataFrame2 (as a string).
- `keys`: Optional. Character vector of variable(s) to merge by. Must exist in BOTH data frames and be in the same order.
- `is_date`: Logical. Indicates whether the time variables are dates (TRUE) or POSIXct (FALSE).
- `MergeBy`: strongDeprecated (since 19.15.0). Use codekeys instead.

**Returns:** A list with: itemmerged_dataframeData frame with closest time matches itemtime_differencesVector of time differences

**See also:** None documented.

## `merge_detail`

**Purpose:** Print a plain-text detail report for one safe_merge result

**Canonical usage**
```r
merge_detail(m, TopN = 10)
```

**Description:** Print static codeknitr::kable() tables (no plots, no htmlwidgets) for the key diagnostic sections of a codelink[=safe_merge]safe_merge() result: the validation checks, unmatched key combinations on each side, overlapping non-key variables, and suspicious duplicate-variable conflicts. Sections with nothing to show are skipped entirely.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `m`: A list returned by codelink[=safe_merge]safe_merge() (must contain a codevalidation element from codelink[=ValidateMerge]ValidateMerge() and a codelog tibble).
- `TopN`: Integer. Maximum number of rows shown for the left-only and right-only unmatched-key previews. Default is code10.

**Returns:** codem, invisibly. Called for its printed output.

**See also:** codelink[=safe_merge]safe_merge(), codelink[=merge_summary_table]merge_summary_table(), codelink[=ExploreMergeValidation]ExploreMergeValidation()

## `merge_summary_table`

**Purpose:** Combine safe_merge logs into a single summary table

**Canonical usage**
```r
merge_summary_table(merge_log, flagged_only = FALSE)
```

**Description:** Bind the one-row verb$log tibbles produced by codelink[=safe_merge]safe_merge() into a single table, one row per merge, optionally filtered to merges that did not pass cleanly. Useful as an end-of-pipeline rollup after a sequence of merges.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `merge_log`: A list of verb$log tibbles from codelink[=safe_merge]safe_merge() results. Full codelink[=safe_merge]safe_merge() result objects are also accepted (the verb$log element is extracted automatically), as is a single log tibble.
- `flagged_only`: Logical. If codeTRUE, only rows with codeStatus != "PASS" are returned. Default is codeFALSE.

**Returns:** A tibble with one row per merge and the same columns as a codelink[=safe_merge]safe_merge() verb$log tibble.

**See also:** codelink[=safe_merge]safe_merge()

## `MergeCodebooks`

**Purpose:** Merge multiple codebooks using harmonization rules

**Canonical usage**
```r
MergeCodebooks(
  codebooks,
  Rules = NULL,
  VariableCol = "Variable",
  warn = TRUE,
  strict = FALSE
)
```

**Description:** Deterministically merge multiple codebooks using optional harmonization rules generated from codeCodebookMergeApp().

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `codebooks`: Named list of codebook data frames.
- `Rules`: Optional harmonization rules generated from codeCodebookMergeApp().
- `VariableCol`: Name of variable identifier column.
- `warn`: Logical; emit warnings.
- `strict`: Logical; stop on unresolved conflicts.

**Returns:** A list containing: describe itemCodebookMerged harmonized codebook itemConflictReportDetected conflicts itemAppliedRulesApplied rules

**See also:** None documented.

## `MergeFragmentedRecords`

**Purpose:** Merge fragmented records into a single observation

**Canonical usage**
```r
MergeFragmentedRecords(
  data,
  id_var = "subject",
  date_var = "date",
  session_var = "session",
  keep_session = TRUE,
  session_name = "first_session",
  n_rows_name = "n_rows_collapsed",
  arrange_desc_session = FALSE,
  empty_strings_to_na = TRUE
)
```

**Description:** Merges multiple rows representing fragments of the same observation (e.g., participant visit, assessment session, study encounter, or questionnaire administration) into a single row by selecting the first non-missing value within each variable.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing fragmented records.
- `id_var`: Character string specifying the participant identifier variable. Default is code"subject".
- `date_var`: Character string specifying the visit or assessment date variable. Default is code"date".
- `session_var`: Character string specifying the session identifier variable used to order fragmented records. Default is code"session".
- `keep_session`: Logical. If codeTRUE, the first session value encountered within each group is retained. Default is codeTRUE.
- `session_name`: Character string specifying the name of the retained session variable. Default is code"first_session".
- `n_rows_name`: Character string specifying the name of the variable recording the number of rows merged. Default is code"n_rows_collapsed".
- `arrange_desc_session`: Logical. If codeTRUE, records are ordered by descending session number before merging. Default is codeFALSE.
- `empty_strings_to_na`: Logical. If codeTRUE, empty character strings are converted to missing values prior to merging. Default is codeTRUE.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A data frame containing one row per unique combination of codeid_var and codedate_var. Additional variables may include: describe itemn_rows_collapsedNumber of fragmented rows merged. itemfirst_sessionFirst session value retained, if codekeep_session = TRUE.

**See also:** None documented.

## `MultivariableRegressionTable`

**Purpose:** Multivariable regression table

**Canonical usage**
```r
MultivariableRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = TRUE,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  FDR = TRUE,
  FDRAlpha = 0.05,
  Method = c("lm", "ridge", "lasso", "elasticnet"),
  CVFolds = 10,
  Lambda = c("lambda.min", "lambda.1se"),
  Seed = 123,
  MissingDataStrategy = c("drop_sparse_impute", "impute", "complete_cases",
    "drop_sparse_complete_cases"),
  MaxMissingPredictor = 0.3,
  ImputeMethod = c("median_mode"),
  MinCompleteCases = NULL,
  outcome_modes = "auto",
  reference_levels = NULL,
  binary_subsets = NULL
)
```

**Description:** Fit one multivariable regression model per outcome and return a stable, label-aware regression object for downstream tables, diagnostics, and plots.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing outcomes, predictors, and covariates.
- `outcome_vars`: Character vector of outcome variable names.
- `predictor_vars`: Character vector of predictor variable names.
- `covariates`: Optional character vector of covariate variable names. Covariates are treated as mandatory adjustments: for penalized methods (code"ridge", code"lasso", code"elasticnet") they are exempted from the penalty (codepenalty.factor = 0), so they are never shrunk or selected out of the model.
- `Standardize`: Logical. If codeTRUE, ordinary models are fit on standardized continuous variables for the primary estimate. Standardized coefficients are always calculated separately regardless of this setting.
- `Relabel`: Logical. If codeTRUE, use variable labels from codesjlabelled when available.
- `TreatOrdinalAs`: How ordinal outcomes and predictors are handled.
- `FDR`: Logical. If codeTRUE, calculate FDR-adjusted p-values for ordinary regression terms.
- `FDRAlpha`: Numeric FDR threshold retained in metadata.
- `Method`: Regression method. One of code"lm", code"ridge", code"lasso", or code"elasticnet".
- `CVFolds`: Number of cross-validation folds for penalized models.
- `Lambda`: Lambda selection rule for penalized models. One of code"lambda.min" or code"lambda.1se".
- `Seed`: Random seed used for deterministic cross-validation folds.
- `MissingDataStrategy`: Missing-data handling strategy. The default, code"drop_sparse_impute", drops sparse predictors and covariates, then imputes remaining predictor missingness.
- `MaxMissingPredictor`: Maximum allowed missingness proportion for predictors and covariates before they are dropped by sparse-drop strategies. Default is code0.30.
- `ImputeMethod`: Imputation method for predictor/covariate missingness. Currently code"median_mode": median for numeric variables and mode for factor, character, and logical variables.
- `MinCompleteCases`: Optional minimum number of modeling rows required after missing-data handling.
- `outcome_modes`: Multi-category outcome strategy. Supply a single code"auto" or a named character vector whose values are code"auto", code"multinomial", code"ordinal", code"one_vs_rest", code"binary_subset", or code"skip". In automatic mode, ordered factors use proportional-odds regression and unordered factors use multinomial regression.
- `reference_levels`: Optional named character vector giving reference levels for categorical outcomes. Unspecified outcomes use their first retained factor level.
- `binary_subsets`: Optional named list. Each outcome assigned code"binary_subset" must have exactly two level names here, ordered as reference then event.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `OutcomeVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `PredictorVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `Covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A named list with stable components: codeModels, codeFormattedTable, codeLargeTable, codeRegressionMatrix, codeVariableImportanceMatrix, codePredictions, codeDiagnostics, codeModelSummary, codeMulticollinearity, codePlots, and codeMetadata. codeFormattedTable is a report-facing codegt table grouped by outcome, matching the style of codelink[=MakeUnivariateRegressionTable]MakeUnivariateRegressionTable(): predictor rows only, a combined verbEstimate (95% CI) cell, and bold significant p-values. codeLargeTable is a data frame holding the full per-term detail (including covariate rows) for programmatic use, plus an codeAliased flag marking perfectly collinear terms the model dropped. codeModelSummary reports per-outcome codeConverged, codeSeparationDetected, and codeAliasedTermCount. For ordinary (code"lm") logistic fits, quasi-complete separation is detected (fitted probabilities pinned at 0/1, exploded standardized coefficients, or non-convergence); the affected model's estimates are blanked (codeNA) and codeConverged is set to codeFALSE so unreliable coefficients do not propagate into tables or plots. codeModelSummary also carries an omnibus model test per outcome (codeModelStat, codeModelStatType, codeModelPValue): an F-test for linear models and a likelihood-ratio test for logistic models (codeNA for penalized fits, which have no valid classical omnibus test). codePlots contains ggplot objects built from the stored result tables and predictions without refitting models; the coefficient heatmap uses robust, clamped fill limits so a single extreme value cannot dominate the scale, and each outcome column is annotated at the top with its omnibus p-value (ordinary models) or cross-validated deviance explained (penalized models) to discourage interpreting coefficients from a model that is not significant overall. Multi-category outcomes add explicit codeOutcomeLevel, codeReferenceLevel, codeContrast, codeComparisonLabel, and codeOutcomeMode fields. Unordered factors use nominal multinomial models by default. Ordered factors use proportional-odds models; their odds ratios describe movement toward a higher category, conditional on the predictors. One-vs-rest models are available for level-specific scientific questions, but their overlapping comparisons should be interpreted with multiplicity in mind. Binary subsets change both the analysis population and estimand. The resolved strategy, reference, engine, class counts, and concise scientific advice are recorded under codeMetadata$Outcomes and codeMetadata$ModelingAdvice.

**See also:** None documented.

## `Pipeline_SOM_MClust`

**Purpose:** SOM + latent profile clustering pipeline (with AHP and distance baselines)

**Canonical usage**
```r
CreateClusterModel_SOM_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize", "explore"),
  k_range = 2:10,
  models = c(1, 2, 3, 6),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  seed_som = 934521L,
  seed_lpa = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_var = NULL,
  lpa_progress = FALSE,
  lpa_em_itmax = 100L,
  lpa_em_tol = 1e-05,
  lpa_timeout_seconds = 120,
  lpa_drop_zero_sd = TRUE,
  lpa_zero_sd_tol = 1e-08,
  skip_model_after_n_failures = 2L,
  slow_fit_seconds = 120,
  min_nodes_per_cluster = 5,
  high_dist_quantile = 0.95,
  low_prob_threshold = 0.7,
  stability_resamples = 0L,
  stability_seed = 934522L,
  stability_progress = FALSE,
  stability_cores = NULL,
  .NodeClusterFn = NULL
)

Pipeline_SOM_MClust(...)

Pipeline_SOMClust(...)

CreateSOMClusterModel(...)
```

**Description:** End-to-end pipeline to: itemize item Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied Z-score object. item Fit a Self-Organizing Map (SOM; kohonen) on complete cases. item Generate aweSOM visualizations (Circular, Line, Cloud) with optional relabeling using variable labels from the original data frame. item Cluster SOM codebook vectors using latent profile analysis (tidyLPA / mclust backend). item In codemethod = "exploratory", fit a grid of models and select a recommended solution using an Analytic Hierarchy Process (AHP)-style index combining AIC, BIC, and Entropy. item In codemethod = "finalize", fit a user-specified model and number of profiles. item Map node-level clusters and posterior probabilities back to individuals. item Store training variable summaries used later to quantify whether projected cohorts fall outside the original training range. This supports a train once, project many clinical phenotyping workflow: learn phenotype structure in a training cohort, then project new cohorts into the fixed phenotype space without reclustering. Ideal use: correlated continuous clinical or biomarker measures where a topology-preserving map is clinically informative before model-based profiles. Missing data: itemize item SOM and clustering are fit only on rows with complete Z-scores. item The returned codeDataWithClusters has exactly the original rows and columns plus one cluster column; rows not used in SOM/LPA get NA. item The returned codeProbFit$individual is also full length, preserving one row per input row with NA posterior probabilities for rows excluded from SOM/LPA. Z-score behavior: itemize item codeZScoreType = "Center and Scale"/"Center Only"/"Scale Only" computes Z-scores from codedf via codeCreateZScoreObject(). item codeZScoreType = "ZScoreObj" projects Z-scores using an external codeZScoreObj via codeProjectZScore(). item codeZScoreType = "PreZScored" uses existing Z-score columns in codedf as-is and does not re-zscore. Readable SOM + Mclust workflow wrapper for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust(). Compatibility wrapper for codelink[=Pipeline_SOM_MClust]Pipeline_SOM_MClust(). Deprecated alias for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateClusterModel_SOM_MClust`, `Pipeline_SOMClust`, `CreateSOMClusterModel`

**Arguments**
- `data`: Data frame containing the variables to be used in SOM and clustering.
- `variables`: Optional character vector of variable names. If NULL, numeric variables are auto-detected using codeSciDataReportR::getNumVars(df, Ordinal = FALSE). In codeZScoreType = "PreZScored", this can also be NULL if you supply codeZScoreVars or if Z-score columns can be auto-detected by prefix.
- `method`: One of code"exploratory" (default) or code"finalize". In code"exploratory", a grid of models is fit and AHP chooses the recommended solution. In code"finalize", the user must specify codefinal_k and codefinal_model.
- `k_range`: Integer vector of numbers of clusters/profiles to consider in exploratory mode. Default code2:10.
- `models`: Integer vector of model specifications for tidyLPA's mclust backend. Model 1 uses equal variance and zero covariance; model 2 uses varying variance and zero covariance; model 3 uses equal variance and equal covariance; and model 6 uses varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are cluster-specific. Supported values and the default are codec(1, 2, 3, 6). Models 4 and 5 require OpenMx and are intentionally unsupported.
- `final_k`: Integer; number of profiles for codemethod = "finalize".
- `final_model`: Integer; model specification for codemethod = "finalize" (should be one of codemodels).
- `ClusterVariableName`: Name of the cluster column in the output. Defaults to code"Cluster". If this column already exists in codedf, it is overwritten (with a message).
- `ZScoreType`: One of: itemize item code"Center and Scale" (default) item code"Center Only" item code"Scale Only" item code"ZScoreObj" (use an existing ZScore object) item code"PreZScored" (use existing Z-score columns in df as-is)
- `ZScoreObject`: Optional ZScoreObj (from codeCreateZScoreObject() or codeProjectZScore()) to use when codeZScoreType = "ZScoreObj".
- `som_xdim, som_ydim`: Optional integers for SOM grid dimensions. If NULL, a square grid with side length codeceiling(n_complete^(1/3)) is used.
- `som_topo`: SOM topology for codekohonen::somgrid(), default code"hexagonal".
- `som_neigh`: SOM neighbourhood function, default code"gaussian".
- `seed_som, seed_lpa`: Integer seeds for SOM and LPA steps (defaults 934521 and 93421).
- `Relabel`: Logical; if TRUE (default), aweSOM plots are relabeled using variable labels from the emphoriginal codedf (via Hmisc or sjlabelled when available) by stripping the Z-score prefix.
- `ZScorePrefix`: Character prefix used for Z-score columns when codeZScoreType = "PreZScored". Default code"Z_".
- `ZScoreVars`: Optional character vector of Z-score column names to use when codeZScoreType = "PreZScored". If NULL, the function attempts to infer them from codevariables or by detecting columns starting with codeZScorePrefix.
- `id_var`: Optional character scalar. If provided and present in codedf, this column is carried into codeProbFit$individual for convenience.
- `lpa_progress`: Logical; if TRUE, print short progress messages while fitting model/profile combinations.
- `lpa_em_itmax`: Integer; maximum number of EM iterations passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_em_tol`: Numeric; EM convergence tolerance passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_timeout_seconds`: Optional timeout in seconds for individual LPA fits. Use NULL to disable timeouts.
- `lpa_drop_zero_sd`: Logical; if TRUE, remove SOM code dimensions with near-zero standard deviation before LPA.
- `lpa_zero_sd_tol`: Numeric tolerance used when codelpa_drop_zero_sd = TRUE.
- `skip_model_after_n_failures`: Optional integer; skip a model family after this many failures.
- `slow_fit_seconds`: Optional runtime threshold used to flag slow LPA fits in diagnostics.
- `min_nodes_per_cluster`: Optional minimum average SOM nodes per cluster considered before attempting a candidate profile count.
- `high_dist_quantile`: Numeric value between 0 and 1 used to define high SOM-distance flags from the training distance distribution. Default is code0.95.
- `low_prob_threshold`: Numeric posterior probability threshold used to flag uncertain phenotype membership. Default is code0.70.
- `stability_resamples`: Number of 90% participant subsample refits used to assess reproducibility for every successful exploratory candidate. Subsamples are drawn without replacement and reuse the reference model's resolved SOM grid dimensions. Defaults to code0 (disabled); use code50 for an exploratory stability screen.
- `stability_seed`: Integer seed for participant subsampling.
- `stability_progress`: Logical; if TRUE, print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `id_col`: strongDeprecated (since 19.15.0). Use codeid_var instead.
- `.NodeClusterFn`: Internal. A function taking the SOM codebook matrix and returning a list with a codenode_cluster integer vector (one label per SOM node) and, optionally, codefit_table, codeahp_best_row, coderecommendation, codebest_fit_name, and codefit_plot. When supplied, the SOM codebook is clustered by that function and the latent-profile grid is not fitted. Used by codelink[=CreateClusterModel_SOM_HDBSCAN]CreateClusterModel_SOM_HDBSCAN(); not part of the user-facing API.
- `...`: Arguments passed to codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Returns:** A list of class code"Pipeline_SOM_MClust" with components: itemize item codemethod, codevars_used, codeZScoreType, codeZScoreObject, codeZScoreVars, codeClusterVariableName item codeDataWithClusters: original codedf with only the cluster column appended item codefit_plot: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model (plus reproducibility when subsample stability is enabled) item codeModelInfo_SOM: list with codesom_model, codesom_codes, codesom_grid, codetraining_variable_summary, codeSOMFit (distance diagnostics, baselines, and per-cluster flags), codeplots (aweSOM plots) item codeModelInfo_MClust: list with codelpa_models, codefit_table, codeAHP information, and codediagnostics for LPA warnings, failures, runtimes, and preprocessing item codeModelInfo_MClust$Stability: subsample replicate, cluster recovery, and summary tables when codestability_resamples > 0 item codeProbFit: list with codenode (node-level posterior probabilities), codeindividual (full-length per-person mapping and probabilities), and probability plots

**See also:** None documented.

## `Pipeline_SOMClust`

**Purpose:** SOM + latent profile clustering pipeline (with AHP and distance baselines)

**Canonical usage**
```r
CreateClusterModel_SOM_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize", "explore"),
  k_range = 2:10,
  models = c(1, 2, 3, 6),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  seed_som = 934521L,
  seed_lpa = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_var = NULL,
  lpa_progress = FALSE,
  lpa_em_itmax = 100L,
  lpa_em_tol = 1e-05,
  lpa_timeout_seconds = 120,
  lpa_drop_zero_sd = TRUE,
  lpa_zero_sd_tol = 1e-08,
  skip_model_after_n_failures = 2L,
  slow_fit_seconds = 120,
  min_nodes_per_cluster = 5,
  high_dist_quantile = 0.95,
  low_prob_threshold = 0.7,
  stability_resamples = 0L,
  stability_seed = 934522L,
  stability_progress = FALSE,
  stability_cores = NULL,
  .NodeClusterFn = NULL
)

Pipeline_SOM_MClust(...)

Pipeline_SOMClust(...)

CreateSOMClusterModel(...)
```

**Description:** End-to-end pipeline to: itemize item Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied Z-score object. item Fit a Self-Organizing Map (SOM; kohonen) on complete cases. item Generate aweSOM visualizations (Circular, Line, Cloud) with optional relabeling using variable labels from the original data frame. item Cluster SOM codebook vectors using latent profile analysis (tidyLPA / mclust backend). item In codemethod = "exploratory", fit a grid of models and select a recommended solution using an Analytic Hierarchy Process (AHP)-style index combining AIC, BIC, and Entropy. item In codemethod = "finalize", fit a user-specified model and number of profiles. item Map node-level clusters and posterior probabilities back to individuals. item Store training variable summaries used later to quantify whether projected cohorts fall outside the original training range. This supports a train once, project many clinical phenotyping workflow: learn phenotype structure in a training cohort, then project new cohorts into the fixed phenotype space without reclustering. Ideal use: correlated continuous clinical or biomarker measures where a topology-preserving map is clinically informative before model-based profiles. Missing data: itemize item SOM and clustering are fit only on rows with complete Z-scores. item The returned codeDataWithClusters has exactly the original rows and columns plus one cluster column; rows not used in SOM/LPA get NA. item The returned codeProbFit$individual is also full length, preserving one row per input row with NA posterior probabilities for rows excluded from SOM/LPA. Z-score behavior: itemize item codeZScoreType = "Center and Scale"/"Center Only"/"Scale Only" computes Z-scores from codedf via codeCreateZScoreObject(). item codeZScoreType = "ZScoreObj" projects Z-scores using an external codeZScoreObj via codeProjectZScore(). item codeZScoreType = "PreZScored" uses existing Z-score columns in codedf as-is and does not re-zscore. Readable SOM + Mclust workflow wrapper for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust(). Compatibility wrapper for codelink[=Pipeline_SOM_MClust]Pipeline_SOM_MClust(). Deprecated alias for codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateClusterModel_SOM_MClust`, `Pipeline_SOM_MClust`, `CreateSOMClusterModel`

**Arguments**
- `data`: Data frame containing the variables to be used in SOM and clustering.
- `variables`: Optional character vector of variable names. If NULL, numeric variables are auto-detected using codeSciDataReportR::getNumVars(df, Ordinal = FALSE). In codeZScoreType = "PreZScored", this can also be NULL if you supply codeZScoreVars or if Z-score columns can be auto-detected by prefix.
- `method`: One of code"exploratory" (default) or code"finalize". In code"exploratory", a grid of models is fit and AHP chooses the recommended solution. In code"finalize", the user must specify codefinal_k and codefinal_model.
- `k_range`: Integer vector of numbers of clusters/profiles to consider in exploratory mode. Default code2:10.
- `models`: Integer vector of model specifications for tidyLPA's mclust backend. Model 1 uses equal variance and zero covariance; model 2 uses varying variance and zero covariance; model 3 uses equal variance and equal covariance; and model 6 uses varying variance and varying covariance. Zero-covariance models assume conditional independence between variables within each cluster. Equal parameters are shared across clusters; varying parameters are cluster-specific. Supported values and the default are codec(1, 2, 3, 6). Models 4 and 5 require OpenMx and are intentionally unsupported.
- `final_k`: Integer; number of profiles for codemethod = "finalize".
- `final_model`: Integer; model specification for codemethod = "finalize" (should be one of codemodels).
- `ClusterVariableName`: Name of the cluster column in the output. Defaults to code"Cluster". If this column already exists in codedf, it is overwritten (with a message).
- `ZScoreType`: One of: itemize item code"Center and Scale" (default) item code"Center Only" item code"Scale Only" item code"ZScoreObj" (use an existing ZScore object) item code"PreZScored" (use existing Z-score columns in df as-is)
- `ZScoreObject`: Optional ZScoreObj (from codeCreateZScoreObject() or codeProjectZScore()) to use when codeZScoreType = "ZScoreObj".
- `som_xdim, som_ydim`: Optional integers for SOM grid dimensions. If NULL, a square grid with side length codeceiling(n_complete^(1/3)) is used.
- `som_topo`: SOM topology for codekohonen::somgrid(), default code"hexagonal".
- `som_neigh`: SOM neighbourhood function, default code"gaussian".
- `seed_som, seed_lpa`: Integer seeds for SOM and LPA steps (defaults 934521 and 93421).
- `Relabel`: Logical; if TRUE (default), aweSOM plots are relabeled using variable labels from the emphoriginal codedf (via Hmisc or sjlabelled when available) by stripping the Z-score prefix.
- `ZScorePrefix`: Character prefix used for Z-score columns when codeZScoreType = "PreZScored". Default code"Z_".
- `ZScoreVars`: Optional character vector of Z-score column names to use when codeZScoreType = "PreZScored". If NULL, the function attempts to infer them from codevariables or by detecting columns starting with codeZScorePrefix.
- `id_var`: Optional character scalar. If provided and present in codedf, this column is carried into codeProbFit$individual for convenience.
- `lpa_progress`: Logical; if TRUE, print short progress messages while fitting model/profile combinations.
- `lpa_em_itmax`: Integer; maximum number of EM iterations passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_em_tol`: Numeric; EM convergence tolerance passed to codemclust::emControl(). Use NULL to leave mclust defaults unchanged.
- `lpa_timeout_seconds`: Optional timeout in seconds for individual LPA fits. Use NULL to disable timeouts.
- `lpa_drop_zero_sd`: Logical; if TRUE, remove SOM code dimensions with near-zero standard deviation before LPA.
- `lpa_zero_sd_tol`: Numeric tolerance used when codelpa_drop_zero_sd = TRUE.
- `skip_model_after_n_failures`: Optional integer; skip a model family after this many failures.
- `slow_fit_seconds`: Optional runtime threshold used to flag slow LPA fits in diagnostics.
- `min_nodes_per_cluster`: Optional minimum average SOM nodes per cluster considered before attempting a candidate profile count.
- `high_dist_quantile`: Numeric value between 0 and 1 used to define high SOM-distance flags from the training distance distribution. Default is code0.95.
- `low_prob_threshold`: Numeric posterior probability threshold used to flag uncertain phenotype membership. Default is code0.70.
- `stability_resamples`: Number of 90% participant subsample refits used to assess reproducibility for every successful exploratory candidate. Subsamples are drawn without replacement and reuse the reference model's resolved SOM grid dimensions. Defaults to code0 (disabled); use code50 for an exploratory stability screen.
- `stability_seed`: Integer seed for participant subsampling.
- `stability_progress`: Logical; if TRUE, print subsample progress messages.
- `stability_cores`: Number of workers for stability refits. codeNULL uses all detected physical cores minus one, capped at the requested resamples and any scheduler limit. Each worker holds a refit in memory, so lower this setting for large SOM or high-dimensional analyses.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `id_col`: strongDeprecated (since 19.15.0). Use codeid_var instead.
- `.NodeClusterFn`: Internal. A function taking the SOM codebook matrix and returning a list with a codenode_cluster integer vector (one label per SOM node) and, optionally, codefit_table, codeahp_best_row, coderecommendation, codebest_fit_name, and codefit_plot. When supplied, the SOM codebook is clustered by that function and the latent-profile grid is not fitted. Used by codelink[=CreateClusterModel_SOM_HDBSCAN]CreateClusterModel_SOM_HDBSCAN(); not part of the user-facing API.
- `...`: Arguments passed to codelink[=CreateClusterModel_SOM_MClust]CreateClusterModel_SOM_MClust().

**Returns:** A list of class code"Pipeline_SOM_MClust" with components: itemize item codemethod, codevars_used, codeZScoreType, codeZScoreObject, codeZScoreVars, codeClusterVariableName item codeDataWithClusters: original codedf with only the cluster column appended item codefit_plot: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model (plus reproducibility when subsample stability is enabled) item codeModelInfo_SOM: list with codesom_model, codesom_codes, codesom_grid, codetraining_variable_summary, codeSOMFit (distance diagnostics, baselines, and per-cluster flags), codeplots (aweSOM plots) item codeModelInfo_MClust: list with codelpa_models, codefit_table, codeAHP information, and codediagnostics for LPA warnings, failures, runtimes, and preprocessing item codeModelInfo_MClust$Stability: subsample replicate, cluster recovery, and summary tables when codestability_resamples > 0 item codeProbFit: list with codenode (node-level posterior probabilities), codeindividual (full-length per-person mapping and probabilities), and probability plots

**See also:** None documented.

## `Plot2GroupStats`

**Purpose:** Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE by p; COLOR by Category (vector or data frame); stable point size; palette via paletteer)

**Canonical usage**
```r
Plot2GroupStats(
  data,
  variables,
  VariableCategories = NULL,
  impClust,
  normalClust,
  group_var,
  missing_threshold = 0.8,
  max_levels = 10,
  label_q = 0.05,
  x_axis = c("signed_logp", "signed_effect", "effect", "logp"),
  sort_by = c("q", "p", "effect", "signed_logp", "signed_effect", "none"),
  mct_args = list(),
  palette = NULL,
  point_size = 3.5
)
```

**Description:** Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE by p; COLOR by Category (vector or data frame); stable point size; palette via paletteer)

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: data.frame
- `variables`: character vector of variables to analyze
- `VariableCategories`: optional: itemize item data frame with columns Variable, Category; OR item vector of categories (named by variable OR unnamed aligned to codeVariables)
- `impClust, normalClust`: labels for the two groups (impClust plotted to the RIGHT for signed axes)
- `group_var`: column name in codeData holding the group labels
- `missing_threshold`: drop vars with > this fraction missing (default 0.80)
- `max_levels`: drop factors with > this many levels (default 10)
- `label_q`: label threshold using q (default 0.05)
- `x_axis`: one of c("signed_logp","signed_effect","effect","logp")
- `sort_by`: one of c("q","p","effect","signed_logp","signed_effect","none")
- `mct_args`: list of extra args to SciDataReportR::MakeComparisonTable(); e.g., AddEffectSize=TRUE
- `palette`: Optional paletteer palette string for category colors. When codeNULL (the default), the SciDataReportR palette is used. Passing a paletteer string such as code"pals::alphabet" still works as before.
- `point_size`: numeric constant for point size (default 3.5)
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `GroupVar`: strongDeprecated (since 19.15.0). Use codegroup_var instead.

**Returns:** list(plot=ggplot, table=gtsummary, pvaltable=data.frame, data_used=tibble)

**See also:** None documented.

## `PlotAnovaRelationshipsMatrix`

**Purpose:** Plot ANOVA Relationships Matrix

**Canonical usage**
```r
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
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** This function plots the relationship between continuous and categorical variables using ANOVA or Kruskal-Wallis tests. It generates a "heatmap" with points colored and shaped based on statistical significance and effect size. When generalized eta-squared is unavailable, raw p-values are colored with codelink[=scale_color_pvalue]scale_color_pvalue(); the FDR-corrected plot uses adjusted p-values.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame containing the variables of interest.
- `CatVars`: Character vector of categorical variable names.
- `ContVars`: Character vector of continuous variable names.
- `covariates`: Optional character vector of covariate names for ANCOVA analysis.
- `Relabel`: Logical indicating whether to relabel variables with their labels (default is TRUE).
- `Parametric`: Logical indicating whether to use parametric (ANOVA) or non-parametric (Kruskal-Wallis) tests (default is TRUE).
- `Ordinal`: Logical, indicating whether ordinal variables should be considered.
- `min_n`: Minimum number of complete observations required for a tested relationship.
- `eps`: Small positive value used to avoid zero-size plotting artifacts.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all p-values at once (historical behavior). code"per_outcome" corrects separately within each outcome: outcomes are the continuous variables (codeContVars).
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing three ggplot objects: p (scatter plot without multiple comparison correction), p_FDR (scatter plot with FDR correction), and pvaltable (data frame of p-values and significance).

**See also:** None documented.

## `PlotAssociations`

**Purpose:** Plot Associations

**Canonical usage**
```r
PlotAssociations(
  data,
  Var1,
  Var2,
  Ordinal = FALSE
)
```

**Description:** This function generates scatter plots or box plots to visualize the relationship between two variables.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame containing the variables of interest.
- `Var1`: The name of the first variable.
- `Var2`: The name of the second variable.
- `Ordinal`: Logical, indicating whether ordinal variables should be included.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A ggplot object representing the relationship between the variables.

**See also:** None documented.

## `PlotBetaProfile`

**Purpose:** Plot standardized beta profiles across continuous predictors

**Canonical usage**
```r
PlotBetaProfile(
  data,
  predictor_vars,
  outcome_var,
  covariates = NULL,
  VariableCategories = NULL,
  Sort = c("original", "pvalue", "fdr", "effect", "within_category_pvalue",
    "within_category_effect"),
  AdjustMethod = "fdr",
  Alpha = 0.05,
  Relabel = TRUE,
  codebook = NULL,
  RemoveXAxisLabels = TRUE,
  InteractiveLabels = TRUE
)
```

**Description:** Fits a separate linear model for every predictor and displays the resulting standardized regression coefficients with 95% confidence intervals. This is the continuous-outcome companion to codelink[=PlotZScore]PlotZScore(): use it to compare the direction and magnitude of associations across a panel of predictors.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the outcome, predictors, and covariates.
- `predictor_vars`: Nonempty character vector of numeric predictor names.
- `outcome_var`: Character string naming a numeric continuous outcome.
- `covariates`: Optional character vector of covariate names. Covariates are included without automatic standardization.
- `VariableCategories`: Optional categories for the predictors. Supply a vector corresponding to codepredictor_vars, a named vector keyed by predictor name, or a data frame with codeVariable and codeCategory columns.
- `Sort`: Variable ordering. One of code"original", code"pvalue", code"fdr", code"effect", code"within_category_pvalue", or code"within_category_effect".
- `AdjustMethod`: Multiple-testing method passed to codelink[stats:p.adjust]stats::p.adjust().
- `Alpha`: Numeric significance threshold recorded in the returned metadata. Significance does not control plot color.
- `Relabel`: Logical. If codeTRUE, display labels are resolved from the supplied codebook, then variable label attributes, with variable names as the fallback.
- `codebook`: Optional data frame containing codeVariable and codeLabel.
- `RemoveXAxisLabels`: Logical. If codeTRUE, x-axis labels are hidden.
- `InteractiveLabels`: Logical. If codeTRUE, the point layer contains a codetext aesthetic for codeplotly::ggplotly(..., tooltip = "text").

**Returns:** A named list with three elements: describe itemcodePlotA ggplot object showing standardized betas and 95% CIs. itemcodeResultsTableA tibble with one successfully analyzed predictor per row and columns codeVariable, codeLabel, codeCategory, codeBeta, codeSE, codeCILow, codeCIHigh, codePValue, codeFDR, codeN, codeR, codeAdjustedR, and codeTooltip. itemcodeMetadataA named list describing the outcome, requested and analyzed predictors, covariates, adjustment method, alpha threshold, sorting mode, confidence level, and number of fitted models.

**See also:** None documented.

## `PlotBlandAltman`

**Purpose:** Plot Bland-Altman Agreement Plot

**Canonical usage**
```r
PlotBlandAltman(
  data,
  Variable1,
  Variable2
)
```

**Description:** Generates a Bland-Altman plot to visualize the agreement between two variables.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the variables to compare.
- `Variable1`: The name of the first variable (as a string) to compare.
- `Variable2`: The name of the second variable (as a string) to compare.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A list containing: itemplotA ggplot2 object of the Bland-Altman plot. itemstatsA list of Bland-Altman statistics from the BlandAltmanLeh package.

**See also:** None documented.

## `PlotCategoricalDistributions`

**Purpose:** Plot categorical distributions

**Canonical usage**
```r
PlotCategoricalDistributions(
  data,
  variables = NULL,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  LabelType = "percent",
  MissingLabel = "Missing"
)
```

**Description:** This function creates plots to visualize the distributions of categorical variables in a dataframe.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataframe containing the variables to be plotted.
- `variables`: Optional. A character vector specifying the names of the categorical variables to be plotted. If NULL, categorical variables are automatically detected.
- `Relabel`: Logical. If TRUE, missing labels in the dataframe are replaced with column names as labels for plotting.
- `Ordinal`: Deprecated logical compatibility option; use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. This categorical plot accepts code"Categorical" or code"Exclude".
- `LabelType`: Character. Either "percent" or "count", indicating what should be shown on the x-axis and inside the bars.
- `MissingLabel`: Character label to use for missing values.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A ggplot object visualizing the distributions of categorical variables.

**See also:** None documented.

## `PlotCatInteractionEffectsMatrix`

**Purpose:** Plot Categorical Interaction Effects Matrix

**Canonical usage**
```r
PlotCatInteractionEffectsMatrix(
  data,
  predictor_vars,
  outcome_vars = NULL,
  xVarLabels = NULL,
  yVarLabels = NULL,
  interVar,
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** This function calculates and visualizes the interaction effects between categorical variables.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset containing the variables.
- `predictor_vars`: A character vector of the names of the x-axis categorical variables.
- `outcome_vars`: A character vector of the names of the y-axis categorical variables. Defaults to NULL, in which case it takes the same values as xVars.
- `xVarLabels`: A character vector of labels for the x-axis variables. Defaults to NULL, in which case it takes the same values as xVars.
- `yVarLabels`: A character vector of labels for the y-axis variables. Defaults to NULL, in which case it takes the same values as yVars.
- `interVar`: The name of the interaction variable.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all interaction p-values at once (historical behavior). code"per_outcome" corrects separately within each y-axis variable (codeoutcome_vars).
- `yVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.

**Returns:** A list containing matrices of interaction coefficients, p-values, ggplot objects for visualizations, and tables of FDR-corrected p-values.

**See also:** None documented.

## `PlotChiSqCovar`

**Purpose:** Plot Chi-Square Tests for Categorical Associations (optionally stratified by covariates)

**Canonical usage**
```r
PlotChiSqCovar(
  data,
  predictor_vars,
  outcome_vars,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  min_n = 4,
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** Conducts Chi-square tests between sets of categorical variables and visualizes the results. NOTE: Chi-square tests do not natively "adjust" for covariates. If codecovars are provided, this function can (optionally) run tests emphwithin strata (each combination of covariate levels), and combine p-values across strata (Fisher's method) for a single summary p-value per pair. If you need true covariate adjustment, use regression-based models (logistic/multinomial).

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data.frame containing the dataset.
- `predictor_vars`: Character vector of x-axis categorical variables.
- `outcome_vars`: Character vector of y-axis categorical variables. If NULL, uses xVars.
- `covariates`: Optional character vector of covariate variables used for stratification (not adjustment).
- `Relabel`: Logical; whether to use variable labels (sjlabelled) in the plot.
- `Ordinal`: Logical; included for backward compatibility (currently unused here).
- `min_n`: Minimum number of complete observations required for a tested association.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `yVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all p-values at once (historical behavior). code"per_outcome" corrects separately within each outcome: outcomes are the y-axis variables (codeoutcome_vars).
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list with: itempggplot for unadjusted p-values itempvaltablewide table of unadjusted p-values itemp_FDRggplot for FDR-adjusted p-values itempvaltable_FDRwide table of FDR-adjusted p-values itemdetailslong table with diagnostics (n, warnings, strata info)

**See also:** None documented.

## `PlotClusterAssignment`

**Purpose:** Plot a two-dimensional cluster review map

**Canonical usage**
```r
PlotClusterMap(
  data,
  x,
  y,
  ClusterVar = "Cluster",
  centroids = NULL,
  title = "Cluster review map",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  noise_label = 0L
)

PlotClusterAssignment(data, x, y, ClusterVar = "Cluster")
```

**Description:** Cluster assignments in a two-dimensional review space. The space is frozen at training time and reused for projection so training and projected cases are directly comparable. It is a display space only and does not affect the clustering itself.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotClusterMap`

**Arguments**
- `data`: Data frame containing coordinates and a cluster column.
- `x, y`: Coordinate variable names.
- `ClusterVar`: Cluster column name.
- `centroids`: Optional data frame of centroid coordinates in the same two columns, overlaid as crosses.
- `title, subtitle, xlab, ylab`: Plot annotations.
- `noise_label`: Cluster value treated as noise, or codeNULL to disable.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterBoxplot`

**Purpose:** Plot cluster boxplots by variable

**Canonical usage**
```r
PlotClusterBoxplot(
  data,
  ClusterVar,
  variables,
  codebook = NULL,
  Scale = FALSE,
  ScoreType = c("auto", "z", "t", "raw"),
  ReferenceLines = c("auto", "z", "t", "none"),
  ClusterLabel = c("n_percent", "n", "none"),
  Relabel = TRUE,
  FillTitle = "Test",
  Palette = NULL,
  YLabel = NULL,
  BoxplotWidth = 0.75,
  OutlierSize = 0.8,
  BaseSize = 14
)
```

**Description:** Create grouped cluster boxplots where clusters are shown on the x-axis and selected variables are shown as filled boxplots within each cluster. This is useful for visualizing cognitive, clinical, biomarker, or domain profiles across discovered clusters or subgroup solutions.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `ClusterVar`: Character string naming the cluster/grouping variable.
- `variables`: Character vector of variable names to plot.
- `codebook`: Optional codebook data frame with columns codeVariable and codeLabel.
- `Scale`: Logical. If codeTRUE, variables are z-scored across all participants before plotting. Default is codeFALSE.
- `ScoreType`: Character. One of code"auto", code"z", code"t", or code"raw". Used for axis labeling and reference-line defaults.
- `ReferenceLines`: Character. One of code"auto", code"z", code"t", or code"none". If code"auto", z-score reference lines are added only when codeScale = TRUE. If code"z", lines are added at code-1, code-0.5, code0, code0.5, and code1. If code"t", lines are added at code40, code45, code50, code55, and code60. If code"none", no reference lines are added so users can customize overlays manually using additional ggplot layers.
- `ClusterLabel`: Character. One of code"n_percent", code"n", or code"none". code"n_percent" adds cluster sample size and percent to the x-axis label. code"n" adds only sample size. code"none" uses only the cluster name/value.
- `Relabel`: Logical. If codeTRUE, variable labels are used when available. Labels are pulled first from codeCodebook, then from variable label attributes. Default is codeTRUE.
- `FillTitle`: Character string used as the fill legend title. Default is code"Test".
- `Palette`: Optional character vector of colors used for the fill scale. If codeNULL, a stable 20-color default palette is used. If more than 20 variables are plotted, colors are interpolated from the default palette so the function does not fail.
- `YLabel`: Optional y-axis label. If codeNULL, an appropriate label is chosen automatically from codeScale and codeScoreType.
- `BoxplotWidth`: Numeric width passed to codeggplot2::geom_boxplot(). Default is code0.75.
- `OutlierSize`: Numeric size of outlier points. Default is code0.8.
- `BaseSize`: Base font size for the plot theme. Default is code14.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `Codebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterCentreHeatmap`

**Purpose:** Plot cluster centre profiles as a heatmap

**Canonical usage**
```r
PlotClusterCentreHeatmap(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre",
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

PlotClusterCentreProfile(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre"
)
```

**Description:** Cluster centres across clustering variables. Values are the centres in the frozen analysis scale, so a centred and scaled model reads directly as standard deviations from the cohort mean. codePlotClusterCentreProfile() shows the same centres as connected lines, which reads more like the SOM line map when variables have a meaningful order.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotClusterCentreProfile`

**Arguments**
- `centers`: Matrix or data frame of cluster centres, one row per cluster.
- `variable_labels`: Optional display labels for the columns.
- `title`: Plot title.
- `value_label`: Legend title describing the centre scale.

**Returns:** A codeggplot object, or codeNULL when centres are unavailable.

**See also:** None documented.

## `PlotClusterCentreProfile`

**Purpose:** Plot cluster centre profiles as a heatmap

**Canonical usage**
```r
PlotClusterCentreHeatmap(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre",
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

PlotClusterCentreProfile(
  centers,
  variable_labels = NULL,
  title = "Cluster centre profiles",
  value_label = "Centre"
)
```

**Description:** Cluster centres across clustering variables. Values are the centres in the frozen analysis scale, so a centred and scaled model reads directly as standard deviations from the cohort mean. codePlotClusterCentreProfile() shows the same centres as connected lines, which reads more like the SOM line map when variables have a meaningful order.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotClusterCentreHeatmap`

**Arguments**
- `centers`: Matrix or data frame of cluster centres, one row per cluster.
- `variable_labels`: Optional display labels for the columns.
- `title`: Plot title.
- `value_label`: Legend title describing the centre scale.

**Returns:** A codeggplot object, or codeNULL when centres are unavailable.

**See also:** None documented.

## `PlotClusterComposition`

**Purpose:** Plot categorical composition by cluster

**Canonical usage**
```r
PlotClusterComposition(
  data,
  variables,
  cluster,
  facet_by = c("variable", "cluster"),
  style = c("stacked", "enrichment")
)
```

**Description:** Plot categorical composition by cluster

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing the categorical variables.
- `variables`: Categorical variable names.
- `cluster`: Cluster assignment vector aligned to codedata.
- `facet_by`: Whether stacked-bar facets represent categorical variables (the default) or clusters.
- `style`: Either code"stacked" for composition bars or code"enrichment" for a cluster-first heatmap of percentage-point differences from the cohort.

**Returns:** A codeggplot object, or codeNULL when no categorical variables are given.

**See also:** None documented.

## `PlotClusterDiagnostic`

**Purpose:** Plot a per-cluster diagnostic value

**Canonical usage**
```r
PlotClusterDiagnostic(individual, value, title = NULL, noise_label = 0L)
```

**Description:** Boxplot of any per-participant diagnostic (posterior probability, distance to centroid, outlier score) split by cluster.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `individual`: Per-participant diagnostic table containing codeCluster.
- `value`: Diagnostic column name.
- `title`: Optional plot title.
- `noise_label`: Cluster value treated as noise, or codeNULL to disable.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterFitReview`

**Purpose:** Plot cluster fit-review metrics

**Canonical usage**
```r
PlotClusterFitReview(
  fit_table,
  x = "Classes",
  metrics = NULL,
  group = NULL,
  title = "Candidate model review"
)
```

**Description:** Plot cluster fit-review metrics

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `fit_table`: Candidate-model fit table.
- `x`: Candidate-count column name.
- `metrics`: Numeric metric columns to display. When codeNULL, a concise set of decision-relevant raw metrics is selected automatically.
- `group`: Optional candidate-family column, such as code"Model" for Mclust or code"Epsilon" for HDBSCAN. When codeNULL, it is inferred when possible.
- `title`: Plot title.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterMap`

**Purpose:** Plot a two-dimensional cluster review map

**Canonical usage**
```r
PlotClusterMap(
  data,
  x,
  y,
  ClusterVar = "Cluster",
  centroids = NULL,
  title = "Cluster review map",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  noise_label = 0L
)

PlotClusterAssignment(data, x, y, ClusterVar = "Cluster")
```

**Description:** Cluster assignments in a two-dimensional review space. The space is frozen at training time and reused for projection so training and projected cases are directly comparable. It is a display space only and does not affect the clustering itself.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotClusterAssignment`

**Arguments**
- `data`: Data frame containing coordinates and a cluster column.
- `x, y`: Coordinate variable names.
- `ClusterVar`: Cluster column name.
- `centroids`: Optional data frame of centroid coordinates in the same two columns, overlaid as crosses.
- `title, subtitle, xlab, ylab`: Plot annotations.
- `noise_label`: Cluster value treated as noise, or codeNULL to disable.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterProfiles`

**Purpose:** Plot labelled numeric profiles by cluster

**Canonical usage**
```r
PlotClusterProfiles(data, ClusterVar, variables, codebook = NULL)
```

**Description:** Plot labelled numeric profiles by cluster

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `ClusterVar`: Character string naming the cluster/grouping variable.
- `variables`: Character vector of variable names to plot.
- `codebook`: Optional codebook data frame with columns codeVariable and codeLabel.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotClusterSilhouette`

**Purpose:** Plot a per-participant silhouette profile

**Canonical usage**
```r
PlotClusterSilhouette(silhouette, title = "Silhouette profile")
```

**Description:** The classic silhouette profile: one bar per participant, sorted within cluster, with the average silhouette width marked. Bars near zero or below sit closer to a neighbouring cluster than their own.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `silhouette`: A codecluster::silhouette() object or a matrix with codecluster, codeneighbor, and codesil_width columns.
- `title`: Plot title.

**Returns:** A codeggplot object, or codeNULL when silhouette widths are unavailable.

**See also:** None documented.

## `PlotContinuousDistributions`

**Purpose:** Plot Continuous Distributions

**Canonical usage**
```r
PlotContinuousDistributions(
  data,
  variables = NULL,
  Fill = NULL,
  Relabel = TRUE,
  FacetLabelStyle = c("both", "label_only", "variable_only", "auto"),
  ncol = 3,
  TreatOrdinalAs = "Categorical"
)
```

**Description:** Creates rain-cloud plots (half-violin + box/median + scatter) for one or more continuous variables, with optional group-wise colouring.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the variables to be plotted.
- `variables`: Character vector of column names to plot.
- `Fill`: Optional column name for grouping.
- `Relabel`: Logical; use variable labels when available.
- `FacetLabelStyle`: One of "both", "label_only", "variable_only", "auto".
- `ncol`: Number of columns in the facet grid.
- `Ordinal`: Deprecated logical compatibility option; use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. code"Continuous" includes their numeric score; code"Exclude" omits them. code"Both" is not meaningful for this plot and errors.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A ggplot object.

**See also:** None documented.

## `PlotCorrelationComparisons`

**Purpose:** Compare correlations between two independent groups

**Canonical usage**
```r
PlotCorrelationComparisons(
  data,
  predictor_vars = NULL,
  outcome_vars = NULL,
  group_var,
  comparison_group = NULL,
  reference_group = NULL,
  covariates = NULL,
  method = "pearson",
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  min_n = 4,
  eps = 1e-12,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  reversal_style = c("outline", "stripe", "none"),
  interactive = c("none", "plotly", "girafe", "both"),
  low_color = "#B2182B",
  mid_color = "white",
  high_color = "#2166AC",
  color_limits = c(-2, 2),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Computes correlations or partial correlations separately within two independent groups, compares corresponding correlations, and visualizes the between-group difference as a heatmap. Tile color represents codeDeltaR = r_comparison - r_reference, significance stars represent the statistical test comparing the two correlations, and striped tiles indicate correlations with opposite signs between groups.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `predictor_vars`: Character vector of predictor variables. If codeNULL, variable selection is inherited from codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap().
- `outcome_vars`: Optional character vector of outcome variables. If codeNULL, the same variables are used on both axes.
- `group_var`: Character string naming the grouping variable.
- `comparison_group`: Optional level of codegroup_var used as the comparison group. Positive DeltaR values indicate a more positive correlation in this group relative to the reference group.
- `reference_group`: Optional level of codegroup_var used as the reference group.
- `covariates`: Optional character vector of covariates used to calculate partial correlations within each group.
- `method`: Correlation method. Either code"pearson" or code"spearman".
- `Relabel`: Logical indicating whether variable labels should be used when available.
- `TreatOrdinalAs`: Passed to codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap().
- `min_n`: Minimum number of complete observations required for an individual correlation.
- `eps`: Variance tolerance passed to codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap().
- `fdr_scope`: Scope for FDR correction of correlation-comparison tests. One of code"matrix", code"per_outcome", or code"per_predictor".
- `reversal_style`: How correlations with opposite signs should be shown. One of code"outline" (default), code"stripe", or code"none". Stripes are a static-only display option and require codeggpattern.
- `interactive`: Optional interactive output. One of code"none" (default), code"plotly", code"girafe", or code"both". Static ggplots are always retained in codeUnadjusted$plot and codeFDRCorrected$plot; requested widgets are added under codeInteractive.
- `low_color`: Color representing negative DeltaR values.
- `mid_color`: Color representing DeltaR = 0.
- `high_color`: Color representing positive DeltaR values.
- `color_limits`: Limits for the DeltaR color scale. The theoretical range is -2 to 2.

**Returns:** A list containing: describe itemCorrelationsThe original codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap() objects for the comparison and reference groups. itemUnadjustedMatrices and heatmap using raw comparison p-values. itemFDRCorrectedMatrices and heatmap using FDR-adjusted comparison p-values. itemResultsA tibble with one row per correlation pair. itemDirectionReversalLogical matrix indicating opposite correlation signs between groups. itemMetadataComparison settings, group information, and the inferential approximation used. itemInteractiveOptional Plotly and/or ggiraph widgets.

**See also:** None documented.

## `PlotCorrelationsHeatmap`

**Purpose:** Plot correlations heatmap

**Canonical usage**
```r
PlotCorrelationsHeatmap(
  data,
  predictor_vars = NULL,
  outcome_vars = NULL,
  covariates = NULL,
  method = "pearson",
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  min_n = 3,
  eps = 1e-12,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Computes correlations or partial correlations and plots a heatmap. Handles: itemize item continuous + categorical covariates item labelled data item non-syntactic names item sparse real-world datasets item ordinal variables item partial correlations via residualization

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: data.frame
- `predictor_vars`: character vector
- `outcome_vars`: character vector
- `covariates`: optional covariates
- `method`: pearson/spearman/kendall
- `Relabel`: use labels
- `Ordinal`: Deprecated logical compatibility option; use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. code"Continuous" includes ordinal scores and code"Exclude" omits them.
- `min_n`: minimum complete rows
- `eps`: variance tolerance
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). With code"matrix", FDR correction is applied across the whole p-value matrix at once (historical behavior). With code"per_outcome", correction is applied separately within each outcome: in this function outcomes are the columns of the p-value matrix, i.e. codeoutcome_vars (codeoutcome_margin = 2).
- `cluster_rows, cluster_columns`: Logical; cluster predictor rows and/or outcome columns from their displayed correlation profiles. Defaults retain the caller-supplied variable order.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `yVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list. codeUnadjusted and codeFDRCorrected each contain coder, codep, codenpairs, and codeplot. The standardized aliases codep (same as codeUnadjusted) and codep_fdr (same as codeFDRCorrected) are also included.

**See also:** None documented.

## `PlotDatasetComparison`

**Purpose:** Plot dataset comparison diagnostics

**Canonical usage**
```r
PlotDatasetComparison(
  CompareObj,
  Plot = c("All", "Checks", "SummaryMetrics", "StructureChanges", "VariableChanges",
    "TopChangedVariables"),
  interactive = TRUE,
  TopN = 10
)
```

**Description:** Create diagnostic plots from a codeCompareDatasets() result object. This function visualizes dataset-version changes, including check status, summary metrics, structure changes, and variable-level value changes.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `CompareObj`: A list returned by codeCompareDatasets().
- `Plot`: Character value specifying which plot to return. Options are code"All", code"Checks", code"SummaryMetrics", code"StructureChanges", code"VariableChanges", and code"TopChangedVariables". Default is code"All".
- `interactive`: Logical; if codeTRUE, plots are converted to interactive codeplotly objects using codeplotly::ggplotly(). Default is codeTRUE.
- `TopN`: Integer number of variables or records to preview in plots and hover text. Default is code10.
- `Interactive`: strongDeprecated (since 19.15.0). Use codeinteractive instead.

**Returns:** If codePlot = "All", a named list of plots. Otherwise, a single plot object. Plot objects are either codeggplot objects or codeplotly htmlwidgets, depending on codeInteractive.

**See also:** None documented.

## `PlotDirectionalHeatmaps`

**Purpose:** Create directional heatmaps across continuous & binary variables

**Canonical usage**
```r
PlotDirectionalHeatmaps(
  data,
  variables = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Combines: itemize item Continuous~Continuous (Pearson/Spearman) item Binary~Binary (Phi; 1 == PositiveLevel) item Binary~Continuous (r_pb; 1 == PositiveLevel) into a single square heatmap with raw and FDR-star overlays.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe.
- `variables`: Character vector of variables to include (subset of codedata columns). The analysis is symmetric: every variable is related to every other, so a single variable set defines both axes. If NULL, uses all detected continuous + binary vars.
- `Relabel`: Logical; use sjlabelled variable labels if present.
- `Ordinal`: Logical; passed to codelink[=PlotPointCorrelationsHeatmap]PlotPointCorrelationsHeatmap() for the binary~continuous block, where it controls whether ordinal variables are treated as continuous. Defaults to codeTRUE.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", threaded through to the three sub-analyses (codelink[=PlotCorrelationsHeatmap]PlotCorrelationsHeatmap(), codelink[=PlotPhiHeatmap]PlotPhiHeatmap(), codelink[=PlotPointCorrelationsHeatmap]PlotPointCorrelationsHeatmap()). Correction is applied within each sub-analysis block (continuous~continuous, binary~binary, binary~continuous), matching historical behavior; each sub-function's documented outcome orientation applies within its block.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `yVars`: strongDeprecated (since 19.15.0). Use codevariables instead. If supplied, the old rectangular x-by-y display is still honored.

**Returns:** list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)

**See also:** None documented.

## `plotForestFromTable`

**Purpose:** Create a Forest Plot from Univariate Regression Tables

**Canonical usage**
```r
PlotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)

plotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)
```

**Description:** This function generates a forest plot from the results of codelink[=MakeUnivariateRegressionTable]MakeUnivariateRegressionTable(). codeplotForestFromTable() was renamed to codePlotForestFromTable() in SciDataReportR 20.5.0 to match the package's verbPlot* naming convention. It remains available as a backwards-compatible synonym.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `PlotForestFromTable`

**Arguments**
- `UnivariateRegressionTables`: Either the full list returned by codelink[=MakeUnivariateRegressionTable]MakeUnivariateRegressionTable() (its codeResults dataframe is used directly), or a dataframe with the codeResults columns. Passing a dataframe lets you filter, reorder, or relabel codeResults before plotting; required columns are codeOutcomeLabel, codeTermLabel, codeEstimate, codeConfLow, codeConfHigh, and codePValue (codeSignificant and codeReferenceValue are recomputed if absent). Lists created by older package versions (without a codeResults element) are still supported.
- `pSize`: Numeric. Size of the points in the plot. Default is 2.
- `Flip`: Logical. If codeFALSE, outcomes are facets and predictors/terms are rows. If codeTRUE, predictors/terms are facets and outcomes are rows.

**Returns:** A ggplot object representing the forest plot.

**See also:** None documented.

## `PlotForestFromTable`

**Purpose:** Create a Forest Plot from Univariate Regression Tables

**Canonical usage**
```r
PlotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)

plotForestFromTable(UnivariateRegressionTables, pSize = 2, Flip = FALSE)
```

**Description:** This function generates a forest plot from the results of codelink[=MakeUnivariateRegressionTable]MakeUnivariateRegressionTable(). codeplotForestFromTable() was renamed to codePlotForestFromTable() in SciDataReportR 20.5.0 to match the package's verbPlot* naming convention. It remains available as a backwards-compatible synonym.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `plotForestFromTable`

**Arguments**
- `UnivariateRegressionTables`: Either the full list returned by codelink[=MakeUnivariateRegressionTable]MakeUnivariateRegressionTable() (its codeResults dataframe is used directly), or a dataframe with the codeResults columns. Passing a dataframe lets you filter, reorder, or relabel codeResults before plotting; required columns are codeOutcomeLabel, codeTermLabel, codeEstimate, codeConfLow, codeConfHigh, and codePValue (codeSignificant and codeReferenceValue are recomputed if absent). Lists created by older package versions (without a codeResults element) are still supported.
- `pSize`: Numeric. Size of the points in the plot. Default is 2.
- `Flip`: Logical. If codeFALSE, outcomes are facets and predictors/terms are rows. If codeTRUE, predictors/terms are facets and outcomes are rows.

**Returns:** A ggplot object representing the forest plot.

**See also:** None documented.

## `PlotInteractionEffectsContinuous`

**Purpose:** Plot Single Interaction Effect

**Canonical usage**
```r
PlotInteractionEffectsContinuous(
  data,
  interVar = NULL,
  outcome_var = NULL,
  predictor_var = NULL,
  covariates = NULL,
  n_lines = 3,
  alpha = 0.6,
  point_size = 2
)
```

**Description:** Creates a scatter plot with regression lines showing the interaction between a predictor and outcome variable, moderated by either a continuous or categorical variable.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the variables to be analyzed
- `interVar`: Character string specifying the interaction variable (moderator)
- `outcome_var`: Character string specifying the outcome variable
- `predictor_var`: Character string specifying the predictor variable
- `covariates`: Character vector of covariate names to include in the model
- `n_lines`: For continuous moderators, number of lines to plot (default: 3 for low/med/high)
- `alpha`: Transparency level for points (default: 0.6)
- `point_size`: Size of points (default: 2)
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `outcomeVar`: strongDeprecated (since 19.15.0). Use codeoutcome_var instead.
- `predictorVar`: strongDeprecated (since 19.15.0). Use codepredictor_var instead.
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A ggplot object showing the interaction effect

**See also:** None documented.

## `PlotInteractionEffectsMatrix`

**Purpose:** Plot Interaction Effects Matrix

**Canonical usage**
```r
PlotInteractionEffectsMatrix(
  data,
  interVar = NULL,
  outcome_vars = NULL,
  predictor_vars = NULL,
  covariates = NULL,
  Relabel = TRUE,
  TreatOrdinalAs = "Exclude",
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** Creates a heatmap visualization of interaction effects between continuous variables, showing whether interactions result in slope reversals or maintain the same direction.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the variables to be analyzed
- `interVar`: Character string specifying the interaction variable (moderator). Can be categorical or continuous.
- `outcome_vars`: Character vector of outcome variable names (displayed on rows)
- `predictor_vars`: Character vector of predictor variable names (displayed on columns)
- `covariates`: Character vector of covariate names to include in the models
- `Relabel`: Logical indicating whether to use variable labels if available (default: TRUE)
- `Ordinal`: strongDeprecated (since 20.20.0). Use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. This interaction matrix accepts code"Exclude", code"Continuous", or code"Categorical". Categorical ordinal outcomes are not supported by the linear models.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `outcomeVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `predictorVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all interaction p-values at once (historical behavior). code"per_outcome" corrects separately within each outcome variable (codeoutcome_vars).
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing: describe itemUnadjustedList with unadjusted results including: itemize item C: Matrix of interaction coefficients item S: Matrix of slope direction indicators (1 = same, -1 = reversed) item P: Matrix of p-values item D: Matrix of interaction coefficient signs item Slope1: Matrix of slopes at low values of interVar (mean - 1SD for continuous, reference group for categorical) item Slope2: Matrix of slopes at high values of interVar (mean + 1SD for continuous, comparison group for categorical) item plot: ggplot object of the heatmap item pvaltable: P-value table in wide format itemFDRCorrectedList with FDR-corrected results (same structure as Unadjusted) itemRelabelLogical indicating whether relabeling was applied itemCovariatesCharacter vector of covariates used iteminterVarThe interaction variable name itemraw_dataProcessed data frame with all calculated values

**See also:** None documented.

## `PlotMergeValidation`

**Purpose:** Plot merge validation diagnostics

**Canonical usage**
```r
PlotMergeValidation(
  MergeObj,
  Plot = c("All", "Checks", "Coverage", "JoinAudit", "Agreement", "Conflicts"),
  interactive = TRUE
)
```

**Description:** Create diagnostic plots from a codeValidateMerge() result object. This function visualizes key merge-audit outputs, including validation check status, key coverage, join-variable auditing, duplicate-variable agreement, and duplicate-variable conflict counts. Use this after running codeValidateMerge() to quickly inspect whether a merged dataset appears trustworthy.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `MergeObj`: A list returned by codeValidateMerge().
- `Plot`: Character value specifying which plot to return. Options are code"All", code"Checks", code"Coverage", code"JoinAudit", code"Agreement", and code"Conflicts". Default is code"All".
- `interactive`: Logical; if codeTRUE, plots are converted to interactive codeplotly objects using codeplotly::ggplotly(). Default is codeTRUE.
- `Interactive`: strongDeprecated (since 19.15.0). Use codeinteractive instead.

**Returns:** If codePlot = "All", a named list of plots. Otherwise, a single plot object. Plot objects are either codeggplot objects or codeplotly htmlwidgets, depending on codeInteractive.

**See also:** None documented.

## `PlotMiningMatrix`

**Purpose:** PlotMiningMatrix

**Canonical usage**
```r
PlotMiningMatrix(
  data,
  outcome_vars,
  predictor_vars = NULL,
  covariates = NULL,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  Parametric = TRUE,
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** Generate a matrix of statistical relationships between variables.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `outcome_vars`: Outcome variables.
- `predictor_vars`: Predictor variables. If NULL, uses OutcomeVars.
- `covariates`: Optional covariates (reserved for future use).
- `Relabel`: Use labels instead of names.
- `TreatOrdinalAs`: How ordinal variables are handled: code"Categorical", code"Continuous", code"Both", or code"Exclude".
- `Parametric`: Use parametric tests.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `OutcomeVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `PredictorVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all pairwise p-values at once (historical behavior, computed on the symmetrized pair table). code"per_outcome" corrects separately within each x-axis variable (codeXVar, ordered by codeoutcome_vars).
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** List with tables and plots.

**See also:** None documented.

## `PlotMissingData`

**Purpose:** Plot Missing Data

**Canonical usage**
```r
PlotMissingData(
  data,
  variables = NULL,
  HoverVars = NULL,
  x_var = NULL,
  facet_by = NULL,
  Relabel = TRUE,
  show_perc = TRUE,
  show_perc_var = TRUE,
  cluster = FALSE
)
```

**Description:** Visualize missing data patterns with variables as rows and observations as columns. Optional hover variables can be included to facilitate quality control workflows when converting the plot to an interactive Plotly figure.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of variables to visualize. If NULL, all columns except codeHoverVars, codex_var, and codefacet_by are used.
- `HoverVars`: Optional character vector of columns to include in hover text. Useful for participant IDs, visit names, dates, sites, etc.
- `x_var`: Optional single column name to use for the x-axis. Numeric and date variables retain their original scale; categorical variables use a discrete axis. Missing x values are displayed as code"Missing".
- `facet_by`: Optional single column name used to create missingness panels. Missing facet values are displayed in a code"Missing" panel.
- `Relabel`: Logical. If TRUE, variable labels are used when available.
- `show_perc`: Logical. If TRUE, overall missingness percentages are shown in the legend.
- `show_perc_var`: Logical. If TRUE, variable-specific missingness percentages are appended to y-axis labels.
- `cluster`: Logical. If TRUE, variables are clustered by missingness pattern.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A ggplot object.

**See also:** None documented.

## `PlotNumInteractionEffectsMatrix`

**Purpose:** Plot Numerical Interaction Effects Matrix

**Canonical usage**
```r
PlotNumInteractionEffectsMatrix(
  data,
  predictor_vars,
  outcome_vars = NULL,
  xVarLabels = NULL,
  yVarLabels = NULL,
  interVar = NULL,
  covariates = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor")
)
```

**Description:** This function calculates interaction effects between numerical variables and plots them as matrices.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset containing the variables.
- `predictor_vars`: A character vector of the names of the x-axis numerical variables.
- `outcome_vars`: A character vector of the names of the y-axis numerical variables. Defaults to NULL.
- `xVarLabels`: A character vector of labels for the x-axis variables. Defaults to NULL.
- `yVarLabels`: A character vector of labels for the y-axis variables. Defaults to NULL.
- `interVar`: The interaction variable. Defaults to NULL.
- `covariates`: A character vector of the names of covariate variables. Defaults to NULL.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `yVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all interaction p-values at once (historical behavior). code"per_outcome" corrects separately within each y-axis variable (codeoutcome_vars).
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing matrices, ggplot objects for visualizations, and tables of p-values.

**See also:** None documented.

## `PlotPartialRegressionScatter`

**Purpose:** Partial Regression Plot

**Canonical usage**
```r
PlotPartialRegressionScatter(
  data,
  IndepVar,
  DepVar,
  covariates = NULL,
  Relabel = TRUE
)
```

**Description:** Generate a partial regression plot for a specified independent and dependent variable while adjusting for covariates. In addition to the figure, key parameters such as the correlation method, whether relabeling was used, the covariates, R-squared, p-value, and sample size are returned.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset to use.
- `IndepVar`: A string specifying the independent variable.
- `DepVar`: A string specifying the dependent variable.
- `covariates`: A character vector of covariate names for adjustment. Defaults to NULL.
- `Relabel`: Logical indicating whether to use labelled names from the data (using sjlabelled::get_label). Defaults to TRUE.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing: itemplotA ggplot2 object representing the partial regression plot. itemmethodThe correlation method (as provided). itemRelabelLogical; whether relabeling was applied. itemCovariatesThe vector of covariates. itemr2The R-squared of the partial regression model. itemp_valueThe p-value for the independent variable coefficient. itemnThe sample size (number of complete cases). itemequationThe regression equation string.

**See also:** None documented.

## `PlotPathway_KT`

**Purpose:** Plot the kynurenine-tryptophan pathway

**Canonical usage**
```r
PlotPathway_KT(
  results_table,
  title = "",
  value_type = "auto",
  metabolite_mapping = NULL,
  use_fdr = FALSE
)

CreatePathwayPlot_KT(...)
```

**Description:** Creates a pathway diagram for the kynurenine-tryptophan metabolic pathway with color-coded fold changes or correlations and significance indicators. codeCreatePathwayPlot_KT() has been superseded by codePlotPathway_KT(). It remains available as a backwards-compatible alias during the pathway plot's planned transition to a metabolomics-focused package.

**Deprecation status:** Current documented interface.

**Related exported aliases:** `CreatePathwayPlot_KT`

**Arguments**
- `results_table`: Data frame with columns: Metabolite, p_value, p_adj, and either "% Change" or "correlation"
- `title`: Character string for plot title
- `value_type`: Character string: "auto", "fold_change", or "correlation"
- `metabolite_mapping`: Named character vector mapping results table names to standard names. For example: c("N'-Formylkynurenine" = "N-Formylkynurenine", "Quinolinic Acid(log10)" = "Quinolinic Acid")
- `use_fdr`: Logical: if TRUE uses FDR-adjusted p-values (p_adj) for significance, if FALSE uses raw p-values. Default is FALSE.
- `...`: Arguments passed to codelink[=PlotPathway_KT]PlotPathway_KT().

**Returns:** A ggplot2 object

**See also:** None documented.

## `plotPCA`

**Purpose:** Plot PCA scores

**Canonical usage**
```r
plotPCA(
  PCAObj,
  Var = NULL,
  t = "NULL",
  HoverVar = NULL,
  HoverVars = NULL,
  Components = NULL,
  Mode = c("auto", "3D", "2D"),
  ColorType = c("auto", "factor", "continuous"),
  Relabel = TRUE,
  Title = NULL
)
```

**Description:** Create an interactive 2D or 3D Plotly scatter plot from a PCA object created by codeCreatePCAObject().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `PCAObj`: A PCA object returned by codeCreatePCAObject(). Must contain codeScores and codeCombinedData.
- `Var`: Optional character string naming a variable in codePCAObj$CombinedData used to color points. Default is codeNULL.
- `t`: Deprecated compatibility argument for color type. Use codeColorType instead. If set to code"Factor", codeVar is treated as categorical. Any other value treats codeVar as continuous. Default is code"NULL".
- `HoverVar`: Optional character string naming one variable in codePCAObj$CombinedData to display in hover text. Retained for backward compatibility. Prefer codeHoverVars for new code. Default is codeNULL.
- `HoverVars`: Optional character vector naming one or more variables in codePCAObj$CombinedData to display in hover text. If codeNULL, row number is shown. Default is codeNULL.
- `Components`: Optional character or numeric vector specifying which score columns to plot. Supply two components for 2D or three components for 3D. If codeNULL, the first three score columns are used.
- `Mode`: Character. Either code"auto", code"3D", or code"2D". If code"auto", the plot dimension is inferred from codeComponents. If codeComponents = NULL, code"auto" defaults to 3D. Default is code"auto".
- `ColorType`: Character. Either code"auto", code"factor", or code"continuous". code"auto" treats character, factor, logical, and labelled variables as categorical and numeric variables as continuous. Default is code"auto".
- `Relabel`: Logical. If codeTRUE, labels attached to hover variables are used in hover text when available. If codeFALSE, raw variable names are used. Default is codeTRUE.
- `Title`: Optional plot title. Default is codeNULL, which produces no title.

**Returns:** A Plotly htmlwidget.

**See also:** None documented.

## `PlotPhiHeatmap`

**Purpose:** Plot Phi Correlations Between Binary Variables

**Canonical usage**
```r
PlotPhiHeatmap(
  data,
  CatVars,
  Relabel = TRUE,
  binary_map = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Computes pairwise phi coefficients between binary categorical variables with explicit 0/1 coding (1 == PositiveLevel from codecreateBinaryMapping()), then renders heatmap-style plots with raw and FDR-adjusted significance.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe.
- `CatVars`: Character vector of binary categorical variable names.
- `Relabel`: Logical; if TRUE, uses sjlabelled variable labels for axes.
- `binary_map`: Optional mapping as returned by codecreateBinaryMapping(). If NULL, a mapping is created internally for codeCatVars.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all p-values at once (historical behavior). code"per_outcome" corrects separately within each y-axis variable (codeYVar); the Phi matrix is symmetric, so this treats each variable's row of tiles as one family.
- `cluster_rows, cluster_columns`: Logical; cluster y-axis rows and/or x-axis columns using displayed Phi-coefficient profiles.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A list with: itemize item codeUnadjusted: list(PvalTable, plot) item codeFDRCorrected: list(PvalTable, plot) item codemethod = "Phi" item codeRelabel item codeBinaryMapping (used)

**See also:** None documented.

## `PlotPointCorrelationsHeatmap`

**Purpose:** Plot Point-Biserial Correlations Between Binary and Continuous Variables

**Canonical usage**
```r
PlotPointCorrelationsHeatmap(
  data,
  CatVars,
  ContVars,
  covariates = NULL,
  Relabel = TRUE,
  Ordinal = TRUE,
  binary_map = NULL,
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Calculates point-biserial correlations (binary vs continuous) with explicit 0/1 coding where 1 == PositiveLevel from codecreateBinaryMapping(), and renders heatmap-style tiles.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe.
- `CatVars`: Character vector of binary categorical variables.
- `ContVars`: Character vector of continuous variables.
- `covariates`: Optional covariates (reserved).
- `Relabel`: Logical; use sjlabelled variable labels for axes.
- `Ordinal`: Logical; reserved for future use.
- `binary_map`: Optional mapping as returned by codecreateBinaryMapping(). If NULL, a mapping is created internally for codeCatVars.
- `fdr_scope`: Either code"matrix" (default) or code"per_outcome", passed to codelink[=ApplyFDRCorrection]ApplyFDRCorrection(). code"matrix" corrects across all p-values at once (historical behavior). code"per_outcome" corrects separately within each continuous variable: outcomes are the continuous variables (codeContVars).
- `cluster_rows, cluster_columns`: Logical; cluster continuous y-axis rows and/or binary x-axis columns using displayed correlation profiles.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list with Unadjusted, FDRCorrected, method ("R_pb"), Relabel, Covariates, BinaryMapping.

**See also:** None documented.

## `PlotPValueComparisons`

**Purpose:** Plot P-Value Comparisons

**Canonical usage**
```r
PlotPValueComparisons(
  data,
  group_var,
  variables = NULL,
  VariableCategories = NULL,
  Relabel = TRUE
)
```

**Description:** This function generates a plot comparing p-values for different variables across two or more groups.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame containing the variables to compare.
- `group_var`: Character string specifying the name of the column in codeData that contains the group labels.
- `variables`: Character vector specifying the names of the columns in codeData to include in the comparison. If codeNULL, all columns except codeGroupVariable are included.
- `VariableCategories`: Character vector specifying the categories for each variable. If codeNULL, no categories are used.
- `Relabel`: Logical indicating whether to replace missing labels with the column names.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `GroupVariable`: strongDeprecated (since 19.15.0). Use codegroup_var instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.

**Returns:** A ggplot object displaying the p-value comparisons.

**See also:** None documented.

## `plotSigAssociations`

**Purpose:** Plot Significant Associations

**Canonical usage**
```r
plotSigAssociations(
  data,
  AnovaMatrixObject,
  PVar = "p",
  Pthresh = 0.05
)
```

**Description:** Generate ggbetweenstats for significant correlations based on a previously generated anova matrix

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset used to generate the scatterplots.
- `AnovaMatrixObject`: The output of the PlotAnovaRelationshipsMatrix function.
- `PVar`: The name of the column used to filter for significance (default is "P").
- `Pthresh`: The significance threshold (default is 0.05).
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A list of scatterplot objects for significant correlations.

**See also:** None documented.

## `plotSigCorrelations`

**Purpose:** Plot Significant Correlations

**Canonical usage**
```r
plotSigCorrelations(
  data,
  CorrelationHeatmapObject,
  PVar = "P",
  Pthresh = 0.05
)
```

**Description:** Generate scatterplots for significant correlations based on a previously generated correlation heatmap.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset used to generate the scatterplots.
- `CorrelationHeatmapObject`: The output of the PlotCorrelationsHeatmap function.
- `PVar`: The name of the column used to filter for significance (default is "P").
- `Pthresh`: The significance threshold (default is 0.05).
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A list of scatterplot objects for significant correlations.

**See also:** None documented.

## `PlotSpiderChart`

**Purpose:** Plot a spider chart across continuous and binary variables

**Canonical usage**
```r
PlotSpiderChart(
  data,
  variables,
  group_var = NULL,
  Relabel = TRUE,
  ContinuousSummary = "mean",
  ContinuousScaling = "zscore",
  Fill = FALSE,
  FillAlpha = 0.2,
  Facet = FALSE,
  VariableOrder = "input",
  VariableCategories = NULL,
  BinaryPositiveValue = 1,
  Palette = NULL,
  LineSize = 1,
  PointSize = 2,
  ShowPoints = FALSE,
  LegendTitle = NULL,
  PlotTitle = NULL,
  Subtitle = NULL,
  Caption = NULL,
  AxisLabelSize = 12,
  AxisTextSize = 10,
  StripTextSize = 11,
  WrapLabels = TRUE,
  LabelWrapWidth = 22,
  LabelRadiusMultiplier = 1.22,
  PlotMarginTop = 40,
  PlotMarginRight = 120,
  PlotMarginBottom = 40,
  PlotMarginLeft = 120,
  interactive = FALSE,
  InteractiveHeight = 700,
  InteractiveWidth = NULL,
  InteractiveAxisMin = NULL,
  InteractiveAxisMax = NULL,
  tooltip_digits = 2
)
```

**Description:** This function summarizes a set of variables and displays them on a spider chart. Continuous variables are plotted as mean z-scores by default using codeCreateZScoreObject(), while binary variables are plotted as percentages. It can overlay groups on one spider chart or facet by group, optionally fill the polygons, relabel spokes using variable labels, wrap long labels, reorder variables to visually emphasize between-group differences, and optionally return an interactive radar chart using plotly.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of variable names to plot.
- `group_var`: Optional grouping variable name. If NULL, one overall summary profile is plotted.
- `Relabel`: Logical; if TRUE, use variable labels when available.
- `ContinuousSummary`: Character; one of code"mean" or code"median".
- `ContinuousScaling`: Character; one of code"zscore", code"none", or code"minmax".
- `Fill`: Logical; if TRUE, add transparent polygon fills in the static ggplot version and filled polygons in the interactive version.
- `FillAlpha`: Numeric transparency for fills.
- `Facet`: Logical; if TRUE and codeGroupVariable is supplied, facet by group instead of overlaying all groups on one spider chart. Ignored when codeinteractive = TRUE.
- `VariableOrder`: Character; one of code"input", code"discrimination", code"hierarchical", code"greedy", or code"category_discrimination".
- `VariableCategories`: Optional character vector of categories for codeVariables. Must be the same length as codeVariables when supplied.
- `BinaryPositiveValue`: Optional positive value to use for non-factor binary variables. Defaults to code1. For factor variables, the second factor level is used.
- `Palette`: Optional character name of an codehcl.colors() palette. When codeNULL (the default), the SciDataReportR palette is used. Passing a name such as code"Dark 3" still works exactly as before.
- `LineSize`: Numeric line width for the static ggplot version.
- `PointSize`: Numeric point size for the static ggplot version.
- `ShowPoints`: Logical; if TRUE, show points at each spoke in the static ggplot version.
- `LegendTitle`: Optional legend title. Defaults to codeGroupVariable.
- `PlotTitle`: Optional plot title.
- `Subtitle`: Optional plot subtitle.
- `Caption`: Optional plot caption.
- `AxisLabelSize`: Numeric axis text size for spoke labels in the static ggplot version.
- `AxisTextSize`: Numeric text size for radial axis labels in the static ggplot version.
- `StripTextSize`: Numeric facet strip text size in the static ggplot version.
- `WrapLabels`: Logical; if TRUE, wrap long spoke labels.
- `LabelWrapWidth`: Numeric wrap width passed to codestringr::str_wrap().
- `LabelRadiusMultiplier`: Numeric multiplier controlling how far labels sit outside the spider in the static ggplot version.
- `PlotMarginTop`: Numeric top plot margin for the static ggplot version.
- `PlotMarginRight`: Numeric right plot margin for the static ggplot version.
- `PlotMarginBottom`: Numeric bottom plot margin for the static ggplot version.
- `PlotMarginLeft`: Numeric left plot margin for the static ggplot version.
- `interactive`: Logical; if TRUE, return an interactive plotly radar chart instead of a static ggplot. Default is codeFALSE.
- `InteractiveHeight`: Numeric height in pixels for the interactive widget.
- `InteractiveWidth`: Optional width passed to plotly layout. Defaults to NULL.
- `InteractiveAxisMin`: Optional numeric minimum for the interactive radial axis. If NULL, auto-detected from the summarized values.
- `InteractiveAxisMax`: Optional numeric maximum for the interactive radial axis. If NULL, auto-detected from the summarized values.
- `tooltip_digits`: Integer number of digits to show in interactive tooltips.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `GroupVariable`: strongDeprecated (since 19.15.0). Use codegroup_var instead.
- `TooltipDigits`: strongDeprecated (since 19.15.0). Use codetooltip_digits instead.
- `MakeInteractive`: strongDeprecated (since 19.15.0). Use codeinteractive instead.

**Returns:** A ggplot object when codeinteractive = FALSE, otherwise a plotly htmlwidget.

**See also:** None documented.

## `PlotSplitViolin`

**Purpose:** Split violin with aligned half-boxplots, significance label, sample sizes, and label-aware title

**Canonical usage**
```r
PlotSplitViolin(
  data,
  Var,
  group_var,
  covariates = NULL,
  nonparametric = FALSE,
  annotation_text = NULL,
  show_ns = FALSE,
  plot_title = NULL,
  use_var_label_as_title = FALSE,
  show_n = TRUE,
  n_position = c("legend", "top"),
  n_size = 3.5,
  left_group = NULL,
  color_palette = NULL,
  box_offset = 0.11,
  box_width = 0.15,
  star_from = c("quantile", "data_max", "whisker"),
  star_quantile = 0.995,
  star_pad = 0.05,
  star_size = 6,
  p_label = c("stars", "p_value", "both"),
  ...
)
```

**Description:** Draws a split (left/right) violin for up to two groups at a single x-position, overlays per-group boxplots aligned with each half, and optionally annotates the plot with a p-value significance label. Supports displaying sample sizes (n) and automatically using variable labels for axis and title when available.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing codeVar, codeGroup, and optional covariates.
- `Var`: Numeric outcome variable (tidy-eval).
- `group_var`: Grouping variable (<= 2 unique values).
- `covariates`: Character vector of covariates (default codeNULL).
- `nonparametric`: Logical. If codeFALSE, uses linear model + emmeans contrast. If codeTRUE, uses Wilcoxon (with residualization if covariates are present).
- `annotation_text`: Optional manual annotation (e.g., "*", "ns").
- `show_ns`: Logical; if codeTRUE, display "ns" for non-significant results.
- `plot_title`: Optional custom plot title.
- `use_var_label_as_title`: Logical; if codeTRUE, uses variable label as title.
- `show_n`: Logical; if codeTRUE, display sample size per group.
- `n_position`: Where to display n: "legend" or "top".
- `n_size`: Text size for n labels when shown on top.
- `left_group`: Optional group to force on left side.
- `color_palette`: Optional named vector of colors.
- `box_offset`: Horizontal offset for boxplots.
- `box_width`: Width of boxplots.
- `star_from`: Position method ("quantile","data_max","whisker").
- `star_quantile`: Quantile used for placement.
- `star_pad`: Padding above anchor.
- `star_size`: Text size for annotation.
- `p_label`: P-value annotation style: significance code"stars" (default), exact code"p_value", or code"both".
- `...`: Additional arguments reserved for future extensions.
- `Group`: strongDeprecated (since 19.15.0). Use codegroup_var instead.
- `covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A ggplot2 object. The exact test result is available from codeattr(plot, "comparison").

**See also:** None documented.

## `PlotSwimmerTransitions`

**Purpose:** Plot swimmer-style transitions for a binary condition over repeated visits

**Canonical usage**
```r
PlotSwimmerTransitions(
  data,
  id_var,
  time_var,
  status_var,
  date_var = NULL,
  participant_subset = NULL,
  max_participants = NULL,
  order_participants_by = c("first_positive", "first_transition", "ever_positive",
    "ever_positive_then_burden", "input_order", "n_visits", "n_positive", "pct_positive"),
  x_axis_type = c("visit", "date", "time_from_baseline"),
  time_from_baseline_unit = c("days", "months", "years"),
  show_transition_points = TRUE,
  show_lines = TRUE,
  show_y_axis_labels = FALSE,
  interactive = FALSE,
  plot_title = NULL,
  x_label = NULL,
  y_label = NULL,
  return_data = FALSE
)
```

**Description:** Create a swimmer-style longitudinal plot for a binary condition measured across repeated visits. The plot shows condition status at each visit, highlights transition points where the condition develops or resolves, and supports multiple participant ordering strategies for exploratory quality control, longitudinal debugging, or figure creation.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing repeated observations per participant.
- `id_var`: Unquoted column name identifying the participant.
- `time_var`: Unquoted column name representing visit order, visit number, or time index.
- `status_var`: Unquoted column name representing the binary condition status. Accepted encodings include numeric 0/1, logical TRUE/FALSE, factor values such as code"Yes" and code"No", and character values such as code"1" and code"0".
- `date_var`: Optional unquoted visit date column. This is required when codex_axis_type = "date" or codex_axis_type = "time_from_baseline".
- `participant_subset`: Optional vector of participant IDs to include.
- `max_participants`: Optional maximum number of participants to display after ordering is applied.
- `order_participants_by`: Character string controlling participant order in the plot. Options are code"first_positive", code"first_transition", code"ever_positive", code"ever_positive_then_burden", code"input_order", code"n_visits", code"n_positive", and code"pct_positive".
- `x_axis_type`: Character string indicating whether the x-axis should use aligned visit number (code"visit"), actual calendar date (code"date"), or elapsed time from each participant's baseline date (code"time_from_baseline").
- `time_from_baseline_unit`: Character string specifying the unit for codex_axis_type = "time_from_baseline". Options are code"days", code"months", and code"years".
- `show_transition_points`: Logical. If codeTRUE, transition visits are highlighted.
- `show_lines`: Logical. If codeTRUE, a swimmer line is drawn across visits within each participant.
- `show_y_axis_labels`: Logical. If codeTRUE, show participant labels on the y-axis. Defaults to codeFALSE.
- `interactive`: Logical. If codeTRUE, return an interactive plotly object. Defaults to codeFALSE.
- `plot_title`: Optional custom plot title.
- `x_label`: Optional x-axis label.
- `y_label`: Optional y-axis label. Defaults to codeNULL.
- `return_data`: Logical. If codeTRUE, return both the plot and the processed data.
- `make_interactive`: strongDeprecated (since 19.15.0). Use codeinteractive instead.

**Returns:** A codeggplot object by default. If codeinteractive = TRUE, returns a codeplotly object. If codereturn_data = TRUE, returns a list with: itemize item codeplot: the codeggplot or codeplotly object item codeplot_data: the processed visit-level plotting data item codeparticipant_summary: the participant-level summary table

**See also:** None documented.

## `PlotTimeDistribution`

**Purpose:** Plot Time Distribution

**Canonical usage**
```r
PlotTimeDistribution(
  data,
  DateVariable = "Date"
)
```

**Description:** This function plots the distribution of time-based data.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The data frame containing the time-based data.
- `DateVariable`: The name of the column in the data frame containing the date information. Default is "Date".
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A ggplot object displaying the distribution of time-based data.

**See also:** None documented.

## `PlotTimeSwimmer`

**Purpose:** Plot longitudinal swimmer timelines

**Canonical usage**
```r
PlotTimeSwimmer(
  data,
  id_var,
  Time,
  State = NULL,
  Event = NULL,
  EventType = NULL,
  TimeScale = c("from_first", "observed", "from_event"),
  EventReference = NULL,
  StateInterval = c("forward", "point"),
  Format = c("state_path", "visit_points", "event_rug", "minimal"),
  SortBy = c("duration", "last_time", "first_time", "state", "id"),
  TimeUnit = c("auto", "days", "weeks", "months", "years", "visits"),
  Relabel = TRUE,
  codebook = NULL,
  LineWidth = 5,
  PointSize = 2.5,
  Alpha = 0.9,
  BaseSize = 13
)
```

**Description:** Create swimmer-style longitudinal timeline plots from long-format visit data. Each row in the input data represents a subject visit or observation. The function can visualize longitudinal state trajectories, visit timing, and events over time.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame in long format with one row per visit or observation.
- `id_var`: Character string naming the participant ID column.
- `Time`: Character string naming the time variable.
- `State`: Optional character string naming a state/group variable used for coloring timelines or visit points.
- `Event`: Optional character string naming a logical or binary event indicator variable.
- `EventType`: Optional character string naming an event category variable. Used for event point shapes/colors.
- `TimeScale`: Character. One of code"from_first", code"observed", or code"from_event". code"observed" uses raw observed time values. code"from_first" normalizes each subject relative to their first observed timepoint. code"from_event" normalizes each subject relative to the first occurrence of codeEventReference.
- `EventReference`: Optional event value used when codeTimeScale = "from_event".
- `StateInterval`: Character. One of code"forward" or code"point". code"forward" extends the current state forward until the next visit. code"point" only colors visit points without extending intervals.
- `Format`: Character. One of code"state_path", code"visit_points", code"event_rug", or code"minimal".
- `SortBy`: Character. One of code"duration", code"last_time", code"first_time", code"state", or code"id".
- `TimeUnit`: Character. One of code"auto", code"days", code"weeks", code"months", code"years", or code"visits".
- `Relabel`: Logical. If codeTRUE, use labels from codeCodebook or variable attributes when available.
- `codebook`: Optional codebook data frame with columns codeVariable and codeLabel.
- `LineWidth`: Numeric line width for swimmer segments. Default is code5.
- `PointSize`: Numeric point size for visit/event points. Default is code2.5.
- `Alpha`: Numeric alpha transparency for swimmer segments. Default is code0.9.
- `BaseSize`: Base font size for the plot theme. Default is code13.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `ID`: strongDeprecated (since 19.15.0). Use codeid_var instead.
- `Codebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A codeggplot object.

**See also:** None documented.

## `PlotVolcanoEffects`

**Purpose:** Plot volcano-style association effects

**Canonical usage**
```r
PlotVolcanoEffects(
  data,
  predictor_vars,
  outcome_var,
  covariates = NULL,
  OutcomeType = c("auto", "continuous", "categorical"),
  EffectMetric = c("auto", "cohens_d", "log2fc"),
  AdjustMethod = "fdr",
  Alpha = 0.05,
  Format = c("tiered", "classic", "fdr_only", "directional", "effect_gradient",
    "minimal", "neon"),
  LabelMode = c("none", "top_n", "raw", "significant", "fdr", "extreme"),
  TopN = 10,
  Relabel = TRUE,
  codebook = NULL,
  InteractiveLabels = TRUE,
  ColorBy = NULL
)
```

**Description:** Screen a set of predictor variables against one outcome and visualize the association results as volcano plots. The function supports continuous outcomes using standardized beta values, and two-group categorical outcomes using either Cohen-style standardized effects or log2 fold change effects. Optional covariates are included in the model for each predictor.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `predictor_vars`: Character vector of predictor variable names to screen.
- `outcome_var`: Character string naming the outcome variable.
- `covariates`: Optional character vector of covariate variable names.
- `OutcomeType`: Outcome type. One of code"auto", code"continuous", or code"categorical". If code"auto", numeric outcomes are treated as continuous and nonnumeric two-level outcomes are treated as categorical.
- `EffectMetric`: Effect metric. One of code"auto", code"cohens_d", or code"log2fc". Used for two-group categorical outcomes. For continuous outcomes, the effect is always a standardized beta from a model using scaled predictor and scaled outcome.
- `AdjustMethod`: Multiple-comparison correction method passed to codestats::p.adjust(). Default is code"fdr".
- `Alpha`: Significance threshold for raw and adjusted p-values. Default is code0.05.
- `Format`: Color format. One of code"tiered", code"classic", code"fdr_only", code"directional", code"effect_gradient", code"minimal", or code"neon".
- `LabelMode`: Labeling mode. One of code"none", code"top_n", code"raw", code"significant", code"fdr", or code"extreme". code"raw" labels variables with codePValue < Alpha; code"significant" is retained as an alias for code"raw". code"fdr" labels variables with codeFDR < Alpha. The default is code"none".
- `TopN`: Number of variables to label when codeLabelMode is code"top_n" or code"extreme". Default is code10.
- `Relabel`: Logical. If codeTRUE, variable labels are used when available. Labels are pulled first from codeCodebook if supplied, then from variable label attributes. Default is codeTRUE.
- `codebook`: Optional codebook data frame with columns codeVariable and codeLabel.
- `InteractiveLabels`: Logical. If codeTRUE, a codetext aesthetic is added for compatibility with codeplotly::ggplotly(tooltip = "text"). Default is codeTRUE.
- `ColorBy`: Optional category mapping used to color points in both codeRawPPlot and codeFDRPlot. Supply either a data frame with codeVariable and codeCategory columns or a named atomic vector whose names are predictor variable names and whose values are categories. Tested predictors without a mapping are shown as code"Unmapped" in grey. When codeNULL (the default), the existing significance-based codeFormat colors are used unchanged.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `xVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `yVar`: strongDeprecated (since 19.15.0). Use codeoutcome_var instead.
- `Covariates`: strongDeprecated (since 19.15.0). Use codecovariates instead.
- `Codebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A named list with codeRawPPlot, codeFDRPlot, and codeResultsTable. codeRawPPlot uses code-log10(PValue) on the y-axis. codeFDRPlot uses code-log10(FDR) on the y-axis. codeResultsTable is a tibble with one row per analyzed predictor. For continuous outcomes it includes codeR (the zero-order Pearson correlation between the predictor and outcome) and codeAdjustedR (the covariate-adjusted partial correlation, codeNA when no covariates are given). For two-group categorical outcomes it includes codeGroup1Level, codeGroup2Level, codeGroup1Mean, and codeGroup2Mean (the raw predictor means within each outcome group). These values are also surfaced in the codeTooltip column used by codeplotly::ggplotly(tooltip = "text").

**See also:** None documented.

## `PlotZScore`

**Purpose:** Plot Z-score group differences with statistical significance

**Canonical usage**
```r
PlotZScore(
  data,
  TargetVar,
  variables,
  VariableCategories = NULL,
  Relabel = TRUE,
  sort = TRUE,
  RemoveXAxisLabels = TRUE,
  TreatOrdinalAs = "Continuous",
  Parametric = TRUE,
  SigP_YCoord = 1.5,
  SigFDR_YCoord = 1.6
)

CreateZScorePlot(...)
```

**Description:** This function generates a Z-score plot to compare multiple variables across different groups. It offers options for parametric or non-parametric tests, ordinal treatment, and custom labeling. Significant p-values and FDR-adjusted p-values are highlighted on the plot. codeCreateZScorePlot() has been superseded by codePlotZScore(). It remains available as a backwards-compatible alias and returns the same scientific visualization.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `CreateZScorePlot`

**Arguments**
- `data`: A dataframe containing the data to be analyzed.
- `TargetVar`: A string specifying the column name of the grouping variable.
- `variables`: A vector of strings specifying the column names of the variables to be analyzed.
- `VariableCategories`: An optional vector categorizing the variables.
- `Relabel`: Logical; if TRUE, variables will be relabeled using their labels from the dataframe.
- `sort`: Logical; if TRUE, variables will be sorted by category and p-value.
- `RemoveXAxisLabels`: Logical; if TRUE, X-axis labels will be removed.
- `Ordinal`: strongDeprecated (since 20.20.0). Use codeTreatOrdinalAs instead.
- `TreatOrdinalAs`: How ordinal variables are handled. This numeric plot accepts code"Continuous" or code"Exclude".
- `Parametric`: Logical; if TRUE, parametric tests (t-test/ANOVA) will be used; otherwise, non-parametric tests (Wilcoxon/Kruskal-Wallis) will be used.
- `SigP_YCoord`: Numeric; the y-coordinate for marking significant p-values.
- `SigFDR_YCoord`: Numeric; the y-coordinate for marking significant FDR-adjusted p-values.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Variables`: strongDeprecated (since 19.15.0). Use codevariables instead.
- `...`: Arguments passed to codelink[=PlotZScore]PlotZScore().

**Returns:** A ggplot object representing the Z-score plot.

**See also:** None documented.

## `PrepNumericData`

**Purpose:** Prepare numeric data safely for analysis

**Canonical usage**
```r

```

**Description:** Safely coerces selected variables in a data frame to numeric format while preserving column names and replacing non-finite values (codeInf, code-Inf, codeNaN) with codeNA.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `variables`: Character vector of variable names to process. Defaults to all columns.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A data frame with selected variables converted to numeric and non-finite values replaced with codeNA.

**See also:** None documented.

## `PrepSPSS`

**Purpose:** Prepare a data frame for SPSS export

**Canonical usage**
```r
PrepSPSS(
  data,
  path = NULL,
  name_map_path = NULL,
  label_map_path = NULL,
  return = c("list", "data", "map"),
  quiet = FALSE,
  show_map = FALSE,
  max_length = 64,
  max_label_length = 120,
  compress = "byte",
  ...
)
```

**Description:** codePrepSPSS() converts column names to SPSS-safe variable names while preserving the original column names as variable labels. It also returns and optionally writes a name map so users can track original and exported names.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame or tibble.
- `path`: Optional file path ending in code.sav. If supplied, the prepared data are written using codehaven::write_sav().
- `name_map_path`: Optional file path for saving the name map as a CSV.
- `label_map_path`: Optional file path for saving the label map (all truncated labels alongside their original text) as a CSV.
- `return`: One of code"list", code"data", or code"map".
- `quiet`: Logical. If codeFALSE, prints a compact summary.
- `show_map`: Logical. If codeTRUE, prints the full original-to-SPSS name map. The default is codeFALSE because large scientific datasets can have thousands of renamed variables.
- `max_length`: Maximum SPSS variable-name length. Defaults to 64.
- `max_label_length`: Maximum value-label length in bytes. Defaults to 120, the SPSS limit enforced by codehaven::write_sav(). Factor levels and value labels longer than this are truncated and recorded in the label map.
- `compress`: Compression type passed to codehaven::write_sav().
- `...`: Additional arguments passed to codehaven::write_sav() if codepath is supplied.

**Returns:** Depending on codereturn, either a list (with elements codedata, codename_map, and codelabel_map), the prepared data frame, or the name map.

**See also:** None documented.

## `Project_SOMClust`

**Purpose:** Project cases through a fitted clustering model

**Canonical usage**
```r
ProjectCluster(object, new_df, ...)

Project_SOMClust(...)

ProjectSOMCluster(...)
```

**Description:** Projects codenew_df through the frozen preprocessing, reduction, and clustering layers stored in codeobject. Dispatch is determined by the fitted model's verbPipeline_* class; callers do not select a method-specific projector. Compatibility wrapper for codelink[=ProjectCluster]ProjectCluster(). Deprecated alias for codelink[=ProjectCluster]ProjectCluster().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `ProjectCluster`, `ProjectSOMCluster`

**Arguments**
- `object`: A finalized object returned by verbCreateClusterModel_*().
- `new_df`: New data to project into the frozen cluster structure.
- `...`: Arguments passed to codelink[=ProjectCluster]ProjectCluster().

**Returns:** A method-specific projection result with the common codeProjectionFit contract and codeProbFit assignment table.

**See also:** None documented.

## `Project_ZScore`

**Purpose:** Project standardized scores onto new data using external parameters

**Canonical usage**
```r
ProjectZScore(
  data,
  variables = NULL,
  parameters,
  ParameterInputType = c("df_parameter", "ZScoreObj", "ExternalDataframe"),
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)

Project_ZScore(...)
```

**Description:** codeProject_ZScore() has been superseded by codeProjectZScore(). It remains available as a backwards-compatible alias and returns the same projected Z-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `ProjectZScore`

**Arguments**
- `data`: Data frame on which to project scores.
- `variables`: Character vector; if NULL, project onto all variables for which parameters exist and that are present in df.
- `parameters`: Source of parameters, interpreted by ParameterInputType: itemize item "df_parameter": a data frame with cols Variable, N, Mean, SD item "ZScoreObj": output object from CreateZScoreObject() item "ExternalDataframe": a raw reference data frame; parameters are estimated via CreateZScoreObject() on that frame.
- `ParameterInputType`: One of "df_parameter", "ZScoreObj", "ExternalDataframe".
- `names_prefix`: Prefix for projected variable names.
- `RetainLabels`: Logical; if TRUE and Hmisc available, copy labels from df to new variables.
- `RenameLabels`: Logical; if TRUE, prefix labels the same way as names.
- `center`: Logical; used when ParameterInputType is "df_parameter" or "ExternalDataframe". Ignored for "ZScoreObj" (it uses stored flags).
- `scale`: Logical; same logic as codecenter.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=ProjectZScore]ProjectZScore().

**Returns:** List with same structure as CreateZScoreObject(): itemize item ZScores: projected standardized variables only item DataWithZ: df with projected scores appended item Parameters: parameter data frame actually used for projection item Center, Scale: flags used

**See also:** None documented.

## `ProjectCluster`

**Purpose:** Project cases through a fitted clustering model

**Canonical usage**
```r
ProjectCluster(object, new_df, ...)

Project_SOMClust(...)

ProjectSOMCluster(...)
```

**Description:** Projects codenew_df through the frozen preprocessing, reduction, and clustering layers stored in codeobject. Dispatch is determined by the fitted model's verbPipeline_* class; callers do not select a method-specific projector. Compatibility wrapper for codelink[=ProjectCluster]ProjectCluster(). Deprecated alias for codelink[=ProjectCluster]ProjectCluster().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `Project_SOMClust`, `ProjectSOMCluster`

**Arguments**
- `object`: A finalized object returned by verbCreateClusterModel_*().
- `new_df`: New data to project into the frozen cluster structure.
- `...`: Arguments passed to codelink[=ProjectCluster]ProjectCluster().

**Returns:** A method-specific projection result with the common codeProjectionFit contract and codeProbFit assignment table.

**See also:** None documented.

## `ProjectPCA`

**Purpose:** Project PCA scores onto new data

**Canonical usage**
```r
ProjectPCA(
  data,
  VarsToReduce = NULL,
  PCAInput,
  InputType = c("PCAObj", "LoadingTable"),
  center = TRUE,
  scale = TRUE
)
```

**Description:** Use an existing PCA solution (either a PCA object from CreatePCAObject or a loading table) to compute principal component scores on a new dataset.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: Data frame on which to project PCA scores.
- `VarsToReduce`: Optional character vector of variable names to use. If NULL, uses all variables that appear in both Data and the PCA solution.
- `PCAInput`: Either: itemize item the full object returned by CreatePCAObject (when InputType = "PCAObj"), or item a loading table like CreatePCAObject()$LoadingTable (when InputType = "LoadingTable").
- `InputType`: One of "PCAObj" or "LoadingTable".
- `center`: Logical; only used when InputType is "LoadingTable". For "PCAObj", the centering choice is taken from PCAInput$ScaleParams$center and this argument is ignored (with a warning if it conflicts).
- `scale`: Logical; only used when InputType is "LoadingTable". For "PCAObj", the scaling choice is taken from PCAInput$ScaleParams$scale and this argument is ignored (with a warning if it conflicts).
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A list with: itemScoresData frame of projected PCA scores. itemCombinedDataOriginal Data with scores appended as new columns. itemLoadingsUsedMatrix of loadings used for projection. itemPCAObjPCA object used (if InputType is "PCAObj"), otherwise NULL. itemVarsUsedVariables used from Data for projection. itemCenterLogical flag indicating whether centering was applied for projection. itemScaleLogical flag indicating whether scaling was applied for projection.

**See also:** None documented.

## `ProjectRCI`

**Purpose:** Project a trained RCI object onto new data

**Canonical usage**
```r
ProjectRCI(
  data,
  Object,
  id_var = NULL
)
```

**Description:** Apply previously learned regression-based Reliable Change Index (RCI) models to a new or expanded longitudinal dataset WITHOUT refitting models.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `Object`: A SciDataReportR_RCI object created with codeCreateRCIObject().
- `id_var`: Optional ID column override.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `ID`: strongDeprecated (since 19.15.0). Use codeid_var instead.

**Returns:** A projected SciDataReportR_RCI object.

**See also:** None documented.

## `ProjectSOMCluster`

**Purpose:** Project cases through a fitted clustering model

**Canonical usage**
```r
ProjectCluster(object, new_df, ...)

Project_SOMClust(...)

ProjectSOMCluster(...)
```

**Description:** Projects codenew_df through the frozen preprocessing, reduction, and clustering layers stored in codeobject. Dispatch is determined by the fitted model's verbPipeline_* class; callers do not select a method-specific projector. Compatibility wrapper for codelink[=ProjectCluster]ProjectCluster(). Deprecated alias for codelink[=ProjectCluster]ProjectCluster().

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `ProjectCluster`, `Project_SOMClust`

**Arguments**
- `object`: A finalized object returned by verbCreateClusterModel_*().
- `new_df`: New data to project into the frozen cluster structure.
- `...`: Arguments passed to codelink[=ProjectCluster]ProjectCluster().

**Returns:** A method-specific projection result with the common codeProjectionFit contract and codeProbFit assignment table.

**See also:** None documented.

## `ProjectZScore`

**Purpose:** Project standardized scores onto new data using external parameters

**Canonical usage**
```r
ProjectZScore(
  data,
  variables = NULL,
  parameters,
  ParameterInputType = c("df_parameter", "ZScoreObj", "ExternalDataframe"),
  names_prefix = "Z_",
  RetainLabels = TRUE,
  RenameLabels = TRUE,
  center = TRUE,
  scale = TRUE
)

Project_ZScore(...)
```

**Description:** codeProject_ZScore() has been superseded by codeProjectZScore(). It remains available as a backwards-compatible alias and returns the same projected Z-score object.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `Project_ZScore`

**Arguments**
- `data`: Data frame on which to project scores.
- `variables`: Character vector; if NULL, project onto all variables for which parameters exist and that are present in df.
- `parameters`: Source of parameters, interpreted by ParameterInputType: itemize item "df_parameter": a data frame with cols Variable, N, Mean, SD item "ZScoreObj": output object from CreateZScoreObject() item "ExternalDataframe": a raw reference data frame; parameters are estimated via CreateZScoreObject() on that frame.
- `ParameterInputType`: One of "df_parameter", "ZScoreObj", "ExternalDataframe".
- `names_prefix`: Prefix for projected variable names.
- `RetainLabels`: Logical; if TRUE and Hmisc available, copy labels from df to new variables.
- `RenameLabels`: Logical; if TRUE, prefix labels the same way as names.
- `center`: Logical; used when ParameterInputType is "df_parameter" or "ExternalDataframe". Ignored for "ZScoreObj" (it uses stored flags).
- `scale`: Logical; same logic as codecenter.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.
- `...`: Arguments passed to codelink[=ProjectZScore]ProjectZScore().

**Returns:** List with same structure as CreateZScoreObject(): itemize item ZScores: projected standardized variables only item DataWithZ: df with projected scores appended item Parameters: parameter data frame actually used for projection item Center, Scale: flags used

**See also:** None documented.

## `ReadSciData`

**Purpose:** Read a scientific data file with optional inspection

**Canonical usage**
```r
ReadSciData(
  path,
  sheet = NULL,
  header_row = NULL,
  col_names = TRUE,
  range = NULL,
  inspect = TRUE,
  print_inspection = interactive(),
  strict = FALSE,
  guess_max = 10000,
  delim = NULL,
  repair_names = TRUE,
  inspect_styles = FALSE,
  fast_delimited = TRUE,
  ...
)
```

**Description:** codeReadSciData() imports common scientific data file formats while preserving original column names and labels as much as possible. It optionally calls codeInspectFile() before import to flag common issues such as multiple sheets, metadata rows, unnamed columns, duplicate column names, and formatting.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `path`: Path to the file.
- `sheet`: Sheet name or index for Excel files.
- `header_row`: Row containing column names for Excel files. If codeNULL, codeReadSciData() may use the probable header row detected by codeInspectFile().
- `col_names`: Passed to file readers where applicable. For Excel files, use codeTRUE, codeFALSE, or a character vector.
- `range`: Optional Excel range passed to codereadxl::read_excel().
- `inspect`: Logical. If codeTRUE, inspect the file before importing.
- `print_inspection`: Logical. If codeTRUE, print compact inspection results when issues are detected. Defaults to codeinteractive().
- `strict`: Logical. If codeTRUE, stop when inspection detects potential issues.
- `guess_max`: Maximum rows used for type guessing where supported.
- `delim`: Delimiter for code.txt files. Defaults to tab.
- `repair_names`: Logical. If codeTRUE, repair blank and duplicate column names after import.
- `inspect_styles`: Logical. If codeTRUE, inspect Excel workbook formatting. This can be slower and noisier for large workbooks.
- `fast_delimited`: Logical. If codeTRUE, use codedata.table::fread() for delimited text files when available and when no extra reader arguments are supplied through code.... This is usually much faster than codereadr for large code.csv, code.tsv, and code.txt files.
- `...`: Additional arguments passed to the underlying reader.

**Returns:** Imported data object. Inspection metadata is attached as the codescidata_inspection attribute when codeinspect = TRUE.

**See also:** codelink[=InspectFile]InspectFile() to inspect a file without importing it, codelink[=CreateVariableTypesTemplate]CreateVariableTypesTemplate() to build a codebook from the imported frame, and codelink[=RevalueData]RevalueData() to apply it.

## `removeString`

**Purpose:** Remove Strings from a Vector

**Canonical usage**
```r
removeString(Orig, Remove)
```

**Description:** This function removes strings from a vector that are present in another vector.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `Orig`: The original vector containing strings.
- `Remove`: The vector of strings to be removed from the original vector.

**Returns:** A vector containing the strings from the original vector that were not present in the removal vector.

**See also:** codelink[=getNumVars]getNumVars(), codelink[=getCatVars]getCatVars(), and codelink[=getBinaryVars]getBinaryVars(), which produce the vectors this usually trims.

## `ReplaceMissingCode`

**Purpose:** Replace Missing Codes with NA

**Canonical usage**
```r
ReplaceMissingCode(
  data,
  codebook
)
```

**Description:** This function replaces specified missing codes in a data frame with codeNA values based on a given variable codebook.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing the data.
- `codebook`: A data frame containing the variable codebook. It must have columns codeVariable and codeMissingCode. codeVariable names the column in codedata; codeMissingCode holds the code, or several codes separated by commas or semicolons, to replace with codeNA. Rows with a missing or blank codeMissingCode are skipped.
- `DataFrame`: strongDeprecated (since 19.15.0). Use codedata instead.
- `VariableCodebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A data frame with specified missing codes replaced by codeNA.

**See also:** codelink[=CreateVariableTypesTemplate]CreateVariableTypesTemplate(), which generates a codebook with a codeMissingCode column ready to fill in, and codelink[=RevalueData]RevalueData(), which applies missing codes as part of the full relabelling workflow.

## `ReplaceMissingLabels`

**Purpose:** Replace Missing Labels in Dataframe Columns

**Canonical usage**
```r

```

**Description:** This function iterates through the columns of a dataframe and assigns the column name as the label to any column that does not have a label.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** The input dataframe with missing labels replaced.

**See also:** None documented.

## `RevalueData`

**Purpose:** Revalue Data

**Canonical usage**
```r
RevalueData(
  data,
  codebook,
  missingVal = -999,
  splitchar = ";",
  on_error = c("stop", "warn")
)
```

**Description:** Revalues variables in a dataset using a VarTypes codebook.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A data.frame or tibble to be revalued.
- `codebook`: A data.frame with columns: Variable, Recode, Code, Type, Label, MissingCode. Only Variable is required. (Backward compatible: if MissingCode is absent/NA, will fall back to Missing.)
- `missingVal`: Default value to treat as missing when VarTypes$MissingCode is absent or NA.
- `splitchar`: Separator used in VarTypes$Code between pairs (default ";").
- `on_error`: Whether to stop at the first variable-level error (the default) or continue and record errors in the returned object.
- `DatatoRevalue`: strongDeprecated (since 19.15.0). Use codedata instead.
- `VarTypes`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A list with: RevaluedData (data), warninglist (character), recodedvars (character), not_in_data (character), and errors (data frame with codeVariable and codeError columns). In the default codeon_error = "stop" mode, an error names the offending variable and preserves the underlying message.

**See also:** None documented.

## `ReValueFactors`

**Purpose:** Revalue Factors

**Canonical usage**
```r
ReValueFactors(
  data,
  codebook
)
```

**Description:** This function revalues factor variables in a dataset according to the specifications provided in a codebook.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: The dataset to be revalued.
- `codebook`: A data frame containing information about the variables and how they should be revalued. It should have columns: Variable (variable names), Recode (yes/no for recoding), and Code (the revalue codes separated by "=" and ","). DO NOT USE COMMAS ANYWHERE ELSE IN THIS COLUMN.
- `DatatoRevalue`: strongDeprecated (since 19.15.0). Use codedata instead.
- `VarTypes`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** The revalued dataset.

**See also:** None documented.

## `reverseFactorLevels`

**Purpose:** Reverse Levels of Categorical Factors

**Canonical usage**
```r

```

**Description:** This function reverses the levels of specified categorical variables in a given dataframe.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe containing the categorical variables to be reversed.
- `variables`: A character vector of column names in the dataframe to reverse levels. These columns must be categorical factors.
- `df`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A dataframe with the levels of the specified factors reversed. Columns not specified in codevariables remain unchanged.

**See also:** None documented.

## `safe_merge`

**Purpose:** Safely merge two data frames with relationship-aware validation

**Canonical usage**
```r
safe_merge(
  df_before,
  df_add,
  by,
  name,
  method = c("exact", "closest_time"),
  time_var_before = NULL,
  time_var_add = NULL,
  min_match_rate = 0.95,
  harmonize_keys = TRUE,
  key_parser = NULL,
  stop_on_failed_numeric = TRUE,
  expected_relationship = c("one-to-one", "many-to-one", "one-to-many", "many-to-many",
    "auto"),
  fail_on_new_duplicate_variables = TRUE,
  fail_on_inherited_duplicate_variables = FALSE,
  ...
)
```

**Description:** codesafe_merge() performs a merge, validates its structure, logs merge metrics, and returns the merged data plus validation results. It is designed for reproducible database construction pipelines where row count, column count, key coverage, duplicate keys, expected merge relationships, and unresolved duplicate variables need to be audited every time the merge is run.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `df_before`: Data frame on the left side of the merge.
- `df_add`: Data frame to add to codedf_before.
- `by`: Character vector of merge key column names.
- `name`: Single character label used in the merge log and summary table.
- `method`: Merge method. code"exact" uses codedplyr::left_join(). code"closest_time" uses codeMerge_ByClosestTime().
- `time_var_before`: Required when codemethod = "closest_time". Time variable in codedf_before.
- `time_var_add`: Required when codemethod = "closest_time". Time variable in codedf_add.
- `min_match_rate`: Numeric between 0 and 1. Merges below this left-key match rate are marked as code"WARNING" unless another structural blocker makes them code"FAIL".
- `harmonize_keys`: Logical. If codeTRUE, keys are harmonized using codeHarmonizeMergeKeys() before merging.
- `key_parser`: Optional parser passed to codeHarmonizeMergeKeys().
- `stop_on_failed_numeric`: Logical passed to codeHarmonizeMergeKeys().
- `expected_relationship`: Character. Expected relationship between codedf_before and codedf_add under codeby. Defaults to code"one-to-one" to preserve strict historical behavior. One of: itemize item code"one-to-one": both sides should be unique by codeby. item code"many-to-one": left side may repeat keys, right side should be unique. item code"one-to-many": left side should be unique, right side may repeat. item code"many-to-many": both sides may repeat. item code"auto": infer and report relationship, but do not enforce it.
- `fail_on_new_duplicate_variables`: Logical. If codeTRUE, unresolved duplicate variable pairs introduced by the current merge cause the merge status to be code"FAIL". Default is codeTRUE.
- `fail_on_inherited_duplicate_variables`: Logical. If codeTRUE, unresolved duplicate variable pairs already present in codedf_before cause the merge status to be code"FAIL". Default is codeFALSE.
- `...`: Additional arguments passed to codeMerge_ByClosestTime() when codemethod = "closest_time".

**Returns:** A list with: itemize item codedata: Merged data frame. item codevalidation: Full validation object from codeValidateMerge(), with additional current-merge-aware duplicate-variable diagnostics. item codelog: One-row tibble containing merge log metrics. item codesummary: A codeknitr::kable() summary table.

**See also:** None documented.

## `scale_color_pvalue`

**Purpose:** Apply an evidence-aware p-value color scale

**Canonical usage**
```r
scale_color_pvalue(
  palette = "inferno",
  direction = -1,
  name = "P-value",
  breaks = c(1, 0.1, 0.05, 0.01, 0.001, 1e-05, 1e-08),
  labels = .format_pvalue_labels,
  limits = c(1e-08, 1),
  na.value = "grey80",
  guide = NULL,
  ...
)
```

**Description:** Maps raw p-values to either the Inferno or Viridis color palette using a threshold-aware transformation. The transformation allocates additional visual resolution around commonly interpreted p-value thresholds of 0.05, 0.01, and 0.001. Values above 0.05 are progressively desaturated to reduce their visual emphasis while retaining a continuous representation of the underlying p-values. The colorbar uses the same warped coordinates, so its raw p-value labels remain separated around these thresholds. Use this scale when p-values are mapped to the codecolor aesthetic.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `palette`: A character value specifying the color palette. Must be code"inferno" or code"viridis". Defaults to code"inferno".
- `direction`: A numeric value controlling the direction of the palette. Use code-1, the default, for darker colors to indicate smaller p-values and stronger statistical evidence. Use code1 to reverse the palette.
- `name`: The legend title. Defaults to code"P-value".
- `breaks`: A numeric vector of raw p-values to display as legend breaks. Defaults to commonly interpreted statistical thresholds and reference values.
- `labels`: A function or character vector used to label the legend breaks. By default, the legend displays raw p-values rather than transformed values.
- `limits`: A numeric vector containing the minimum and maximum p-values represented by the scale. Values outside these limits are squished to the nearest limit. Defaults to codec(1e-8, 1).
- `na.value`: The color assigned to missing p-values. Defaults to code"grey80".
- `guide`: A guide function or guide name. The default, codeNULL, uses a taller colorbar so the threshold-aware breaks remain legible. Supply code"colourbar" for ggplot2's standard-sized colorbar or another guide to override this behavior.
- `...`: Additional arguments passed to codelink[ggplot2:continuous_scale]ggplot2::continuous_scale().

**Returns:** A continuous ggplot2 color scale.

**See also:** None documented.

## `scale_color_SciData`

**Purpose:** SciDataReportR discrete color scale

**Canonical usage**
```r
scale_color_SciData(...)
```

**Description:** Apply the SciDataReportR qualitative palette to the codecolor aesthetic in a ggplot. Use this scale for categorical groups represented by points, lines, outlines, or other color-mapped geometries.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `...`: Additional arguments passed to codeggplot2::scale_color_manual().

**Returns:** A ggplot2 discrete color scale.

**See also:** None documented.

## `scale_fill_pvalue`

**Purpose:** Apply an evidence-aware p-value fill scale

**Canonical usage**
```r
scale_fill_pvalue(
  palette = "inferno",
  direction = -1,
  name = "P-value",
  breaks = c(1, 0.1, 0.05, 0.01, 0.001, 1e-05, 1e-08),
  labels = .format_pvalue_labels,
  limits = c(1e-08, 1),
  na.value = "grey80",
  guide = NULL,
  ...
)
```

**Description:** Maps raw p-values to either the Inferno or Viridis color palette using a threshold-aware transformation. The transformation allocates additional visual resolution around commonly interpreted p-value thresholds of 0.05, 0.01, and 0.001. Values above 0.05 are progressively desaturated to reduce their visual emphasis while retaining a continuous representation of the underlying p-values. The colorbar uses the same warped coordinates, so its raw p-value labels remain separated around these thresholds. Use this scale when p-values are mapped to the codefill aesthetic.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `palette`: A character value specifying the color palette. Must be code"inferno" or code"viridis". Defaults to code"inferno".
- `direction`: A numeric value controlling the direction of the palette. Use code-1, the default, for darker colors to indicate smaller p-values and stronger statistical evidence. Use code1 to reverse the palette.
- `name`: The legend title. Defaults to code"P-value".
- `breaks`: A numeric vector of raw p-values to display as legend breaks. Defaults to commonly interpreted statistical thresholds and reference values.
- `labels`: A function or character vector used to label the legend breaks. By default, the legend displays raw p-values rather than transformed values.
- `limits`: A numeric vector containing the minimum and maximum p-values represented by the scale. Values outside these limits are squished to the nearest limit. Defaults to codec(1e-8, 1).
- `na.value`: The fill color assigned to missing p-values. Defaults to code"grey80".
- `guide`: A guide function or guide name. The default, codeNULL, uses a taller colorbar so the threshold-aware breaks remain legible. Supply code"colourbar" for ggplot2's standard-sized colorbar or another guide to override this behavior.
- `...`: Additional arguments passed to codelink[ggplot2:continuous_scale]ggplot2::continuous_scale().

**Returns:** A continuous ggplot2 fill scale.

**See also:** None documented.

## `scale_fill_SciData`

**Purpose:** SciDataReportR discrete fill scale

**Canonical usage**
```r
scale_fill_SciData(...)
```

**Description:** Apply the SciDataReportR qualitative palette to the codefill aesthetic in a ggplot. Use this scale for categorical groups represented by bars, boxes, violins, areas, or other filled geometries.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `...`: Additional arguments passed to codeggplot2::scale_fill_manual().

**Returns:** A ggplot2 discrete fill scale.

**See also:** None documented.

## `SciDataPalette`

**Purpose:** SciDataReportR qualitative color palette

**Canonical usage**
```r
SciDataPalette(n = NULL, names = TRUE)
```

**Description:** Return colors from the SciDataReportR qualitative color system. The palette is designed for categorical scientific graphics and begins with the package anchor colors, Navy and Orange.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `n`: Number of colors to return. If codeNULL, return the full palette.
- `names`: Logical indicating whether to preserve the color names.

**Returns:** A character vector of hexadecimal color values.

**See also:** None documented.

## `ScreenBiomarkerPerformance`

**Purpose:** Screen biomarker performance across outcomes

**Canonical usage**
```r
ScreenBiomarkerPerformance(
  data,
  outcome_vars,
  biomarker_vars,
  covariates = NULL,
  PositiveLevel = NULL,
  OutcomeType = c("auto", "binary", "continuous"),
  Validation = c("none", "bootstrap", "cross_validation"),
  BootstrapR = 500,
  CVFolds = 10,
  CIBootstrapR = 200,
  CILevel = 0.95,
  HeatmapMetric = "AdjustedAUC",
  Seed = 123,
  Relabel = TRUE,
  codebook = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
```

**Description:** Applies codelink[=EvaluateBiomarkerPerformance]EvaluateBiomarkerPerformance() across many candidate biomarkers and one or more binary or continuous outcomes. Each biomarker-outcome pair uses all complete observations available for that pair and the requested covariates. The function returns comparison tables and screening plots, including an interactive-ready heatmap whose cells contain hover text.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame.
- `outcome_vars`: Character vector of outcome variable names.
- `biomarker_vars`: Character vector of biomarker variable names.
- `covariates`: Optional character vector of covariate variable names.
- `PositiveLevel`: Optional positive level. Supply one value for all binary outcomes or a named character vector with names matching codeoutcome_vars.
- `OutcomeType`: One of code"auto", code"binary", or code"continuous", applied to all outcomes unless auto-detection is used.
- `Validation`: One of code"none", code"bootstrap", or code"cross_validation".
- `BootstrapR`: Number of bootstrap resamples for internal validation.
- `CVFolds`: Number of cross-validation folds.
- `CIBootstrapR`: Number of bootstrap resamples for performance confidence intervals.
- `CILevel`: Confidence level. Default is code0.95.
- `HeatmapMetric`: Performance metric shown by the heatmap fill. Default is code"AdjustedAUC". Common alternatives are code"AUC", code"DeltaAUC", code"AdjustedR2", and code"DeltaR2".
- `Seed`: Random seed used for resampling.
- `Relabel`: Logical indicating whether labels should be used when available.
- `codebook`: Optional data frame with columns codeVariable and codeLabel.

**Returns:** A named list with codePerformanceTable, codeRegressionTable, codeThresholdTable, codeFailureTable, codeEvaluations, codePlots, and codeMetadata. For binary outcomes, codePlots also includes codeBiomarkerPanels (raw outcome-stratified distributions for continuous biomarkers) and codeROCFacets (adjusted ROC curves), each annotated with pair-specific performance.

**See also:** None documented.

## `SummarizeTransitions`

**Purpose:** Summarize participant transitions for a binary longitudinal condition

**Canonical usage**
```r
SummarizeTransitions(
  data,
  id_var,
  time_var,
  status_var,
  date_var = NULL,
  participant_subset = NULL,
  max_participants = NULL,
  order_participants_by = c("first_positive", "first_transition", "ever_positive",
    "ever_positive_then_burden", "input_order", "n_visits", "n_positive", "pct_positive"),
  x_axis_type = c("visit", "date", "time_from_baseline"),
  time_from_baseline_unit = c("days", "months", "years")
)
```

**Description:** Create participant-level and condition-level summary tables for a binary condition observed across repeated visits. This function uses the same transition logic and data preparation workflow as codePlotSwimmerTransitions() so that plotting and summary outputs remain aligned.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `data`: A data frame containing repeated observations per participant.
- `id_var`: Unquoted column name identifying the participant.
- `time_var`: Unquoted column name representing visit order, visit number, or time index.
- `status_var`: Unquoted column name representing the binary condition status.
- `date_var`: Optional unquoted visit date column. This is required when codex_axis_type = "date" or codex_axis_type = "time_from_baseline".
- `participant_subset`: Optional vector of participant IDs to include.
- `max_participants`: Optional maximum number of participants to retain after ordering is applied.
- `order_participants_by`: Character string controlling participant order. Options are code"first_positive", code"first_transition", code"ever_positive", code"ever_positive_then_burden", code"input_order", code"n_visits", code"n_positive", and code"pct_positive".
- `x_axis_type`: Character string indicating whether longitudinal ordering should follow aligned visit number (code"visit"), actual date (code"date"), or elapsed time from baseline (code"time_from_baseline").
- `time_from_baseline_unit`: Character string specifying the unit for codex_axis_type = "time_from_baseline". Options are code"days", code"months", and code"years".

**Returns:** A list with: itemize item codeparticipant_summary: participant-level summary table item codecondition_summary: one-row tibble with overall counts item codePlots: figures for the two summaries, described below

**See also:** codelink[=PlotSwimmerTransitions]PlotSwimmerTransitions() for the participant-level swimmer plot.

## `UnivariateRegressionTable`

**Purpose:** Univariate Regression Table

**Canonical usage**
```r
MakeUnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical"
)

UnivariateRegressionTable(
  data,
  outcome_vars,
  predictor_vars,
  covariates = NULL,
  Standardize = FALSE,
  Method = c("auto", "lm", "logistic"),
  LogisticExponentiate = TRUE,
  ReturnModels = FALSE
)
```

**Description:** Creates a list of univariate regression tables with variable labels and standardized coefficients (if specified). codeUnivariateRegressionTable() was renamed to codeMakeUnivariateRegressionTable() in SciDataReportR 20.5.0 to match the package's verbMake* naming convention. It remains available as a backwards-compatible synonym.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** `MakeUnivariateRegressionTable`

**Arguments**
- `data`: Dataframe containing the variables
- `outcome_vars`: Character vector of outcome variable names
- `predictor_vars`: Character vector of predictor variable names
- `covariates`: Character vector of covariate variable names (default: NULL)
- `Standardize`: Logical indicating whether to standardize numeric variables (default: FALSE)
- `Method`: Character. Regression method to use. code"auto" detects linear regression for numeric outcomes and logistic regression for two-level outcomes. code"lm" and code"logistic" force one model family for all outcomes.
- `LogisticExponentiate`: Logical. If codeTRUE, logistic regression estimates are exponentiated and reported as odds ratios.
- `ReturnModels`: Logical. If codeTRUE, return fitted model objects in codeModelSummaries. Default is codeFALSE to keep large screening runs lighter.
- `Relabel`: Logical; if TRUE (default), display attached variable labels.
- `TreatOrdinalAs`: How ordinal outcomes and predictors are handled.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.
- `OutcomeVars`: strongDeprecated (since 19.15.0). Use codeoutcome_vars instead.
- `PredictorVars`: strongDeprecated (since 19.15.0). Use codepredictor_vars instead.
- `Covars`: strongDeprecated (since 19.15.0). Use codecovariates instead.

**Returns:** A list containing: itemize item FormattedTable: A codegt table with formatted regression results item LargeTable: A codegt table with unformatted regression results item Results: A tidy dataframe with one row per estimated term. Columns: codeOutcome, codeOutcomeLabel, codeOutcomeFamily, codeEffectType, codePredictor, codePredictorLabel, codeTerm, codeLevel, codeTermLabel, codeN, codeEstimate, codeStdError, codeConfLow, codeConfHigh, codePValue, codeSignificant, and codeReferenceValue. This dataframe can be filtered and passed directly to codelink[=PlotForestFromTable]PlotForestFromTable(). item ModelSummaries: A list of fitted model objects when codeReturnModels = TRUE, otherwise codeNULL item Metadata: Outcome families and analysis settings

**See also:** codelink[=PlotForestFromTable]PlotForestFromTable() to visualize codeResults, codelink[=MultivariableRegressionTable]MultivariableRegressionTable() for mutually adjusted models, and codelink[=ApplyFDRCorrection]ApplyFDRCorrection() for multiple-comparison correction.

## `UpdateCodebook`

**Purpose:** Update an existing codebook based on a given dataframe

**Canonical usage**
```r
UpdateCodebook(
  data,
  codebook,
  RemoveMissing = TRUE,
  ReplaceLabels = FALSE
)
```

**Description:** This function updates a codebook by adding missing variables, removing variables that no longer exist in the dataframe (if specified), and optionally replacing outdated labels.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A dataframe for which the codebook needs to be updated.
- `codebook`: A dataframe representing the existing codebook with at least a 'Variable' column.
- `RemoveMissing`: Logical; if TRUE, removes variables from the codebook that are not in the dataframe.
- `ReplaceLabels`: Logical; if TRUE, replaces outdated labels in the codebook with new ones from the dataframe.
- `Dataframe`: strongDeprecated (since 19.15.0). Use codedata instead.
- `Codebook`: strongDeprecated (since 19.15.0). Use codecodebook instead.

**Returns:** A list containing: itemize item codeUpdatedCodebook: The updated codebook dataframe. item codeNewVariables: Variables present in the dataframe but missing from the original codebook. item codeNotExistingVariables: Variables in the codebook that are not present in the dataframe. item codeMismatchedLabels: A dataframe of variables with mismatched labels between the codebook and the dataframe.

**See also:** None documented.

## `UpdateDataDictionary`

**Purpose:** Update an existing data dictionary with new variables and types

**Canonical usage**
```r
UpdateDataDictionary(OldDataDictionary, NewDataFrame)
```

**Description:** This function updates an existing data dictionary with new variables and their corresponding types.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `OldDataDictionary`: The existing data dictionary to be updated.
- `NewDataFrame`: The new data frame containing additional variables.

**Returns:** A list containing the new variables and the updated data dictionary.

**See also:** None documented.

## `use_EDATemplate`

**Purpose:** Use the EDATemplate Quarto Template

**Canonical usage**
```r
use_EDATemplate(filename = "Reports/EDA_Report.qmd")
```

**Description:** Copies the EDATemplate template to the working directory.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `filename`: The name to save the Quarto file as (default: "Reports/EDA_Report.qmd").

**Returns:** Invisibly returns codefilename after copying the template. Called for its side effect of creating the report file.

**See also:** None documented.

## `use_GetStartedScript`

**Purpose:** Use the Get Started Script Template

**Canonical usage**
```r
use_GetStartedScript(filename = "Scripts/script_GetStarted.R")
```

**Description:** Copies the get-started R script template to the working directory.

**Deprecation status:** Current documented interface.

**Related exported aliases:** None.

**Arguments**
- `filename`: The name to save the R script as (default: "Scripts/script_GetStarted.R").

**Returns:** Invisibly returns codefilename after copying the template. Called for its side effect of creating the script file.

**See also:** None documented.

## `ValidateMerge`

**Purpose:** Validate a merge between two source data frames and a merged result

**Canonical usage**
```r
ValidateMerge(
  LeftData,
  RightData,
  MergedData,
  keys,
  expected_relationship = c("one-to-one", "many-to-one", "one-to-many", "many-to-many",
    "auto")
)
```

**Description:** codeValidateMerge() audits a completed merge by comparing the left/source data, right/add-on data, and merged result. It checks key coverage, key uniqueness, expected merge relationship, row inflation, overlapping non-key variables, unresolved code.x/code.y or verb_x/verb_y columns, and duplicated-variable conflicts.

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `LeftData`: A data frame used as the left side of the merge.
- `RightData`: A data frame used as the right side of the merge.
- `MergedData`: A data frame produced by merging codeLeftData and codeRightData.
- `keys`: Character vector of merge key column names.
- `Keys`: Deprecated. Use codekeys.
- `expected_relationship`: Character. Expected relationship between codeLeftData and codeRightData under codekeys. Defaults to code"one-to-one" to preserve strict historical behavior. One of: itemize item code"one-to-one": left keys and right keys should both be unique. item code"many-to-one": left keys may repeat, right keys should be unique. This is common when merging participant-level data into a longitudinal participant-visit master table by coderecord_id. item code"one-to-many": left keys should be unique, right keys may repeat. item code"many-to-many": both sides may repeat. This is allowed only when explicitly expected. item code"auto": infer and report relationship, but do not fail solely because of duplicate keys or relationship type.

**Returns:** A list containing merge validation summaries, checks, relationship audits, duplicate-key audits, coverage tables, overlap audits, duplicated variable audits, conflict tables, and suggested actions.

**See also:** None documented.

## `windsorize`

**Purpose:** Winsorize a numeric vector using SD or IQR thresholds

**Canonical usage**
```r
windsorize(
  data,
  sdlim = 2.5,
  iqrlim = 1.5,
  method = "sd",
  side = "both"
)
```

**Description:** This function performs winsorization on a numeric vector by capping extreme values at calculated lower and upper thresholds. Thresholds can be based on either standard deviation (assuming approximate normality) or interquartile range (robust to skewed distributions).

**Deprecation status:** Contains deprecated compatibility interface(s); use the current usage below.

**Related exported aliases:** None.

**Arguments**
- `data`: A numeric vector to be winsorized.
- `sdlim`: Numeric. Number of standard deviations for the "sd" method.
- `iqrlim`: Numeric. Multiplier for the IQR when method = "iqr" (default 1.5).
- `method`: Character string specifying the method: "sd" (default) or "iqr".
- `side`: Character string specifying which tail(s) to winsorize: "both" (default), "right" for high values only, or "left" for low values only.
- `Data`: strongDeprecated (since 19.15.0). Use codedata instead.

**Returns:** A numeric vector with values winsorized to the specified thresholds.

**See also:** None documented.

