# SciDataReportR 21.0.0

# SciDataReportR 20.24.0

## Clustering API is standardized around fitted pipelines and one projector

* **Breaking unreleased clustering-output cleanup.** Finalized fits and
  projections now return `DataWithClusters` (formerly `df_with_clusters`) and
  `ClusterVariableName` (formerly `ClusterName`). The mistyped
  `ModelInfo_Mclust`, public `complete_rows`, `CandidateAudit`, and injected
  `.scidr_rowid` columns have been removed. `ModelInfo_MClust` is the canonical
  Mclust-specific layer; `ModelInfo` remains the generic compatibility alias.

* SOM + Mclust now reports `MinProfileNodeN`/`MaxProfileNodeN` alongside their
  explicit proportions, and `BLRTStatistic`/`BLRTPValue`. The candidate plot
  includes an unselected BLRT p-value panel with a neutral 0.05 reference line.

* Clustering constructors now use `CreateClusterModel_<reduction>_<method>()`,
  for example `CreateClusterModel_PCA_MClust()` and
  `CreateClusterModel_SOM_HDBSCAN()`. Fitted objects use matching
  `Pipeline_<reduction>_<method>` classes.

* `ProjectCluster(object, new_df)` is the single public projection interface.
  It dispatches from the fitted pipeline class and retains each method's frozen
  preprocessing, reduction, support checks, and projection diagnostics.

* `CreateSOMClusterModel()`, `ProjectSOMCluster()`, `Pipeline_SOMClust()`, and
  `Project_SOMClust()` remain for established SOM compatibility. The first two
  emit lifecycle deprecation warnings and delegate to the canonical API.

## The SciData palette is now the default categorical color scheme

* Every plot the package produces now draws categorical color from
  `SciDataPalette()`. Previously each function chose its own: some fell back to
  ggplot2's default hue palette, others hardcoded a tableau-style vector, an
  `hcl.colors()` name, or a `paletteer` string. Clusters, groups, cohorts, and
  category levels now look like one system across the whole package.

* **Only categorical color changed.** Continuous and diverging scales
  (correlation and effect-size heatmaps, the p-value scales) are untouched, as
  are the meaning-carrying scales: PASS/WARNING/FAIL status, the signed
  significance ladders, RCI classification, volcano significance tiers, the
  grey used for missing data, and the grey reserved for the density-based
  `Noise` cluster.

* **Fixes a latent failure.** `CreatePCATable()` and `CreateZScorePlot()` shared
  a 20-color vector that went straight into `scale_color_manual()`, so both
  failed with `Insufficient values in manual scale` once a plot had more than 20
  variable categories. Package plots now extend the palette instead of erroring,
  however many categories a variable has.

* Labels drawn inside filled shapes now pick black or white by background
  luminance. `PlotCategoricalDistributions()` draws its labels inside the bars,
  and the palette's darker colors would have made them unreadable.

* Palette arguments keep their existing signatures. `PlotClusterBoxplot(Palette)`,
  `PlotSplitViolin(color_palette)`, `PlotSpiderChart(Palette)` and
  `Plot2GroupStats(palette)` all still accept exactly what they accepted before;
  only the default changed. `PlotSpiderChart()` and `Plot2GroupStats()` now
  default to `NULL`, meaning "use the package palette" - passing `"Dark 3"` or
  `"pals::alphabet"` explicitly works as it always did.

* `SciDataPalette()` itself is unchanged, including erroring when asked for more
  than its 34 colors.

## Cluster occupancy figures removed from the clustering pipelines

* The `occupancy` bar chart is no longer attached to `ModelInfo$plots` by
  `CreateClusterModel_MClust()`, `CreateClusterModel_KMeans()`,
  `CreateClusterModel_PCA_MClust()`, `CreateClusterModel_PCA_KMeans()`,
  `CreateClusterModel_HDBSCAN()`, `CreateClusterModel_Gower_PAM()`,
  `CreateClusterModel_LatentClass()`, `CreateClusterModel_MCA_MClust()`, or
  `CreateSOMClusterModel()`. Cluster sizes are already reported in the fit
  tables and in `ProbFit$individual`, so the figure added a panel without
  adding information.

* `PlotClusterOccupancy()` is unchanged and still exported, so the plot can
  still be drawn directly from a model's `DataWithClusters`.

* The training-versus-projected occupancy comparison in the `Project*`
  functions is unaffected: it compares two distributions rather than
  describing one, which is a projection-quality check rather than a
  restatement of the cluster sizes.

## Symmetric p-value matrices are corrected once per pair

* `ApplyFDRCorrection()` treated a symmetric matrix - the shape produced
  whenever a variable set is correlated against itself, as
  `PlotCorrelationsHeatmap()` does when `outcome_vars` is left `NULL` - as a
  family of `n * (n - 1)` tests, when only `n * (n - 1) / 2` were run. New
  `symmetric` argument, defaulting to `"auto"`, detects the case and corrects
  each pair once, mirroring the adjusted values back so the matrix stays
  symmetric. Pass `symmetric = FALSE` for the old behavior or
  `symmetric = TRUE` to require it.

* Benjamini-Hochberg is invariant to exact duplication, so `method = "fdr"`
  results - the default, and what every plotting function uses - are unchanged.
  `"bonferroni"`, `"holm"` and `"hochberg"` were coming out exactly twice as
  large as they should have been, and are now correct.

* The diagonal of a symmetric matrix holds self-comparisons that were never
  tested. It is now excluded from the family and returned as `NA`; set
  `include_diagonal = TRUE` to keep the old behavior. This was the damaging
  case: where the diagonal carried real values, a self-correlation p-value of 0
  entered the family as the most significant test in it and pulled every
  off-diagonal result down, making findings look stronger than they were.

## Misspelled variables are now an error instead of a silent result change

* The matrix and heatmap functions used to handle an unknown variable name five
  different ways. `PlotCorrelationsHeatmap()`, `PlotAnovaRelationshipsMatrix()`,
  `PlotChiSqCovar()`, `PlotMiningMatrix()`, `PlotNumInteractionEffectsMatrix()`
  and `PlotCatInteractionEffectsMatrix()` dropped the name silently, returning a
  smaller results matrix with no warning that anything had gone missing;
  `PlotPhiHeatmap()` and `PlotPointCorrelationsHeatmap()` failed with base R's
  `"undefined columns selected"`, which never named the offending variable; and
  `PlotDirectionalHeatmaps()` and `PlotSpiderChart()` each raised their own
  differently worded error. All of them now stop with one message that names both
  the variable and the argument it was supplied to, for example
  `Variables not found in `data`: NOPE (supplied to `predictor_vars`)`.

* **A misspelled covariate was the more serious case.**
  `PlotCorrelationsHeatmap()`, `PlotChiSqCovar()`,
  `PlotAnovaRelationshipsMatrix()`, `PlotMiningMatrix()` and
  `PlotNumInteractionEffectsMatrix()` silently ignored a covariate that did not
  match a column, so the function computed and reported **unadjusted** statistics
  while the caller believed they were adjusted. Unknown covariates now error.

* Validation happens against the data frame as supplied, before ordinal handling
  renames anything, so variables legitimately removed by `TreatOrdinalAs` are
  unaffected. Auto-detected variable sets (`variables = NULL`) still work as
  before, and `PlotCorrelationsHeatmap()` still reports covariates dropped by
  ordinal handling through `CovariatesMissing`.

* New internal helpers `ScidrValidateVariables()` and `ScidrValidateVariable()`
  provide the single message format, replacing the
  `unique(intersect(as.character(vars), names(Data)))` idiom that caused the
  silent dropping.

## Significance stars agree across the package

* Six separate star implementations disagreed on whether the cut points were
  inclusive: a p-value of exactly 0.001 earned `***` in
  `MultivariableRegressionTable()` and `PlotInteractionEffectsMatrix()` but only
  `**` in `PlotSplitViolin()`, while every plot caption stated the strict `<`
  form regardless. All of them now route through one `ScidrPValueStars()` using
  inclusive upper bounds (`p <= 0.05` is `*`), which matches
  `rstatix::add_significance()` — already used by `PlotPhiHeatmap()`,
  `PlotPointCorrelationsHeatmap()` and `PlotAnovaRelationshipsMatrix()` — and is
  verified against it at every boundary.

* Captions and legend labels now state `<=` to match the thresholds actually
  applied. This affects `geom_starcaption()`, `MakePairwiseHeatmap()` captions,
  `PlotMiningMatrix()` legend labels, and the `Plot2GroupStats()` shape legend
  (whose `cut(right = TRUE)` breaks were always inclusive despite reading `<`).

## Other fixes

* `PlotDirectionalHeatmaps()` had the deprecated `yVars` in formal position 3, so
  `PlotDirectionalHeatmaps(df, vars, TRUE)` bound `TRUE` to `yVars` rather than
  `Relabel`. All deprecated formals are now trailing, as elsewhere in the package.

* `PlotDirectionalHeatmaps()`'s `Ordinal` argument was documented as "reserved for
  future use" and as not affecting the computed tiles, but it is passed to
  `PlotPointCorrelationsHeatmap()` and does change how ordinal variables are
  treated in the binary~continuous block. The documentation now describes what it
  actually does.

# SciDataReportR 20.23.0

* Clustering pipelines no longer return a flat `plots` list. Figures are stored
  beside the object they describe, following the layout `CreateSOMClusterModel()`
  already used: `fit_plot` reviews candidates, `ModelInfo$plots` describes the
  selected solution, `ModelInfo$FitDiagnostics$plots` describes how individual
  training cases sit inside it, `ProbFit$plots` describes membership confidence,
  `Stability$plots` describes bootstrap reproducibility, and `ProjectionFit$plots`
  describes projected cases against the frozen training reference. The previous
  `plots` list mixed aliases, duplicated entries under two names, and stored
  interactive widgets alongside `ggplot` objects; it has been removed.

* Every method now carries figures appropriate to that method rather than a
  shared lowest common denominator. K-means and Gower/PAM gain per-participant
  silhouette profiles, elbow and average-silhouette curves; Mclust gains BIC,
  ICL, entropy, and classification-uncertainty maps; HDBSCAN gains density-grid
  review and per-cluster persistence; latent class analysis gains item response
  profiles; and the PCA/MCA pipelines lead with scree and loadings. All methods
  gain a frozen two-dimensional review map, cluster centre profiles, and a
  distance histogram/boxplot/ECDF triad against a frozen high-distance cutoff.

* `Project*Cluster()` results gain `ProjectionFit`, which triages every projected
  case against the frozen training reference into `Good fit`,
  `Uncertain membership`, `Poor fit to training structure`, or
  `Potential novel phenotype`, and adds `Projection_Fit_Class` to
  `DataWithClusters`. This extends the triage `ProjectSOMCluster()` already
  performed to all clustering methods.

* `CreateSOMClusterModel()` and `ProjectSOMCluster()` keep their existing figures
  unchanged. The SOM model gains `ModelInfo_SOM$plots$occupancy` and
  `Stability$plots`; the projection exposes its diagnostics as `ProjectionFit`
  (still also available as `SOMProj`) and gains `som_grid_map`.

* `CreateSOMClusterModel(method = "finalize", stability_resamples > 0)` now
  refits through the exploratory path by a direct recursive call rather than
  evaluating a reconstructed call in the caller's frame, which failed when
  arguments were supplied as local variables.

* Bootstrap stability resamples now preserve observed factor levels, and the MCA
  pipeline freezes its reduction on observed categories only. Previously every
  MCA stability replicate could fail when a resample dropped a rare category, and
  a declared-but-unobserved factor level made projection fail outright.
  `ProjectCluster()` now warns and treats categories absent from the
  training model as incomplete cases instead of erroring.

* New exported figure helpers: `PlotClusterOccupancy()`, `PlotClusterMap()`,
  `PlotClusterSilhouette()`, `PlotClusterCentreHeatmap()`,
  `PlotClusterCentreProfile()`, `PlotClusterComposition()`, and
  `PlotClusterDiagnostic()`. `PlotProjectionDiagnostics()` has been renamed to
  `PlotClusterDiagnostic()`.

# SciDataReportR 20.22.0

# SciDataReportR 20.21.0

# SciDataReportR 20.20.0

# SciDataReportR 20.19.0

# SciDataReportR 20.18.0

# SciDataReportR 20.17.0

# SciDataReportR 20.16.0

# SciDataReportR 20.15.0

* `MakeComparisonTable()` now reports absolute Cohen's d for two-group parametric comparisons, including covariate-adjusted d from estimated marginal means, and Cohen's f for multi-group omnibus comparisons. Effect-size captions identify the scales present, provide qualified magnitude guides, and warn when d and f should not be compared numerically.

# SciDataReportR 20.14.0

# SciDataReportR 20.13.0

# SciDataReportR 20.12.0

* `PlotVolcanoEffects()` now reports more detail in its `ResultsTable` and point tooltips. Continuous outcomes gain `R` (zero-order Pearson correlation) and `AdjustedR` (covariate-adjusted partial correlation). Two-group categorical outcomes gain `Group1Level`, `Group2Level`, `Group1Mean`, and `Group2Mean` (raw predictor means within each group). These values also appear in the `Tooltip` column so they show up when the plot is passed to `plotly::ggplotly(tooltip = "text")`.
* New `FreezeTableHeader()` wraps a `gtsummary` table (or data frame) so its header row stays frozen while scrolling long tables in HTML Quarto/R Markdown output.
* `MakeComparisonTable()` now uses safe internal names for its pairwise columns. Previously contrast labels such as `"3 - 1"` produced column names like `pw_X3...1`, whose trailing `...1` collided with tidyverse name-repair and broke downstream tibble round-trips (for example `gtsummary::as_kable_extra()`, and therefore `FreezeTableHeader()`). Displayed contrast headers are unchanged.
* `MakeUnivariateRegressionTable()` now extracts regression results directly from model coefficient tables and formats the final display with `gt`, avoiding the slow per-model `gtsummary::tbl_regression()` path. Fitted model objects are now skipped by default; set `ReturnModels = TRUE` to include them in `ModelSummaries`.
* `ReadSciData()` now uses `data.table::fread()` for ordinary delimited files when available, which substantially speeds up large `.csv`, `.tsv`, and `.txt` imports. Set `fast_delimited = FALSE` to force the previous `readr` path.

# SciDataReportR 20.11.0

# SciDataReportR 20.10.0

# SciDataReportR 20.9.0

# SciDataReportR 20.8.0

# SciDataReportR 20.7.0

# SciDataReportR 20.6.0

# SciDataReportR 20.5.0

* `UnivariateRegressionTable()` was renamed to `MakeUnivariateRegressionTable()` and `plotForestFromTable()` was renamed to `PlotForestFromTable()` to match the package's naming conventions. The old names remain available as backwards-compatible synonyms (soft-deprecated).
* `MakeUnivariateRegressionTable()` now returns a `Results` element: a tidy dataframe of estimates, confidence intervals, and p-values, one row per term.
* `PlotForestFromTable()` now plots from the `Results` dataframe, and also accepts a (filtered) `Results` dataframe directly. Objects created by older package versions still plot via a fallback extraction.

# SciDataReportR 20.4.0

# SciDataReportR 20.3.0

# SciDataReportR 20.2.0

# SciDataReportR 20.1.0

# SciDataReportR 20.0.0

# SciDataReportR 19.14.0

# SciDataReportR 19.13.0

# SciDataReportR 19.12.0

# SciDataReportR 19.11.0

# SciDataReportR 19.10.0

# SciDataReportR 19.9.0

# SciDataReportR 19.8.0

# SciDataReportR 19.7.0

# SciDataReportR 19.6.0

# SciDataReportR 19.5.0

# SciDataReportR 19.4.0

# SciDataReportR 19.3.0

# SciDataReportR 19.2.0

# SciDataReportR 19.1.0

# SciDataReportR 19.0.0

# SciDataReportR 18.0.0

# SciDataReportR 17.3.0

# SciDataReportR 17.2.0

# SciDataReportR 17.1.0

# SciDataReportR 17.0.0

## Dependency stability

- Internalized the minimal half-violin geom used by `PlotSplitViolin()` so the
  split-violin workflow no longer depends on the archived `gghalves` package.

# SciDataReportR 16.25.0

## Documentation and infrastructure positioning

- Repositioned the package as scientific workflow infrastructure for
  reproducible life science data reporting.
- Expanded the README with workflow families, a visualization gallery,
  implemented methods, function dependency chains, and future heatmap/volcano
  plot roadmap items.
- Reorganized pkgdown navigation and reference topics around workflows,
  visualization functions, projection chains, and reusable infrastructure.
- Clarified the downstream relationship between `PlotCorrelationsHeatmap()`,
  `add_r_and_stars()`, and `geom_starcaption()`.
- Added a Quarto getting-started article based on the R/Medicine workflow narrative.

## Workflow-oriented public names

- Added canonical workflow-oriented names for reusable objects, models,
  projections, plots, and metadata workflows.
- Added `CreatePCAObject()`, `CreateMCAObject()`, `PlotZScore()`,
  `PlotPathway_KT()`, `MakeDataDictionary()`, `ProjectZScore()`,
  `ProjectSOMCluster()`, `CreateZScoreObject()`, `CreateMScoreObject()`,
  `CreateNormativeTScoreModel()`, and `CreateSOMClusterModel()`.
- Preserved backward compatibility: the older public names remain exported as
  wrappers and do not emit lifecycle warnings.
- Documented `PlotPathway_KT()` as a SciDataReportR visualization that is
  expected to move to a future metabolomics-focused package.

# SciDataReportR 16.24.0

# SciDataReportR 16.23.0

# SciDataReportR 16.22.0

# SciDataReportR 16.21.0

# SciDataReportR 16.20.0

# SciDataReportR 16.19.0

# SciDataReportR 16.18.0

# SciDataReportR 16.17.0

# SciDataReportR 16.16.0

# SciDataReportR 16.15.0

# SciDataReportR 16.14.0

# SciDataReportR 16.13.0

# SciDataReportR 16.12.0

# SciDataReportR 16.11.0

# SciDataReportR 16.10.0

# SciDataReportR 16.9.0

# SciDataReportR 16.8.0

# SciDataReportR 16.7.0

# SciDataReportR 16.6.0

# SciDataReportR 16.5.0

# SciDataReportR 16.4.0

# SciDataReportR 16.3.0

# SciDataReportR 16.2.0

# SciDataReportR 16.1.0

# SciDataReportR 16.0.0

# SciDataReportR 15.12.0

# SciDataReportR 15.11.0

# SciDataReportR 15.10.0

# SciDataReportR 15.9.0

# SciDataReportR 15.8.0

# SciDataReportR 15.7.0

# SciDataReportR 15.6.0

# SciDataReportR 15.5.0

# SciDataReportR 15.4.0

# SciDataReportR 15.3.0

# SciDataReportR 15.2.0

# SciDataReportR 15.1.0

# SciDataReportR 15.0.0

# SciDataReportR 14.19.0

# SciDataReportR 14.18.0

# SciDataReportR 14.17.0

# SciDataReportR 14.16.0

# SciDataReportR 14.15.0

# SciDataReportR 14.14.0

# SciDataReportR 14.13.0

# SciDataReportR 14.12.0

# SciDataReportR 14.11.0

# SciDataReportR 14.10.0

# SciDataReportR 14.9.0

# SciDataReportR 14.8.0

# SciDataReportR 14.7.0

# SciDataReportR 14.6.0

# SciDataReportR 14.5.0

# SciDataReportR 14.4.0

# SciDataReportR 14.3.0

# SciDataReportR 14.2.0

# SciDataReportR 14.1.0

# SciDataReportR 14.0.0

# SciDataReportR 13.6.0

# SciDataReportR 13.5.0

# SciDataReportR 13.4.0

# SciDataReportR 13.3.0

# SciDataReportR 13.2.0

# SciDataReportR 13.1.0

# SciDataReportR 13.0.0

# SciDataReportR 12.22.0

# SciDataReportR 12.21.0

# SciDataReportR 12.20.0

# SciDataReportR 12.19.0

# SciDataReportR 12.18.0

# SciDataReportR 12.17.0

# SciDataReportR 12.16.0

# SciDataReportR 12.15.0

# SciDataReportR 12.14.0

# SciDataReportR 12.13.0

# SciDataReportR 12.12.0

# SciDataReportR 12.11.0

# SciDataReportR 12.10.0

# SciDataReportR 12.9.0

# SciDataReportR 12.8.0

# SciDataReportR 12.7.0

# SciDataReportR 12.6.0

# SciDataReportR 12.5.0

# SciDataReportR 12.4.0

# SciDataReportR 12.3.0

# SciDataReportR 12.2.0

# SciDataReportR 12.1.0

# SciDataReportR 12.0.0

# SciDataReportR 11.5.0

# SciDataReportR 11.4.0

# SciDataReportR 11.3.0

# SciDataReportR 11.2.0

# SciDataReportR 11.1.0

# SciDataReportR 11.0.0

# SciDataReportR 10.34.0

# SciDataReportR 10.33.0

# SciDataReportR 10.32.0

# SciDataReportR 10.31.0

# SciDataReportR 10.30.0

# SciDataReportR 10.29.0

# SciDataReportR 10.28.0

# SciDataReportR 10.27.0

# SciDataReportR 10.26.0

# SciDataReportR 10.25.0

# SciDataReportR 10.24.0

# SciDataReportR 10.23.0

# SciDataReportR 10.22.0

# SciDataReportR 10.21.0

# SciDataReportR 10.20.0

# SciDataReportR 10.19.0

# SciDataReportR 10.18.0

# SciDataReportR 10.17.0

# SciDataReportR 10.16.0

# SciDataReportR 10.15.0

# SciDataReportR 10.14.0

# SciDataReportR 10.13.0

# SciDataReportR 10.12.0

# SciDataReportR 10.11.0

# SciDataReportR 10.10.0

# SciDataReportR 10.9.0

# SciDataReportR 10.8.0

# SciDataReportR 10.7.0

# SciDataReportR 10.6.0

# SciDataReportR 10.5.0

# SciDataReportR 10.4.0

# SciDataReportR 10.3.0

# SciDataReportR 10.2.0

# SciDataReportR 10.1.0

# SciDataReportR 10.0.0

# SciDataReportR 9.17.0

# SciDataReportR 9.16.0

# SciDataReportR 9.15.0

# SciDataReportR 9.14.0

# SciDataReportR 9.13.0

# SciDataReportR 9.12.0

# SciDataReportR 9.11.0

# SciDataReportR 9.10.0

# SciDataReportR 9.9.0

# SciDataReportR 9.8.0

# SciDataReportR 9.7.0

# SciDataReportR 9.6.0

# SciDataReportR 9.5.0

# SciDataReportR 9.4.0

# SciDataReportR 9.3.0

# SciDataReportR 9.2.0

# SciDataReportR 9.1.0

# SciDataReportR 9.0.0

# SciDataReportR 8.23.0

# SciDataReportR 8.22.0

# SciDataReportR 8.21.0

# SciDataReportR 8.20.0

# SciDataReportR 8.19.0

# SciDataReportR 8.18.0

# SciDataReportR 8.17.0

# SciDataReportR 8.16.0

# SciDataReportR 8.15.0

# SciDataReportR 8.14.0

# SciDataReportR 8.13.0

# SciDataReportR 8.12.0

# SciDataReportR 8.11.0

# SciDataReportR 8.10.0

# SciDataReportR 8.9.0

# SciDataReportR 8.8.0

# SciDataReportR 8.7.0

# SciDataReportR 8.6.0

* Initial CRAN submission.
