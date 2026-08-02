# Changelog

## SciDataReportR 20.21.0

## SciDataReportR 20.20.0

## SciDataReportR 20.19.0

## SciDataReportR 20.18.0

## SciDataReportR 20.17.0

## SciDataReportR 20.16.0

## SciDataReportR 20.15.0

- [`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
  now reports absolute Cohen’s d for two-group parametric comparisons,
  including covariate-adjusted d from estimated marginal means, and
  Cohen’s f for multi-group omnibus comparisons. Effect-size captions
  identify the scales present, provide qualified magnitude guides, and
  warn when d and f should not be compared numerically.

## SciDataReportR 20.14.0

## SciDataReportR 20.13.0

## SciDataReportR 20.12.0

- [`PlotVolcanoEffects()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotVolcanoEffects.md)
  now reports more detail in its `ResultsTable` and point tooltips.
  Continuous outcomes gain `R` (zero-order Pearson correlation) and
  `AdjustedR` (covariate-adjusted partial correlation). Two-group
  categorical outcomes gain `Group1Level`, `Group2Level`, `Group1Mean`,
  and `Group2Mean` (raw predictor means within each group). These values
  also appear in the `Tooltip` column so they show up when the plot is
  passed to `plotly::ggplotly(tooltip = "text")`.
- New
  [`FreezeTableHeader()`](https://rdastgh1.github.io/SciDataReportR/reference/FreezeTableHeader.md)
  wraps a `gtsummary` table (or data frame) so its header row stays
  frozen while scrolling long tables in HTML Quarto/R Markdown output.
- [`MakeComparisonTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeComparisonTable.md)
  now uses safe internal names for its pairwise columns. Previously
  contrast labels such as `"3 - 1"` produced column names like
  `pw_X3...1`, whose trailing `...1` collided with tidyverse name-repair
  and broke downstream tibble round-trips (for example
  [`gtsummary::as_kable_extra()`](https://www.danieldsjoberg.com/gtsummary/reference/as_kable_extra.html),
  and therefore
  [`FreezeTableHeader()`](https://rdastgh1.github.io/SciDataReportR/reference/FreezeTableHeader.md)).
  Displayed contrast headers are unchanged.
- [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  now extracts regression results directly from model coefficient tables
  and formats the final display with `gt`, avoiding the slow per-model
  [`gtsummary::tbl_regression()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_regression.html)
  path. Fitted model objects are now skipped by default; set
  `ReturnModels = TRUE` to include them in `ModelSummaries`.
- [`ReadSciData()`](https://rdastgh1.github.io/SciDataReportR/reference/ReadSciData.md)
  now uses
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  for ordinary delimited files when available, which substantially
  speeds up large `.csv`, `.tsv`, and `.txt` imports. Set
  `fast_delimited = FALSE` to force the previous `readr` path.

## SciDataReportR 20.11.0

## SciDataReportR 20.10.0

## SciDataReportR 20.9.0

## SciDataReportR 20.8.0

## SciDataReportR 20.7.0

## SciDataReportR 20.6.0

## SciDataReportR 20.5.0

- [`UnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  was renamed to
  [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  and
  [`plotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md)
  was renamed to
  [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md)
  to match the package’s naming conventions. The old names remain
  available as backwards-compatible synonyms (soft-deprecated).
- [`MakeUnivariateRegressionTable()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeUnivariateRegressionTable.md)
  now returns a `Results` element: a tidy dataframe of estimates,
  confidence intervals, and p-values, one row per term.
- [`PlotForestFromTable()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotForestFromTable.md)
  now plots from the `Results` dataframe, and also accepts a (filtered)
  `Results` dataframe directly. Objects created by older package
  versions still plot via a fallback extraction.

## SciDataReportR 20.4.0

## SciDataReportR 20.3.0

## SciDataReportR 20.2.0

## SciDataReportR 20.1.0

## SciDataReportR 20.0.0

## SciDataReportR 19.14.0

## SciDataReportR 19.13.0

## SciDataReportR 19.12.0

## SciDataReportR 19.11.0

## SciDataReportR 19.10.0

## SciDataReportR 19.9.0

## SciDataReportR 19.8.0

## SciDataReportR 19.7.0

## SciDataReportR 19.6.0

## SciDataReportR 19.5.0

## SciDataReportR 19.4.0

## SciDataReportR 19.3.0

## SciDataReportR 19.2.0

## SciDataReportR 19.1.0

## SciDataReportR 19.0.0

## SciDataReportR 18.0.0

## SciDataReportR 17.3.0

## SciDataReportR 17.2.0

## SciDataReportR 17.1.0

## SciDataReportR 17.0.0

### Dependency stability

- Internalized the minimal half-violin geom used by
  [`PlotSplitViolin()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSplitViolin.md)
  so the split-violin workflow no longer depends on the archived
  `gghalves` package.

## SciDataReportR 16.25.0

### Documentation and infrastructure positioning

- Repositioned the package as scientific workflow infrastructure for
  reproducible life science data reporting.
- Expanded the README with workflow families, a visualization gallery,
  implemented methods, function dependency chains, and future
  heatmap/volcano plot roadmap items.
- Reorganized pkgdown navigation and reference topics around workflows,
  visualization functions, projection chains, and reusable
  infrastructure.
- Clarified the downstream relationship between
  [`PlotCorrelationsHeatmap()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotCorrelationsHeatmap.md),
  [`add_r_and_stars()`](https://rdastgh1.github.io/SciDataReportR/reference/add_r_and_stars.md),
  and
  [`geom_starcaption()`](https://rdastgh1.github.io/SciDataReportR/reference/geom_starcaption.md).
- Added a Quarto getting-started article based on the R/Medicine
  workflow narrative.

### Workflow-oriented public names

- Added canonical workflow-oriented names for reusable objects, models,
  projections, plots, and metadata workflows.
- Added
  [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md),
  [`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md),
  [`PlotZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotZScore.md),
  [`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md),
  [`MakeDataDictionary()`](https://rdastgh1.github.io/SciDataReportR/reference/MakeDataDictionary.md),
  [`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md),
  [`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md),
  [`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md),
  [`CreateMScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMScoreObject.md),
  [`CreateNormativeTScoreModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateNormativeTScoreModel.md),
  and
  [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md).
- Preserved backward compatibility: the older public names remain
  exported as wrappers and do not emit lifecycle warnings.
- Documented
  [`PlotPathway_KT()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotPathway_KT.md)
  as a SciDataReportR visualization that is expected to move to a future
  metabolomics-focused package.

## SciDataReportR 16.24.0

## SciDataReportR 16.23.0

## SciDataReportR 16.22.0

## SciDataReportR 16.21.0

## SciDataReportR 16.20.0

## SciDataReportR 16.19.0

## SciDataReportR 16.18.0

## SciDataReportR 16.17.0

## SciDataReportR 16.16.0

## SciDataReportR 16.15.0

## SciDataReportR 16.14.0

## SciDataReportR 16.13.0

## SciDataReportR 16.12.0

## SciDataReportR 16.11.0

## SciDataReportR 16.10.0

## SciDataReportR 16.9.0

## SciDataReportR 16.8.0

## SciDataReportR 16.7.0

## SciDataReportR 16.6.0

## SciDataReportR 16.5.0

## SciDataReportR 16.4.0

## SciDataReportR 16.3.0

## SciDataReportR 16.2.0

## SciDataReportR 16.1.0

## SciDataReportR 16.0.0

## SciDataReportR 15.12.0

## SciDataReportR 15.11.0

## SciDataReportR 15.10.0

## SciDataReportR 15.9.0

## SciDataReportR 15.8.0

## SciDataReportR 15.7.0

## SciDataReportR 15.6.0

## SciDataReportR 15.5.0

## SciDataReportR 15.4.0

## SciDataReportR 15.3.0

## SciDataReportR 15.2.0

## SciDataReportR 15.1.0

## SciDataReportR 15.0.0

## SciDataReportR 14.19.0

## SciDataReportR 14.18.0

## SciDataReportR 14.17.0

## SciDataReportR 14.16.0

## SciDataReportR 14.15.0

## SciDataReportR 14.14.0

## SciDataReportR 14.13.0

## SciDataReportR 14.12.0

## SciDataReportR 14.11.0

## SciDataReportR 14.10.0

## SciDataReportR 14.9.0

## SciDataReportR 14.8.0

## SciDataReportR 14.7.0

## SciDataReportR 14.6.0

## SciDataReportR 14.5.0

## SciDataReportR 14.4.0

## SciDataReportR 14.3.0

## SciDataReportR 14.2.0

## SciDataReportR 14.1.0

## SciDataReportR 14.0.0

## SciDataReportR 13.6.0

## SciDataReportR 13.5.0

## SciDataReportR 13.4.0

## SciDataReportR 13.3.0

## SciDataReportR 13.2.0

## SciDataReportR 13.1.0

## SciDataReportR 13.0.0

## SciDataReportR 12.22.0

## SciDataReportR 12.21.0

## SciDataReportR 12.20.0

## SciDataReportR 12.19.0

## SciDataReportR 12.18.0

## SciDataReportR 12.17.0

## SciDataReportR 12.16.0

## SciDataReportR 12.15.0

## SciDataReportR 12.14.0

## SciDataReportR 12.13.0

## SciDataReportR 12.12.0

## SciDataReportR 12.11.0

## SciDataReportR 12.10.0

## SciDataReportR 12.9.0

## SciDataReportR 12.8.0

## SciDataReportR 12.7.0

## SciDataReportR 12.6.0

## SciDataReportR 12.5.0

## SciDataReportR 12.4.0

## SciDataReportR 12.3.0

## SciDataReportR 12.2.0

## SciDataReportR 12.1.0

## SciDataReportR 12.0.0

## SciDataReportR 11.5.0

## SciDataReportR 11.4.0

## SciDataReportR 11.3.0

## SciDataReportR 11.2.0

## SciDataReportR 11.1.0

## SciDataReportR 11.0.0

## SciDataReportR 10.34.0

## SciDataReportR 10.33.0

## SciDataReportR 10.32.0

## SciDataReportR 10.31.0

## SciDataReportR 10.30.0

## SciDataReportR 10.29.0

## SciDataReportR 10.28.0

## SciDataReportR 10.27.0

## SciDataReportR 10.26.0

## SciDataReportR 10.25.0

## SciDataReportR 10.24.0

## SciDataReportR 10.23.0

## SciDataReportR 10.22.0

## SciDataReportR 10.21.0

## SciDataReportR 10.20.0

## SciDataReportR 10.19.0

## SciDataReportR 10.18.0

## SciDataReportR 10.17.0

## SciDataReportR 10.16.0

## SciDataReportR 10.15.0

## SciDataReportR 10.14.0

## SciDataReportR 10.13.0

## SciDataReportR 10.12.0

## SciDataReportR 10.11.0

## SciDataReportR 10.10.0

## SciDataReportR 10.9.0

## SciDataReportR 10.8.0

## SciDataReportR 10.7.0

## SciDataReportR 10.6.0

## SciDataReportR 10.5.0

## SciDataReportR 10.4.0

## SciDataReportR 10.3.0

## SciDataReportR 10.2.0

## SciDataReportR 10.1.0

## SciDataReportR 10.0.0

## SciDataReportR 9.17.0

## SciDataReportR 9.16.0

## SciDataReportR 9.15.0

## SciDataReportR 9.14.0

## SciDataReportR 9.13.0

## SciDataReportR 9.12.0

## SciDataReportR 9.11.0

## SciDataReportR 9.10.0

## SciDataReportR 9.9.0

## SciDataReportR 9.8.0

## SciDataReportR 9.7.0

## SciDataReportR 9.6.0

## SciDataReportR 9.5.0

## SciDataReportR 9.4.0

## SciDataReportR 9.3.0

## SciDataReportR 9.2.0

## SciDataReportR 9.1.0

## SciDataReportR 9.0.0

## SciDataReportR 8.23.0

## SciDataReportR 8.22.0

## SciDataReportR 8.21.0

## SciDataReportR 8.20.0

## SciDataReportR 8.19.0

## SciDataReportR 8.18.0

## SciDataReportR 8.17.0

## SciDataReportR 8.16.0

## SciDataReportR 8.15.0

## SciDataReportR 8.14.0

## SciDataReportR 8.13.0

## SciDataReportR 8.12.0

## SciDataReportR 8.11.0

## SciDataReportR 8.10.0

## SciDataReportR 8.9.0

## SciDataReportR 8.8.0

## SciDataReportR 8.7.0

## SciDataReportR 8.6.0

- Initial CRAN submission.
