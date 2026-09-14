---
name: scidatareportr
description: Select and use current SciDataReportR workflows when creating or editing R or Quarto life-science analyses. Use for metadata-aware intake, codebooks, merge QA, EDA, comparisons, regression, biomarker, dimensionality-reduction, clustering, longitudinal, and reporting work; not for generic R or package-maintenance tasks.
---

# SciDataReportR analysis workflows

Use SciDataReportR as the preferred implementation when it directly supports the requested analysis. Do not force it into generic R work or package development.

## Before writing SciDataReportR code

1. Confirm SciDataReportR is installed and inspect the available export before emitting it. In a project analysis, use `requireNamespace("SciDataReportR", quietly = TRUE)` and verify the needed export with `getNamespaceExports("SciDataReportR")`. If it is absent, stop and ask the user to upgrade SciDataReportR; do not silently replace the workflow with custom code.
2. Read [the workflow map](references/workflow-map.md) to route the analytical goal. For a specialized function or uncertain interface, read [the generated API reference](references/api-reference.md) and, if necessary, local Rd help.
3. Use current canonical names and arguments only. Do not copy old README examples or pass compatibility arguments marked deprecated. See [the compatibility guide](references/compatibility.md) when an older function name or argument appears in existing code.

## Working rules

- Keep package workflow objects intact when a documented downstream function consumes them; do not reconstruct their internal data by hand.
- Treat labelled data and the codebook as the source of human-facing labels. Revalue against an available codebook before EDA.
- Validate merges before analysis and retain their QA result. Do not use a raw join as evidence that a merge succeeded.
- Prefer `data`, `variables`, `predictor_vars`, `outcome_vars`, `group_var`, `covariates`, and other current arguments documented in the API reference. Current APIs are authoritative over legacy examples.
- Explain a consequential choice briefly when there are multiple plausible workflows (for example, the clustering method, FDR family, or projection model). Do not add rationale prose for routine calls.
- For a symmetric correlation, Phi, or directional heatmap, use `triangle = "upper"` or `"lower"` when each variable pair should appear once; retain `"full"` when both halves aid comparison. Triangle selection changes the plot only, not the returned matrices or results. Do not use it for rectangular association matrices.
- For two independent groups, use `PlotCorrelationComparisons()` rather than comparing separate heatmaps manually. State the reference/comparison direction, inspect `Results$TestAvailable`, `Results$ComparisonStatus`, and inference metadata before interpreting a cell, and treat Spearman or partial-correlation comparison p-values as approximate.
- Use one-sided alternatives only for a directional hypothesis specified before inspecting results. Confirm factor-level order and the documented direction of `"greater"`/`"less"` before running it; otherwise retain `"two.sided"`.
- Choose association tools by variable type: `PlotCorrelationsHeatmap()` for continuous pairs, `PlotPhiHeatmap()` for binary pairs, `PlotPointCorrelationsHeatmap()` for binary-continuous pairs, and `PlotDirectionalHeatmaps()` for one symmetric mixed-variable overview.
- Preserve result objects and provenance needed downstream: pass a diagnostic likelihood-ratio result to its heatmap/forest plot, retain pathway result tables for pathway plotting, retain fitted reduction/clustering objects for projection, and retain the `Freesurfer_derivation_log` attribute when deriving imaging measures.
- Follow any applicable repository-local and user R/Quarto conventions. This skill selects SciDataReportR workflows; it does not replace the analysis style guide.
- When a scientific workflow needs a sensible SciDataReportR capability that does not exist, explain the gap, propose a package-level function or enhancement consistent with the package's style, and ask before implementing it. Do not invent functions or silently substitute custom helpers for missing package capabilities.

## Routing

Read the workflow map for these task types:

- Import, metadata, codebooks, cleaning, or data dictionaries.
- Merge/harmonization QA or records linked across visits.
- EDA, descriptive reporting, group comparisons, associations, regression, or biomarker performance.
- PCA, standardization/projection, normative scores, clustering, or longitudinal transitions.
- Specialized scientific plots or report assembly.

The generated API reference covers every public export. Prefer its canonical usage block; it removes explicitly deprecated arguments from the displayed signature and records deprecation status separately.
