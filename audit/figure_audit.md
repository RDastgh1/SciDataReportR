# Reference-page figure audit

Scope: every example figure rendered into `docs/reference/` (177 referenced PNGs
across 76 pages), reviewed visually. Grouped by severity. Nothing has been
changed — this is a list of suggested changes only.

---

## 1. Rendering bugs — the figure is wrong or unreadable as drawn

| # | Where | Problem | Suggested fix |
|---|---|---|---|
| 1.1 | `PlotClusterDiagnostic-1`, `PlotProjectionDiagnostics-1` | Every y-axis tick label reads `1`. Posterior-max values are ~0.9999x and get rounded to a single digit, so three different gridlines all say "1". | Format the axis with enough significant digits (`scales::label_number(accuracy = 1e-4)`), or plot `1 - posterior_max` on a log scale. |
| 1.2 | `AssemblePlots-1` | Panel A collapses: p-value annotations overprint each other (`p = 0.142.12e-09…`), x tick labels merge (`E2E2E2E3E2E4…`), n-labels overlap. Panel B's marginal axis reads `0102030`. | Render the example at a larger `fig.width`/`fig.height`, or assemble two plots that survive half-width (no ggstatsplot annotations). |
| 1.3 | `PlotPathway_KT-1`, `-2` | Node labels overflow their boxes and collide: "Kynurenic Acid"/"3-Hydroxykynurenine"/"Anthranilic Acid" overlap on one row; "N-Formylkynurenine" spills outside its box. | Wrap labels, reduce `cex`, size boxes to `strwidth()`, and widen horizontal spacing between sibling nodes. |
| 1.4 | `PlotDirectionalHeatmaps-1`, `-2` | Three stacked legends exceed the panel height — the top legend is clipped, showing only `-0.5 / -1.0` with no title. | Collect the guides horizontally, or drop to a single scale and facet by statistic type. |
| 1.5 | `windsorize-1` | Facet strips print the label twice: "raw / raw", "winsorized / winsorized". | The labeller emits name + label even when they are identical; suppress the second line when `name == label`. |
| 1.6 | `PlotPartialRegressionScatter-1` | Subtitle reads `Residual AXL = 0 + 0 * age_resid` — both coefficients rounded to zero. Also the x-axis is labelled "Age" but the values are residuals (−40…40). | Use `signif(, 3)` for the fitted equation; label the axis "Age \| Adiponectin (residual)". |
| 1.7 | `CreatePCAObject-1`, `CreateMCAObject-1` | Y-axis says "Cumulative Var" / "cumulative percentage of variance", but the bars are the *per-component* proportion. Only the overlaid line is cumulative. | Retitle the axis ("Proportion of variance"), or add a secondary axis / annotate the line as the cumulative series. Label the dashed threshold line. |
| 1.8 | `PlotAssociations-1` | The whole panel background is painted with the `NA` fill colour, so the figure reads as a red plot with blue bars. `NA` also appears as a legend key. | Drop `NA` from the fill scale (`na.translate = FALSE`) or set `na.value = "grey90"`; the panel background should never inherit a data colour. |

## 2. Legends with keys that draw nothing

A recurring pattern: the legend lists every factor level, but only the levels
present in the data get a glyph, leaving blank rows the reader can't interpret.

- `PlotTimeSwimmer-1` — `EventType` legend has a blank `NA` row.
- `PlotSwimmerTransitions-1` — `Status` legend has a blank `Missing` row.
- `PlotMergeValidation-1` — `Status` legend shows `WARNING` and `FAIL` with no swatch.
- `PlotPValueComparisons-1` — `Sig` legend lists `ns / * / ** / *** / ****`; only two have glyphs.
- `PlotAnovaRelationshipsMatrix-1`, `PlotChiSqCovar-1` — same, plus the legend is titled `p<.05` while listing five significance bands.

**Suggested fix:** drop unused levels from the guide, or keep them and draw a
greyed-out glyph so the key is legible. Retitle the significance guides
"Significance" rather than `p<.05`.

## 3. Figures that mislead

| # | Where | Problem | Suggested fix |
|---|---|---|---|
| 3.1 | `PlotClusterMap-1`, `-2` | The example maps clusters onto two raw variables, so clusters 1/4 and 2/3 sit completely on top of each other — the clustering looks like it failed. The same model plotted in the frozen review space (`CreateKMeansClusterModel-6`) separates cleanly. | Default the example (and ideally the function) to the model's `ReviewSpace`, or state in the caption that raw-variable axes are not the review space. |
| 3.2 | `PlotPhiHeatmap-1`, `PlotChiSqCovar-1` | Self-associations (Sex×Sex, Diagnosis×Diagnosis) are included at φ = 1 / V ≈ 1 and dominate the colour scale, flattening every real off-diagonal value to near-white. | Mask the diagonal (`NA`) as `PlotDirectionalHeatmaps` already does with grey tiles. |
| 3.3 | `PlotVolcanoEffects-1`, `-3` | Titled "Volcano plot using raw p-values" but the colour legend still reads "FDR significant". | Make the legend text follow the p-source: "p < .05" on the raw panel, "FDR significant" on the adjusted panel. |
| 3.4 | `CreateLatentClassClusterModel-1`, `PlotClusterFitReview-1` | Free y-scales on near-constant metrics: Entropy is drawn across 0.9986–0.9994 and "Min Cluster N" across 79.95–80.05, implying large differences where the metric is effectively flat. | Drop panels whose range is below a tolerance, or annotate them ("constant at 80"). |
| 3.5 | `PlotClusterFitReview-1`, `CreateLatentClassClusterModel-1` | X-axis shows 2.0, 2.5, 3.0, 3.5, 4.0 for an integer class count. | `scale_x_continuous(breaks = unique(fit_table[[x]]))` in `PlotClusterFitReview` (`R/ClusterPlots.R:617`). |
| 3.6 | `PlotTimeSwimmer-3` | Two legends both contain an entry called "Flare" with different colours (black point = `EventType`, salmon = `State`). The cluster palette also differs from `PlotTimeSwimmer-1`. | Rename the guides ("Event" vs "State at visit") and fix one palette across the function's panels. |

## 4. Examples whose data defeats the figure

| # | Where | Problem | Suggested fix |
|---|---|---|---|
| 4.1 | `CreateNormativeTScoreModel-3` (and siblings) | The example has n = 10 (8 reference, 2 clinical). T-scores should sit near 50 ± 10; the axis runs to −250 and the histogram has 4 bars of height 1 and 3. | Simulate a reference sample of a few hundred so the normative model produces a plausible distribution. |
| 4.2 | `plotSigCorrelations-1` | The example correlates `mpg` with `mpg`: r = 1, p = 0, a perfect diagonal line. The "Adjusted for" caption is also truncated/empty. | Correlate two different variables. |
| 4.3 | `PlotDatasetComparison-3` | Fully empty panel: four zeros with a y-axis from −0.05 to 0.05. | Use comparison data that actually adds/removes records, or emit the "no changes" placeholder the merge-validation plots use. |
| 4.4 | `PlotMergeValidation-4`, `-5` | Blank "No unresolved duplicate variables detected" placeholders. Correct behaviour, but as documentation they show nothing. | Use an example with an unresolved duplicate so the real figure is visible; keep one placeholder to document the empty state. |
| 4.5 | `PlotPointCorrelationsHeatmap-1` | All values are ≈ 0, so the heatmap is effectively blank white. | Pick variables with real point-biserial signal. |
| 4.6 | `PlotInteractionEffectsMatrix-1` | 2 non-NA cells in a 9 × 10 grid; the rest renders as blank white with no tile outline, so it doesn't read as a matrix. | Shade `ns` cells light grey so the grid is visible, and subset the example to variables with interactions. |
| 4.7 | `PlotMiningMatrix-1` | Effect-size scale fixed 0–1 while every point sits below 0.25, so all markers are the same yellow. Marker size is also inverse to significance (`***` renders smaller than `*`). | Let the effect-size scale follow the data range; make the size/shape ramp monotone in significance. |
| 4.8 | `PlotClusterOccupancy-1`, `CreateSOMClusterModel-3` | All four bars are identical (80, 25% each) because the simulated data is perfectly balanced. | Fine as a sanity check, but consider an unbalanced example so the plot demonstrates something. |

## 5. Redundant / near-duplicate figures

- `PlotCorrelationsHeatmap-1` vs `-2` — visually indistinguishable (raw vs FDR); only a handful of stars differ. Add a subtitle naming the correction, or show a single panel.
- `PlotDirectionalHeatmaps-1` vs `-2`, `PlotVolcanoEffects-1` vs `-2` — same issue.
- `PlotSpiderChart-1`, `-2`, `-4` — three near-identical radars differing only in axis ordering. Reduce to two and title them with the option being demonstrated.
- `CreateKMeansClusterModel-2`, `-3` (elbow, silhouette-by-k) reproduce panels already shown in `-1` (`fit_plot`). Consider showing the facetted review plus one detail curve, not both.

## 6. Label / styling consistency

- **Codebook labels are applied inconsistently.** Heatmaps and distribution plots use pretty labels ("AXL receptor tyrosine kinase"), while `PlotSplitViolin` (`Ab_42`), `PlotAssociations-2` (`age`), `Plot2GroupStats`, `PlotPValueComparisons` (`age`, `AXL`), and `PlotClusterComposition` (`CatVar1`) use raw column names. `Plot2GroupStats-1` even truncates one (`ACE_CD143_Angiotensin_Converti`). Route every axis/strip label through the same codebook lookup.
- **`ahp index`** is lower-cased in every fit-review facet while its neighbours are Title Case. `.ClusterMetricLabel` (`R/ClusterPlots.R`) prettifies generically; add an explicit map for `ahp_index` → "AHP index", `WSS`, `MinClusterN` → "Min cluster N".
- **Theme is not uniform.** Most plots use `theme_bw()`; `PlotPValueComparisons`, `PlotAnovaRelationshipsMatrix`, `PlotChiSqCovar`, `PlotPhiHeatmap`, `scale_color_pvalue` still render on the default grey theme.
- **Cluster palettes differ across the family.** Silhouette/map/composition use the ggplot hue palette, occupancy and persistence use a single blue, centre heatmap uses red/blue, projection-fit uses green/orange. Pick one cluster colour scale and reuse it.
- **Meaningless y-axes are labelled.** The raincloud plots (`PlotContinuousDistributions`, `windsorize`) show 0.8–1.4 tick labels that are layout offsets, not data. Blank them (`theme(axis.text.y = element_blank())`).
- **Facet strips truncate**: `Gamma_Interferon_induced_Monokin`, `ated apoptosis-inducing ligand re` in `PlotContinuousDistributions`. Wrap strip text.
- **Missing axis labels**: `PlotZScore-1` has no x tick labels at all (you cannot tell which variable is which) and two unexplained blue/green point series; `PlotTimeDistribution-1` has no x axis title or labels and shows time as `2024.25`; `PlotSplitViolin-1` has no x axis and an unexplained `***`.
- `PlotCategoricalDistributions-1`: small slices (E2E2, E2E4, E4E4) are unlabelled and there is no legend, so they are unidentifiable.
- `PlotDatasetComparison-1` / `PlotMergeValidation-1`: bar length encodes nothing (every bar is full width). A status tile grid with the status text inside would carry the same information more honestly.
- `MakePairwiseHeatmap-1`: the heavy black border on the Cortisol row is unexplained — add a caption for what it marks.
- `PlotForestFromTable-1`: everything is grey with very low contrast and only three widely-spaced rows. Colour by significance and tighten row spacing.
- `IQROutliers-1`: legend reads `FALSE`/`TRUE`; use `No`/`Yes`. Top outlier sits above the last y break.
- `geom_starcaption-1`: 90°-rotated column labels push the plot into the top third; angle at 45° or wrap.
- `CodebookMergeApp-harmonization.png`: the screenshot is clipped mid-column (`Resolutio…`, `Categoric…`). Re-capture at a wider viewport.

## 7. Housekeeping

- **109 orphaned PNGs** in `docs/reference/` are no longer referenced by any HTML page — leftovers from earlier builds (e.g. `CreateKMeansClusterModel-9…14`, all of `ProjectKMeansCluster-3…7`, `plotForestFromTable-1.png` from the old lower-case alias). They inflate the repo and the site bundle. Delete `docs/reference/` and rebuild, or prune the unreferenced files.
- No broken image references were found — every `<img src>` in `docs/` resolves.
