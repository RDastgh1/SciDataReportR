# Onboarding Series — Slide Outlines

Five Quarto Beamer decks introducing the lab pipeline and `SciDataReportR` to a
new postdoc. Outlines only — slide titles, content stubs, and the code each
slide should demo. Prose to be written by RD.

Running data: `SampleData` + `SampleVariableTypes` for Decks 1–3 and most of 4.
`SimulatedPhenotypeData` + `SimulatedPhenotypeVariableTypes` for clustering.
External/simulated data pulled in where noted.

`# TODO:` marks something only RD can fill in (server paths, IRB numbers, study
names, real screenshots).

---

## Shared setup

Proposed file layout:

```
talks/onboarding-series/
  _brand.yml            # JHU heritage blue palette
  _common.qmd           # shared setup chunk, included by each deck
  preamble.tex          # metropolis tweaks, frame numbering, code font
  01-how-we-work.qmd
  02-codebook.qmd
  03-describe-and-screen.qmd
  04-model-and-project.qmd
  05-longitudinal-and-reporting.qmd
  figures/              # cached/static images
```

Add `^talks$` to `.Rbuildignore` so the decks don't ship with the package.

Shared YAML block (same in all five, only `title`/`subtitle` change):

```yaml
---
title: "Deck N — Title"
subtitle: "SciDataReportR onboarding series"
author: "Raha Dastgheyb"
institute: "Johns Hopkins Neurology"
date: today
format:
  beamer:
    aspectratio: 169
    theme: metropolis
    themeoptions:
      - progressbar=frametitle
      - numbering=fraction
    colortheme: default
    include-in-header: preamble.tex
    slide-level: 2
    fig-align: center
    fig-width: 7
    fig-height: 3.6
    code-line-numbers: false
    highlight-style: github
execute:
  echo: true
  warning: false
  message: false
  cache: true
  freeze: auto
---
```

`preamble.tex` sets metropolis accent to JHU heritage blue (`#002D72`),
secondary `#68ACE5`, and shrinks `\ttfamily` for code frames.

Every code-and-output slide follows the same two-column pattern so the
visual rhythm is predictable:

```
:::: {.columns}
::: {.column width="45%"}
<code chunk, echo: true, eval: false>
:::
::: {.column width="55%"}
<same call, echo: false, eval: true>
:::
::::
```

---

# Deck 1 — How We Work

*No package content. This is orientation: what the lab studies, where data
lives, and how a project is laid out. ~20 slides, ~25 min.*

## 1.1 Title

## 1.2 What this series covers
Five decks, one sentence each. Where to find them afterward.

## 1.3 What the lab studies
`# TODO:` current active studies, cohorts, and the questions behind them.
One slide, names + one-line aim each — enough that the methods later have
somewhere to land.

## 1.4 The shape of our data
Recurring structural features across studies: participants × visits,
labelled variables, mixed continuous/ordinal/categorical, wide feature
blocks (proteomics, imaging, cognition), missingness that means something.
This slide motivates every design decision in the package.

## 1.5 The problem this pipeline solves
Three recurring costs: re-writing EDA per study, metadata drifting away
from the data, and "can you re-run this on the new cohort?" being a
week of work instead of an afternoon.

## 1.6 Data governance — where things live
`# TODO:` raw data location, who can access it, how it gets there.
Rules: raw is read-only, never edited in place; nothing identifiable leaves
the secure location; analysis reads from a de-identified copy.

## 1.7 De-identification and what you may share
`# TODO:` de-ID standard used, date-shifting policy, what can go in a
manuscript figure vs. a supplement vs. a repo.
Explicit: never commit data to git, ever.

## 1.8 Project anatomy
```r
library(SciDataReportR)
CreateProjectFolders()
```
Show the resulting tree. Raw / Clean / Output separation and why the
boundary between them is one function call (`RevalueData()`).

## 1.9 `.Rproj` and `here()`
Never `setwd()`. Every path is `here::here("Data", "Clean", ...)`.
Why: the analysis has to run on her machine, RD's machine, and a cluster.

## 1.10 Version control: git *and* dated filenames
Both, on purpose. Git is the audit trail; `_2026_06_17.qmd` suffixes are
informal checkpoints when a shared analysis is in flight. Not redundant —
don't consolidate them.

## 1.11 What goes in the repo
Code, codebooks, `.qmd`, rendered HTML. Not data, not `.RData`, not
outputs with participant-level rows. Show a `.gitignore` starter.

## 1.12 Style rules that actually bite (1/2)
- `<-` for assignment, `=` only for arguments
- `%>%`, never `|>` — deliberate, don't modernize
- Every chunk gets a `#| label:`
- Modern chunk options (`#| echo: false`), not legacy inline

## 1.13 Style rules that actually bite (2/2)
- Naming: `df_` for data frames, `vars_` for name vectors, PascalCase for
  helper functions
- No `StepN:` anywhere — topical names
- No `# ──` banner dividers; headings carry structure
- Don't delete commented-out code
- `# TODO:` / `# QUESTION:` so notes are greppable

## 1.14 Figure captions go in headings
The one that surprises everyone: no `ggtitle()`, no `labs(title=)`, no
`fig-cap`. Context lives in the Markdown heading above the chunk.
Side-by-side: wrong vs. right.

## 1.15 The standard analysis skeleton
Ordered list, matching the house convention:
LoadPackages → LoadData → provenance callout → `RevalueData()` → EDA →
codebook display → codebook updates → dimension reduction → topic-named
analysis sections → `sessionInfo()` with `save.image()` commented out.

## 1.16 Two ways to start
```r
use_GetStartedScript()   # bare project bootstrap
use_EDATemplate()        # full EDA .qmd scaffold
```
Note which one to use when.

## 1.17 Reproducibility hygiene
`sessionInfo()` at the end of every file. `renv` — `# TODO:` confirm whether
the lab uses it. Seeds set explicitly in anything stochastic.

## 1.18 Where to look things up
pkgdown site, `vignette(package = "SciDataReportR")`, `?FunctionName`,
GitHub issues. Note that older function names (`CreatePCATable()`,
`Make_DataDictionary()`, `CreateZScorePlot()`) still work as aliases, so old
lab code will look different from the docs.

## 1.19 Day-one checklist
Install R + RStudio + Quarto → `devtools::install_github()` → request data
access `# TODO:` → clone the template repo → run `use_EDATemplate()` on
`SampleData` → render it.

## 1.20 Exercise
Set up a project folder, install the package, render the EDA template
against `SampleData` unmodified. Deliverable: the rendered HTML.

---

# Deck 2 — The Codebook Is the Analysis

*The central deck. If only one lands, make it this one. ~22 slides, ~30 min.*

## 2.1 Title

## 2.2 The claim
The variable-types file is not documentation. It is the analysis
specification, and every function downstream reads from it. Hand-coding
around it is the single most common way to fight the package.

## 2.3 What the codebook controls
One slide, four columns: labels (everything human-facing), types
(continuous / categorical / binary / ordinal), missing codes, and
inclusion flags (Table 1, Exclude, Category, Domain). Arrows from each to
the functions it governs.

## 2.4 Start from the data
```r
data("SampleData")
df_Raw <- SampleData
VariableTypes_Template <- CreateVariableTypesTemplate(df_Raw)
```
Show what the generated template looks like — mostly guesses that need
human review.

## 2.5 The manual edit
Screenshot the CSV open in Excel/RStudio. Walk column by column:
Variable, Label, Type, Category, Domain, Exclude, missing-code columns.
`# TODO:` real screenshot from a lab codebook (de-identified).

## 2.6 Writing good labels
Labels must stand alone without the variable name for context.
Bad: `"with them"`. Good: `"Felt that family was emotionally supportive"`.
Reason: the label is what appears on a poster axis, and no one will fix it
at submission time.

## 2.7 Types are decisions, not facts
Ordinal is the interesting case — a 0–4 symptom scale can be modeled either
way. `TreatOrdinalAs = "Continuous"` vs `"Categorical"` changes the answer.
Make the choice once, in the codebook, and record why.

## 2.8 Missing codes and the 999 trap
```r
df_Revalued <- RevalueData(df_Raw, VariableTypes)$RevaluedData
```
Show `SampleData$age` before and after: `999` is a sentinel, and a mean age
computed before `RevalueData()` is silently wrong. This is the slide that
scares people into using the codebook.

## 2.9 Check the warnings
`RevalueData()` returns more than the data frame. Show the full return
object and what to inspect — unmatched variables, unmapped levels.

## 2.10 Where `RevalueData()` sits
Restate the boundary: Raw → `RevalueData()` → Clean. Nothing upstream of
that call touches an analysis. Everything downstream assumes it ran.

## 2.11 No codebook? That's fine
Legitimate case: a one-off dataset with no codebook gets recoded manually.
That is not a shortcut or a violation — but say so in the file.

## 2.12 Seeing the codebook in a report
```r
df_DataDictionary <- MakeDataDictionary(df_Revalued)
FormattedDataDictionary(df_Revalued)
```
Read-only display for collaborators. Note `UpdateDataDictionary()`.

## 2.13 Merging: verify, never assume
Framing slide. Every join is a claim about record structure, and the claim
is usually slightly wrong.

## 2.14 Validating a merge
```r
df_Merged <- safe_merge(df_A, df_B, by = "ID")
ValidateMerge(df_A, df_B, df_Merged, keys = "ID",
              expected_relationship = "one-to-one")
```
Show `PlotMergeValidation()` output. Emphasize `expected_relationship` —
you state your expectation and the function checks it.
`# TODO:` build a small two-frame example with a deliberate duplicate key.

## 2.15 Exploring a bad merge
`ExploreMergeValidation()`, `merge_summary_table()`, `merge_detail()`.
What to document: unmatched IDs, duplicated keys, rows gained, columns
gained.

## 2.16 Longitudinal joins
```r
Merge_ByClosestTime(df_Clinical, df_Labs,
                    TimeVar1 = "VisitDate", TimeVar2 = "DrawDate",
                    keys = "ID", is_date = TRUE)
```
Also `MergeFragmentedRecords()` for split records.
Full longitudinal treatment is Deck 5.

## 2.17 Comparing two versions of a dataset
```r
CompareDatasets(df_v1, df_v2)
ExploreDatasetComparison(...)
PlotDatasetComparison(...)
```
Use case: the data manager sends a re-export and you need to know what moved.

## 2.18 Missingness
```r
PlotMissingData(df_Revalued, Relabel = TRUE)
```
Reading the plot: variable-level vs. participant-level, MCAR/MAR/MNAR in
one sentence each, and which variables are now off the table.

## 2.19 The codebook is a living document
As you derive variables — recodes, composites, cluster assignments — they
go into the codebook too:
```r
AddToCodebook(...); UpdateCodebook(...)
CombineCodebooks(...); MergeCodebooks(...)
```
Mention `CodebookMergeApp()` for the interactive case.

## 2.20 Other cleaning tools
`ReplaceMissingCode()`, `ReplaceMissingLabels()`, `ReValueFactors()`,
`ConvertOrdinalToNumeric()`, `reverseFactorLevels()`, `windsorize()`,
`IQROutliers()`. One line each — a menu, not a lesson.

## 2.21 Importing from elsewhere
`ReadSciData()`, `InspectFile()`, `PrepSPSS()`. Flag the SPSS label
truncation limits (variable and value labels get cut; `label_map` recovers
the full text).

## 2.22 Exercise
Take `SampleData`, deliberately break the codebook (mislabel a type, drop a
missing code), re-run `RevalueData()`, and describe what went wrong
downstream and how you'd have caught it.

---

# Deck 3 — Describe and Screen

*From clean data to a defensible shortlist. ~22 slides, ~30 min.*

## 3.1 Title

## 3.2 Framing: screen wide, model narrow
The two-stage logic of most of our papers. Screening output is a hypothesis
list, not a results section. State this early and repeat it.

## 3.3 Getting variable sets
```r
vars_Continuous  <- getNumVars(df_Revalued, Ordinal = FALSE)
vars_Categorical <- getCatVars(df_Revalued)
vars_Binary      <- getBinaryVars(df_Revalued)
```
And pulling sets from the codebook:
```r
vars_Table1 <- VariableTypes %>% filter(Table1 == 1) %>% pull(Variable)
```
This is where the codebook pays for itself.

## 3.4 Continuous distributions
```r
PlotContinuousDistributions(df_Revalued, vars_Continuous[1:12], ncol = 3)
```
What to look for: skew, floor/ceiling effects, impossible values that
survived cleaning, bimodality.

## 3.5 Categorical distributions
```r
PlotCategoricalDistributions(df_Revalued, vars_Categorical)
InspectCategoricalSummary(...)
```
Sparse cells and why they break tests later.

## 3.6 Summary table
```r
CreateSummaryTable(df_Revalued, vars_Continuous, Relabel = TRUE)
```
Note `Relabel = TRUE` — the codebook labels appear automatically.

## 3.7 Table 1
```r
MakeTable1(df_Revalued, TreatOrdinalAs = "Continuous")
```
Stratified version. Point out that the ordinal choice from Deck 2 shows up
here as different reported statistics.

## 3.8 Group comparison tables
```r
MakeComparisonTable(df_Revalued, ...)
MakeFacetCatComparisonTable(...)
CreateStatisticsTable(df_Revalued, TargetVar = "Diagnosis")
```
When to use which.

## 3.9 One relationship at a time
```r
PlotAssociations(df_Revalued, "age", "Adiponectin")
PlotAssociations(df_Revalued, "Diagnosis", "Ab_42")
```
Same function, different variable types, appropriate plot and test chosen
automatically. Show both outputs side by side.

## 3.10 Two-group focus plots
```r
Plot2GroupStats(...)
PlotSplitViolin(...)
```
Publication-ready single comparisons.

## 3.11 Scaling up: correlation heatmaps
```r
correlation_result <- PlotCorrelationsHeatmap(
  df_Revalued,
  xVars = vars_Continuous[1:5],
  yVars = vars_Continuous[20:40],
  method = "pearson", covars = NULL, Relabel = TRUE, Ordinal = FALSE
)
```
Emphasize: returns a *structured object*, not just a plot. Show `names()`
of the result.

## 3.12 Adjusting while screening
`covars =` in the heatmap; `PlotChiSqCovar()`, `PlotPartialRegressionScatter()`.
Partial correlations as a screening-stage adjustment.

## 3.13 Mixed-type screening
```r
PlotMiningMatrix(df_Revalued,
                 outcome_vars = vars_Outcomes,
                 predictor_vars = vars_Predictors,
                 covariates = NULL, TreatOrdinalAs = "Categorical")
```
Handles continuous, categorical, and ordinal predictors in one matrix.
This is the workhorse for a first pass on a new dataset.

## 3.14 Multiplicity
```r
ApplyFDRCorrection(pmat, fdr_scope = "per_outcome")
```
The three scopes — `"matrix"`, `"per_outcome"`, `"per_predictor"` — and how
to justify the choice. Point out that `PlotMiningMatrix()` and
`MultivariableRegressionTable()` take `fdr_scope` / `FDR` directly.
Slide should make the reviewer-facing argument, not just the syntax.

## 3.15 Direction, not just significance
```r
PlotDirectionalHeatmaps(...)
PlotVolcanoEffects(...)
```
Effect size and sign carry the biology; p-values alone don't.

## 3.16 Other screening views
`MakePairwiseHeatmap()`, `PlotPhiHeatmap()`, `PlotPointCorrelationsHeatmap()`,
`PlotAnovaRelationshipsMatrix()`, `PlotPValueComparisons()`.
Menu slide, one line each.

## 3.17 Interaction screening
```r
PlotInteractionEffectsMatrix(...)
PlotNumInteractionEffectsMatrix(...)
PlotCatInteractionEffectsMatrix(...)
PlotInteractionEffectsContinuous(...)
```
Caution slide: interactions are underpowered; screening them multiplies the
multiplicity problem.

## 3.18 Colors and house style
```r
SciDataPalette(); scale_color_SciData(); scale_fill_SciData()
scale_color_pvalue(); scale_fill_pvalue()
```
Consistent p-value color scales across figures so a reader learns the
mapping once.

## 3.19 Assembling figures
```r
AssemblePlots(Plots, ncol = 2, CollectLegend = TRUE,
              Theme = ggplot2::theme_minimal())
```
Multi-panel figures with a shared legend. Reminder: panel context goes in
the surrounding heading, not in the plot.

## 3.20 What a screening result is worth
Restate the framing. A screening hit is an entry on a list to model
properly in Deck 4. Show an example of what *not* to write in an abstract.

## 3.21 Exercise
Screen `SampleData` biomarkers against `Diagnosis` with `PlotMiningMatrix()`,
apply FDR per outcome, produce a shortlist of five, and write one paragraph
justifying the FDR scope you chose.

---

# Deck 4 — Model and Project

*The payoff deck. Formal models, then the reusable-transformation pattern.
~26 slides, ~35 min.*

## 4.1 Title

## 4.2 Where we are
Pipeline strip: clean → described → screened → **modeled**. Shortlist from
Deck 3 is the input.

## 4.3 Univariate regression tables
```r
MakeUnivariateRegressionTable(df_Revalued,
                              outcome_vars = vars_Outcomes,
                              predictor_vars = vars_Shortlist,
                              covariates = c("age", "sex"),
                              Standardize = TRUE, FDR = TRUE)
```
One model per outcome–predictor pair, tabulated. Note the
`UnivariateRegressionTable()` alias in older lab code.

## 4.4 Reading the table
Standardized vs. raw coefficients, CIs, adjusted p-values, N per model
(which varies with missingness — call this out).

## 4.5 Multivariable models
```r
MultivariableRegressionTable(df_Revalued,
                             outcome_vars = "Outcome",
                             predictor_vars = vars_Shortlist,
                             covariates = c("age", "sex"),
                             Standardize = TRUE, Method = "lm", FDR = TRUE)
```

## 4.6 When predictors outnumber sensible degrees of freedom
```r
MultivariableRegressionTable(..., Method = "lasso",
                             CVFolds = 10, Lambda = "lambda.min", Seed = 123)
```
Ridge / lasso / elastic net live in the same function. Set `Seed` or CV
results move between runs.

## 4.7 Forcing covariates into a penalized model
The key detail: covariates passed via `covariates =` are forced in with
`penalty.factor = 0`, so age and sex are never shrunk out. Explain why
that matters for interpretation.

## 4.8 Forest plots
```r
PlotForestFromTable(regression_table, ...)
```
Takes the table object from either regression function. Show the
table → figure handoff explicitly.

## 4.9 Biomarker performance
```r
EvaluateBiomarkerPerformance(df_Revalued,
                             outcome_var = "Diagnosis",
                             biomarker_var = "Ab_42",
                             covariates = c("age", "sex"),
                             ThresholdMethod = "youden",
                             Validation = "bootstrap", Seed = 123)
```
Walk the returned object: ROC, AUC with CI, threshold, calibration.

## 4.10 Thresholds are a clinical choice
`"youden"` vs. `"sensitivity"` vs. `"specificity"` vs. `"custom"` — the
right one depends on the cost of a miss. One slide, no code.

## 4.11 Validation
`Validation = "bootstrap"` / `"cross_validation"`. Optimism and why an
apparent AUC from the training data is not a result.

## 4.12 Screening many biomarkers
```r
ScreenBiomarkerPerformance(...)
```
Plus `PlotBlandAltman()` for method-comparison / agreement questions.

## 4.13 Transition: from testing to structure
Framing slide. Everything so far asked "is X related to Y?" The rest of the
deck asks "what structure is in this feature block?"

## 4.14 PCA
```r
pca_obj <- CreatePCAObject(df_Revalued, variables = vars_Features)
plotPCA(pca_obj)
ExtractPCAComponentSummary(pca_obj)
```
Scree plot, loadings lollipop, how many components to keep.
(`CreatePCATable()` is the old name.)

## 4.15 Interpreting components
Loadings → a scientific name for the component. Warn against
over-narrating PC3.

## 4.16 Categorical structure
```r
CreateMCAObject(...)
```
MCA as the categorical analogue. When it's the right tool.

## 4.17 Standardizing within a cohort
```r
z_obj <- CreateZScoreObject(df_Revalued, variables = vars_Features)
PlotZScore(...)
CalcZScore(...); CalcMScore(...); CreateMScoreObject(...)
```

## 4.18 Normative modeling
```r
model_T <- CreateNormativeTScoreModel(df_Revalued,
                                      test_var = "...", count_var = "...",
                                      covariates = c("age", "sex", "education"),
                                      reference_var = "Group",
                                      reference_value = "Control",
                                      include_practice_effect = FALSE)
ApplyNormativeTScores(df_New, model_T)
```
Demographically-adjusted T-scores against a reference group. Note the
practice-effect option for repeated cognitive testing.

## 4.19 Reliable change
```r
rci <- CreateRCIObject(df_Revalued, variables = vars_Cognition,
                       DataFormat = "wide", id_var = "ID",
                       Method = "regression",
                       BaselineSpecifier = "_v1", FollowupSpecifier = "_v2")
ProjectRCI(rci, df_New)
```
Has this person actually changed, or is it measurement noise?

## 4.20 **The Create / Project pattern**
The idea the whole series builds toward. Every `Create*Object()` returns a
fitted transformation; every `Project*()` applies it, unchanged, to new
data. `CreatePCAObject()` → `ProjectPCA()`. `CreateZScoreObject()` →
`ProjectZScore()`. `CreateRCIObject()` → `ProjectRCI()`.
Diagram: cohort A trains, cohort B is scored. No refitting, no leakage.

## 4.21 Why it matters
Three consequences: no data leakage into a validation cohort, the
transformation is a saveable object you can archive with the paper, and
"re-run this on the new cohort" becomes one function call.

## 4.22 Clustering — switch to `SimulatedPhenotypeData`
```r
data("SimulatedPhenotypeData")
data("SimulatedPhenotypeVariableTypes")
```
Why a different dataset: phenotype structure is the point, and it's
simulated so the ground truth is known.

## 4.23 The pipeline menu
`CreateClusterModel_KMeans()`, `_MClust()`, `_PCA_KMeans()`,
`_PCA_MClust()`, `_HDBSCAN()`, `_SOM_MClust()`, `_SOM_HDBSCAN()`,
`_Gower_PAM()` (mixed types), `_LatentClass()`, `_MCA_MClust()`
(categorical). Table: input type → reduction → cluster method.

## 4.24 Fitting one
```r
cluster_model <- CreateClusterModel_SOM_MClust(df_Phenotype,
                                               variables = vars_Features)
```
Show the SOM → MClust logic. `Pipeline_SOMClust()` for the wrapper.

## 4.25 Is the clustering real?
```r
PlotClusterDiagnostic(); PlotClusterSilhouette(); PlotClusterFitReview()
```
The honest slide: k-means always returns k clusters. Stability and
diagnostics before interpretation.

## 4.26 Describing clusters
```r
PlotClusterProfiles(); PlotClusterCentreHeatmap(); PlotClusterCentreProfile()
PlotClusterComposition(); PlotClusterBoxplot(); PlotClusterMap()
PlotSpiderChart()
```
Turning cluster indices into phenotypes with names.

## 4.27 Projecting clusters onto a new cohort
```r
ProjectCluster(cluster_model, df_NewCohort)
```
The pattern again, at the level of a whole pipeline. S3 dispatch means the
call is identical regardless of which pipeline was fitted.

## 4.28 Feeding clusters back
New cluster assignment → add to the data → **add to the codebook**
(`AddToCodebook()`) → it's now a variable like any other and can go into
Deck 3's screening tools. Close the loop back to Deck 2.

## 4.29 Exercise
Fit two different cluster pipelines to `SimulatedPhenotypeData`, compare
diagnostics, pick one, project it onto a held-out split, and report whether
the phenotypes replicate.

---

# Deck 5 — Longitudinal Data and Reporting

*Visit structure, change over time, and getting work out the door.
~18 slides, ~25 min.*

## 5.1 Title

## 5.2 Why longitudinal gets its own deck
Repeated visits break the one-row-per-participant assumption behind most of
Decks 3 and 4. Wide vs. long, and which functions expect which.

## 5.3 Our visit structure
`# TODO:` actual visit schedule for the lab's cohorts — planned windows,
what "baseline" means, how much drift is tolerated.

## 5.4 Seeing the time structure
```r
PlotTimeDistribution(df_Long, ...)
```
Actual visit timing vs. protocol. Nearly always messier than the protocol
implies — this slide sets expectations honestly.

## 5.5 Aligning measurements in time
```r
Merge_ByClosestTime(df_Clinical, df_Labs,
                    TimeVar1 = "VisitDate", TimeVar2 = "DrawDate",
                    keys = "ID", is_date = TRUE)
```
The unspoken decision: how far apart is too far? Document the window.

## 5.6 Fragmented records
```r
MergeFragmentedRecords(...)
```
One participant, several partial rows. Common with REDCap exports.

## 5.7 Status over time
```r
SummarizeTransitions(df_Long, id_var = "ID", time_var = "Visit",
                     status_var = "Symptomatic",
                     order_participants_by = "first_positive")
```
Counts of transitions, time to first event, persistence.

## 5.8 Swimmer plots
```r
PlotSwimmerTransitions(df_Long, id_var = "ID", time_var = "Visit",
                       status_var = "Symptomatic",
                       x_axis_type = "time_from_baseline",
                       time_from_baseline_unit = "months")
```
Reading it: one row per participant, color = state. Ordering choices change
the story the figure tells — pick one and say why.

## 5.9 Events and states
```r
PlotTimeSwimmer(df_Long, id_var = "ID", Time = "Date",
                State = "Status", Event = "Hospitalization",
                TimeScale = "from_event", Format = "state_path")
```
`TimeScale = "from_event"` for anchoring on infection/treatment date.

## 5.10 Change scores done properly
Callback to `CreateRCIObject()` / `ProjectRCI()` from Deck 4 — the right way
to ask whether an individual changed. Contrast with naive follow-up-minus-
baseline.

## 5.11 What we are *not* covering
Mixed models / GEE are outside the package. `# TODO:` name the lab's
current approach (`lme4`? `nlme`?) and where the template code lives.
Being explicit about the boundary prevents her looking for a function that
isn't there.

## 5.12 Domain-specific derivations
```r
DeriveFreesurferVolumes(...)
calculate_pathway_results(); PlotPathway_KT(); CreatePathwayPlot_KT()
```
`# TODO:` confirm which of these she'll actually need.

## 5.13 Report-ready tables
```r
CreateSummaryTable(); CreateStatisticsTable(); FreezeTableHeader()
InsertValues()
```
`InsertValues()` for inline numbers in prose — the way to never have a stale
number in a manuscript.

## 5.14 Rendering and sharing
Quarto output targets: HTML for collaborators, `docx` when a clinician needs
to track-change it, PDF for supplements. `#| column: screen` for wide
matrices. `freeze` so a re-render doesn't recompute everything.

## 5.15 What ships with an analysis
The rendered report, the code, the codebook version used, `sessionInfo()`,
and the provenance callout naming the data file. Anything less and the
result isn't reproducible in six months.

## 5.16 Contributing back
When a helper gets reused twice, it becomes a package function. Point at
`CONTRIBUTING.md`, roxygen docs, `devtools::check()`, the vignette
convention. Frame it as expected, not optional.

## 5.17 The whole pipeline on one slide
Single diagram tying all five decks together. This is the slide to print and
tape above her monitor.

## 5.18 Capstone exercise
`# TODO:` pick a small real (de-identified) dataset. Full pass: codebook →
clean → validate merges → screen → model → report. Deliverable is a
rendered `.qmd` that RD can run unmodified.

---

## Open questions for RD

1. Deck 1 slides 1.3, 1.6, 1.7 and Deck 5 slide 5.3 need lab-specific
   content that isn't in the repo.
2. Does the lab use `renv`? Changes slide 1.17 and the day-one checklist.
3. Is `InspectMergeKeys()` / `HarmonizeMergeKeys()` still exported? It has an
   `.Rd` and an `.R` file but doesn't appear in `NAMESPACE` — slide 2.14
   currently avoids it.
4. Deck 5 slide 5.11 — confirm the lab's mixed-model approach so the
   boundary is stated accurately.
5. Deck 4 is the longest at ~26 slides. Splitting regression/biomarkers from
   dimension-reduction/clustering would make six decks; happy either way.
