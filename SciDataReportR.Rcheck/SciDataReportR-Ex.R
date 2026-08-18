pkgname <- "SciDataReportR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SciDataReportR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AddToCodebook")
### * AddToCodebook

flush(stderr()); flush(stdout())

### Name: AddToCodebook
### Title: Add a new variable to a codebook
### Aliases: AddToCodebook

### ** Examples





cleanEx()
nameEx("ApplyFDRCorrection")
### * ApplyFDRCorrection

flush(stderr()); flush(stdout())

### Name: ApplyFDRCorrection
### Title: Apply multiple-comparison correction across a p-value matrix
### Aliases: ApplyFDRCorrection

### ** Examples

pm <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
             nrow = 2,
             dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3")))

# One family across the whole matrix (classic behavior)
ApplyFDRCorrection(pm)

# Correct each outcome (column) separately
ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 2)

# Vector input with explicit outcome grouping
ApplyFDRCorrection(c(0.01, 0.04, 0.02, 0.03),
                   fdr_scope = "per_outcome",
                   outcome_ids = c("y1", "y1", "y2", "y2"))

# A symmetric matrix: six filled cells, but only three pairs tested
vars <- c("mmse", "trails", "digit_span")
pm_sym <- matrix(NA_real_, 3, 3, dimnames = list(vars, vars))
pm_sym[lower.tri(pm_sym)] <- c(0.001, 0.020, 0.040)
pm_sym[upper.tri(pm_sym)] <- t(pm_sym)[upper.tri(pm_sym)]
pm_sym

# Detected automatically: three tests in the family, diagonal excluded
ApplyFDRCorrection(pm_sym)

# Bonferroni, corrected once per pair and then against all six cells
ApplyFDRCorrection(pm_sym, method = "bonferroni")
ApplyFDRCorrection(pm_sym, method = "bonferroni", symmetric = FALSE)




cleanEx()
nameEx("ApplyNormativeTScores")
### * ApplyNormativeTScores

flush(stderr()); flush(stdout())

### Name: ApplyNormativeTScores
### Title: Apply a normative T-score model to new data
### Aliases: ApplyNormativeTScores

### ** Examples

# A reference group and a clinical group tested on Trail Making A
set.seed(206)
n_reference <- 220
n_clinical <- 80

df <- tibble::tibble(
  Group = c(rep("Reference", n_reference), rep("Clinical", n_clinical)),
  Age = round(c(
    stats::rnorm(n_reference, 52, 12),
    stats::rnorm(n_clinical, 58, 12)
  )),
  Education = factor(sample(
    c("HighSchool", "College", "Graduate"), n_reference + n_clinical,
    replace = TRUE
  )),
  Sex = factor(sample(c("F", "M"), n_reference + n_clinical, replace = TRUE)),
  Visit = sample(1:3, n_reference + n_clinical, replace = TRUE)
)
df$TrailsA <- round(1000 * exp(stats::rnorm(
  nrow(df),
  mean = log(28) + 0.011 * (df$Age - 52) - 0.05 * (df$Visit - 1) +
    ifelse(df$Group == "Clinical", 0.35, 0),
  sd = 0.22
)))

# Fit the norms on the reference group only
norm_obj <- CreateNormativeTScoreModel(
  data = df,
  test_var = "TrailsA",
  count_var = "Visit",
  covariates = c("Age", "Education", "Sex"),
  reference_var = "Group",
  reference_value = "Reference",
  include_practice_effect = TRUE,
  reverse_score = TRUE,
  convert_seconds = TRUE,
  log_transform = TRUE,
  return_plots = FALSE
)

# Score everyone through the same model
scored_df <- ApplyNormativeTScores(
  data = df,
  normative_obj = norm_obj
)

# Before: raw completion times, grouped by clinical status
attr(scored_df$NormRaw, "label") <- "Trail Making A completion time (ms)"
PlotContinuousDistributions(
  scored_df, variables = "NormRaw", Fill = "Group", ncol = 1
)

# After: demographically adjusted T-scores. The reference group is centered
# near 50, and the clinical shift is now directly interpretable in SD units.
attr(scored_df$NormT, "label") <- "Demographically adjusted Trail Making A T-score"
PlotContinuousDistributions(
  scored_df, variables = "NormT", Fill = "Group", ncol = 1
)




cleanEx()
nameEx("AssemblePlots")
### * AssemblePlots

flush(stderr()); flush(stdout())

### Name: AssemblePlots
### Title: Assemble ggplot objects into a unified multi-panel figure
### Aliases: AssemblePlots

### ** Examples

data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

p_Categorical <- PlotAssociations(df_Labelled, "Diagnosis", "Genotype")
p_Continuous <- PlotAssociations(df_Labelled, "age", "AXL")

# Keep each association plot's own statistical annotation and legend.
AssemblePlots(
  list(
    "Diagnosis and genotype" = p_Categorical,
    "Age and AXL" = p_Continuous
  ),
  ncol = 2,
  CollectLegend = FALSE,
  Labels = c("A", "B")
)




cleanEx()
nameEx("CodebookMergeApp")
### * CodebookMergeApp

flush(stderr()); flush(stdout())

### Name: CodebookMergeApp
### Title: Interactive codebook harmonization dashboard
### Aliases: CodebookMergeApp

### ** Examples

## Not run: 
##D data(SampleVariableTypes)
##D 
##D # Two overlapping codebooks to harmonize before a deterministic merge
##D cb_a <- SampleVariableTypes[1:12, c("Variable", "Label", "Type")]
##D cb_b <- cb_a[-(1:2), ]
##D cb_b$Type[cb_b$Type == "Double"] <- "numeric"
##D 
##D # Launch the interactive harmonization dashboard
##D CodebookMergeApp(
##D   codebooks = list(CohortA = cb_a, CohortB = cb_b),
##D   VariableCol = "Variable"
##D )
## End(Not run)



cleanEx()
nameEx("CompareDatasets")
### * CompareDatasets

flush(stderr()); flush(stdout())

### Name: CompareDatasets
### Title: Compare two versions of a dataset
### Aliases: CompareDatasets

### ** Examples





cleanEx()
nameEx("ConvertOrdinalToNumeric")
### * ConvertOrdinalToNumeric

flush(stderr()); flush(stdout())

### Name: ConvertOrdinalToNumeric
### Title: Prepare ordinal variables for analysis
### Aliases: ConvertOrdinalToNumeric

### ** Examples





cleanEx()
nameEx("CreateClusterModel_Gower_PAM")
### * CreateClusterModel_Gower_PAM

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_Gower_PAM
### Title: Fit a projectable Gower-distance PAM model for mixed clinical
###   data
### Aliases: CreateClusterModel_Gower_PAM

### ** Examples




cleanEx()
nameEx("CreateClusterModel_HDBSCAN")
### * CreateClusterModel_HDBSCAN

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_HDBSCAN
### Title: Fit a projectable HDBSCAN model
### Aliases: CreateClusterModel_HDBSCAN

### ** Examples




cleanEx()
nameEx("CreateClusterModel_KMeans")
### * CreateClusterModel_KMeans

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_KMeans
### Title: Fit a projectable K-means clustering model
### Aliases: CreateClusterModel_KMeans

### ** Examples




cleanEx()
nameEx("CreateClusterModel_LatentClass")
### * CreateClusterModel_LatentClass

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_LatentClass
### Title: Fit a projectable latent class model for categorical measures
### Aliases: CreateClusterModel_LatentClass

### ** Examples




cleanEx()
nameEx("CreateClusterModel_MCA_MClust")
### * CreateClusterModel_MCA_MClust

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_MCA_MClust
### Title: Fit MCA followed by Mclust for nominal categorical data
### Aliases: CreateClusterModel_MCA_MClust

### ** Examples




cleanEx()
nameEx("CreateClusterModel_MClust")
### * CreateClusterModel_MClust

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_MClust
### Title: Fit a projectable Gaussian-mixture clustering model
### Aliases: CreateClusterModel_MClust

### ** Examples




cleanEx()
nameEx("CreateClusterModel_PCA_KMeans")
### * CreateClusterModel_PCA_KMeans

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_PCA_KMeans
### Title: Fit PCA followed by K-means
### Aliases: CreateClusterModel_PCA_KMeans

### ** Examples




cleanEx()
nameEx("CreateClusterModel_PCA_MClust")
### * CreateClusterModel_PCA_MClust

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_PCA_MClust
### Title: Fit PCA followed by Mclust
### Aliases: CreateClusterModel_PCA_MClust

### ** Examples




cleanEx()
nameEx("CreateClusterModel_SOM_MClust")
### * CreateClusterModel_SOM_MClust

flush(stderr()); flush(stdout())

### Name: CreateClusterModel_SOM_MClust
### Title: SOM + latent profile clustering pipeline (with AHP and distance
###   baselines)
### Aliases: CreateClusterModel_SOM_MClust Pipeline_SOM_MClust
###   Pipeline_SOMClust CreateSOMClusterModel

### ** Examples




cleanEx()
nameEx("CreateMCAObject")
### * CreateMCAObject

flush(stderr()); flush(stdout())

### Name: CreateMCAObject
### Title: Create a reusable MCA object and visualizations
### Aliases: CreateMCAObject CreateMCATable

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

mca <- CreateMCAObject(
  Labelled,
  VarsToReduce = c("Diagnosis", "Genotype")
)

# Variance explained by MCA dimensions
mca$p_scree

# Variable loadings across MCA dimensions
mca$Lollipop



cleanEx()
nameEx("CreateNormativeTScoreModel")
### * CreateNormativeTScoreModel

flush(stderr()); flush(stdout())

### Name: CreateNormativeTScoreModel
### Title: Create normative T-scores from a regression model
### Aliases: CreateNormativeTScoreModel CreateNormativeTScores

### ** Examples

# A reference sample large enough to estimate the covariate effects
set.seed(4127)
n_Reference <- 240
n_Clinical <- 60
n_Total <- n_Reference + n_Clinical

df <- tibble::tibble(
  Group = c(rep("Reference", n_Reference), rep("Clinical", n_Clinical)),
  Age = round(stats::rnorm(n_Total, mean = 55, sd = 12)),
  Education = factor(sample(
    c("High School", "College", "Graduate"), n_Total, replace = TRUE
  )),
  Sex = factor(sample(c("F", "M"), n_Total, replace = TRUE)),
  Visit = sample(1:3, n_Total, replace = TRUE)
)

# Trail Making A: slower with age, faster with practice, slower if impaired
df$TrailsA <- 1000 * exp(
  3.35 +
    0.011 * (df$Age - 55) -
    0.04 * (df$Visit - 1) +
    0.30 * (df$Group == "Clinical") +
    stats::rnorm(n_Total, sd = 0.16)
)

out <- CreateNormativeTScoreModel(
  data = df,
  test_var = "TrailsA",
  count_var = "Visit",
  covariates = c("Age", "Education", "Sex"),
  reference_var = "Group",
  reference_value = "Reference",
  include_practice_effect = TRUE,
  reverse_score = TRUE,
  convert_seconds = TRUE,
  log_transform = TRUE,
  return_plots = TRUE
)

out$data
out$model

# Raw, transformed, and normed score distributions
out$plots$raw
out$plots$scaled
out$plots$tscore

# T-score diagnostics by reference group, practice count, and covariates
out$plots$reference
out$plots$practice
out$plots$Age
out$plots$Education
out$plots$Sex



cleanEx()
nameEx("CreatePCAObject")
### * CreatePCAObject

flush(stderr()); flush(stdout())

### Name: CreatePCAObject
### Title: Create a reusable PCA object and visualizations
### Aliases: CreatePCAObject CreatePCATable

### ** Examples

PCA <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

PCA_imputed <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3,
  MissingDataStrategy = "impute",
  ImputeMethod = "median",
  MaxMissingForImputation = 0.20
)

PCA_raw_names <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  Relabel = FALSE,
  numComponents = 3
)

PCA_omics <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  Mode = "omics",
  backend = "prcomp",
  maxComponents = 5,
  maxScreeComponents = 5
)

# Display the scree plot from the returned object
PCA$p_scree

# Lollipop-style loading plot of the component loadings
PCA$Lollipop

# Color the lollipop loadings by a variable grouping
PCA_grouped <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  VariableCategories = c("Performance", "Performance", "Size", "Performance",
                         "Drivetrain", "Size", "Performance", "Config",
                         "Config", "Drivetrain", "Config"),
  numComponents = 3
)
PCA_grouped$Lollipop




cleanEx()
nameEx("CreateRCIObject")
### * CreateRCIObject

flush(stderr()); flush(stdout())

### Name: CreateRCIObject
### Title: Create a Reliable Change Index (RCI) object
### Aliases: CreateRCIObject

### ** Examples

set.seed(20260803)
rci_data <- data.frame(
  id = rep(1:30, each = 2),
  visit = rep(c("Baseline", "Followup"), 30),
  Score = round(rnorm(60, mean = 50, sd = 10), 1)
)

# Use a +/-1 cutoff here so all three change classifications are visible.
# The default Confidence = 0.95 retains the conventional +/-1.96 cutoff.
rci <- CreateRCIObject(
  data = rci_data,
  variables = "Score",
  DataFormat = "long",
  id_var = "id",
  VisitColumn = "visit",
  BaselineVisit = "Baseline",
  Confidence = 0.68
)

# Individual trajectories across visits
rci$Plots$Spaghetti$Score

# Participant-level reliable-change values
rci$Plots$Waterfall$Score

# Baseline values relative to reliable change
rci$Plots$Quadrant$Score



cleanEx()
nameEx("CreateStatisticsTable")
### * CreateStatisticsTable

flush(stderr()); flush(stdout())

### Name: CreateStatisticsTable
### Title: Create Statistics Table
### Aliases: CreateStatisticsTable

### ** Examples




cleanEx()
nameEx("CreateSummaryTable")
### * CreateSummaryTable

flush(stderr()); flush(stdout())

### Name: CreateSummaryTable
### Title: Create Summary Table
### Aliases: CreateSummaryTable

### ** Examples




cleanEx()
nameEx("CreateVariableTypesTemplate")
### * CreateVariableTypesTemplate

flush(stderr()); flush(stdout())

### Name: CreateVariableTypesTemplate
### Title: Create a Template for Variable Types
### Aliases: CreateVariableTypesTemplate

### ** Examples





cleanEx()
nameEx("CreateZScoreObject")
### * CreateZScoreObject

flush(stderr()); flush(stdout())

### Name: CreateZScoreObject
### Title: Calculate Z-scores (or standardized scores) and return data +
###   parameters
### Aliases: CreateZScoreObject CalcZScore

### ** Examples





cleanEx()
nameEx("DeriveFreesurferVolumes")
### * DeriveFreesurferVolumes

flush(stderr()); flush(stdout())

### Name: DeriveFreesurferVolumes
### Title: Derive Freesurfer bilateral totals and ICV-adjusted measures
### Aliases: DeriveFreesurferVolumes

### ** Examples

## Not run: 
##D fs_derived <- DeriveFreesurferVolumes(df_freesurfer)
##D 
##D df_freesurfer <- dplyr::bind_cols(
##D   df_freesurfer,
##D   DeriveFreesurferVolumes(df_freesurfer)
##D )
##D 
##D attr(fs_derived, "Freesurfer_derivation_log")
## End(Not run)




cleanEx()
nameEx("ExploreDatasetComparison")
### * ExploreDatasetComparison

flush(stderr()); flush(stdout())

### Name: ExploreDatasetComparison
### Title: Explore dataset comparison results interactively
### Aliases: ExploreDatasetComparison

### ** Examples

data(SampleData)

old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data[-c(1, 2), ]
new_data$MMP7[1:12] <- new_data$MMP7[1:12] * 1.15
new_data$tau[20:35] <- new_data$tau[20:35] + 5
new_data$QualityReview <- ifelse(seq_len(nrow(new_data)) %% 3 == 0, "Review", "Pass")
new_data$Smoker <- NULL
new_data <- rbind(
  new_data,
  transform(new_data[1, ], id = max(old_data$id) + 1, QualityReview = "New record")
)

comparison <- CompareDatasets(old_data, new_data, keys = "id")

# Render the dashboard in an HTML report, Shiny app, Viewer, or HTML file
dashboard <- ExploreDatasetComparison(comparison, TopN = 8)



cleanEx()
nameEx("ExploreMergeValidation")
### * ExploreMergeValidation

flush(stderr()); flush(stdout())

### Name: ExploreMergeValidation
### Title: Explore merge validation results interactively
### Aliases: ExploreMergeValidation

### ** Examples

set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Render the interactive dashboard in an HTML report, Shiny app, Viewer,
# or standalone HTML file.
dashboard <- ExploreMergeValidation(validation)



cleanEx()
nameEx("ExtractPCAComponentSummary")
### * ExtractPCAComponentSummary

flush(stderr()); flush(stdout())

### Name: ExtractPCAComponentSummary
### Title: Extract PCA component summaries
### Aliases: ExtractPCAComponentSummary

### ** Examples





cleanEx()
nameEx("FormattedDataDictionary")
### * FormattedDataDictionary

flush(stderr()); flush(stdout())

### Name: FormattedDataDictionary
### Title: Create a formatted data dictionary table
### Aliases: FormattedDataDictionary

### ** Examples





cleanEx()
nameEx("FreezeTableHeader")
### * FreezeTableHeader

flush(stderr()); flush(stdout())

### Name: FreezeTableHeader
### Title: Freeze the header row of a long table when scrolling
### Aliases: FreezeTableHeader

### ** Examples





cleanEx()
nameEx("IQROutliers")
### * IQROutliers

flush(stderr()); flush(stdout())

### Name: IQROutliers
### Title: Detect outliers using the Tukey IQR rule and visualize results
### Aliases: IQROutliers

### ** Examples

# Synthetic assay data with a few extreme values injected into each batch
set.seed(2024)
lab_data <- data.frame(
  SampleID      = paste0("S", 1:60),
  Batch         = rep(c("Batch1", "Batch2", "Batch3"), each = 20),
  Concentration = c(
    rnorm(19, 100, 10), 180,   # Batch1: one high outlier
    rnorm(19, 105, 10), 15,    # Batch2: one low outlier
    rnorm(18, 98, 10), 175, 190 # Batch3: two high outliers
  )
)

# Flag outliers within each batch
result <- IQROutliers(lab_data, "Concentration",
                      id_var = "SampleID", group = "Batch")
result$outlierdf

# Display the diagnostic plot (flagged outliers shown in red)
result$p

# Without a grouping variable (single combined boxplot)
result_all <- IQROutliers(lab_data, "Concentration",
                          id_var = "SampleID", group = NULL)
result_all$outlierdf




cleanEx()
nameEx("InspectCategoricalSummary")
### * InspectCategoricalSummary

flush(stderr()); flush(stdout())

### Name: InspectCategoricalSummary
### Title: Inspect categorical variables
### Aliases: InspectCategoricalSummary

### ** Examples

df <- data.frame(
  group = factor(c("A", "B", "A", NA), levels = c("A", "B", "C")),
  flag = c(TRUE, FALSE, TRUE, NA),
  text = c("low", "high", "low", "mid")
)

result <- InspectCategoricalSummary(df, Plot = FALSE)
result$Summary

plot_result <- InspectCategoricalSummary(df, variables = "group", Plot = TRUE)
plot_result$Plot



cleanEx()
nameEx("KeepEnv")
### * KeepEnv

flush(stderr()); flush(stdout())

### Name: KeepEnv
### Title: Keep selected objects in an environment and remove everything
###   else
### Aliases: KeepEnv

### ** Examples

# An analysis environment: a few results among many intermediates
env_Analysis <- new.env()
local({
  df_Raw <- data.frame(id = 1:5, value = rnorm(5))
  df_Clean <- df_Raw[!is.na(df_Raw$value), ]
  tmp_merge <- df_Clean
  i <- 3
  scratch_vector <- 1:100
  model_Final <- lm(value ~ id, data = df_Clean)
  df_Results <- data.frame(term = "id", estimate = coef(model_Final)[2])
}, envir = env_Analysis)

ls(env_Analysis)

# Check what would go, before anything is removed
preview <- KeepEnv(
  Keep = c("df_Results", "model_Final"),
  Env = env_Analysis,
  DryRun = TRUE
)
preview$removed

# Then do it for real, and save a workspace holding only what matters
KeepEnv(c("df_Results", "model_Final"), Env = env_Analysis)
ls(env_Analysis)

save(list = ls(env_Analysis), envir = env_Analysis,
     file = file.path(tempdir(), "analysis_results.RData"))

# `Invert = TRUE`: drop the objects named, keep the rest
env_Other <- new.env()
assign("df_Huge", data.frame(x = 1:10), envir = env_Other)
assign("df_Small", data.frame(x = 1:2), envir = env_Other)
KeepEnv("df_Huge", Env = env_Other, Invert = TRUE)
ls(env_Other)




cleanEx()
nameEx("MakeComparisonTable")
### * MakeComparisonTable

flush(stderr()); flush(stdout())

### Name: MakeComparisonTable
### Title: Make comparison table with covariate adjustment, effect sizes,
###   and pairwise contrasts
### Aliases: MakeComparisonTable

### ** Examples





cleanEx()
nameEx("MakeDataDictionary")
### * MakeDataDictionary

flush(stderr()); flush(stdout())

### Name: MakeDataDictionary
### Title: Create a data dictionary for a data frame
### Aliases: MakeDataDictionary Make_DataDictionary

### ** Examples

df <- tibble::tibble(
  group = factor(c("A", "B", "A")),
  status = c("yes", "no", "yes")
)
if (requireNamespace("codebook", quietly = TRUE)) {
  MakeDataDictionary(df)
}



cleanEx()
nameEx("MakeFacetCatComparisonTable")
### * MakeFacetCatComparisonTable

flush(stderr()); flush(stdout())

### Name: MakeFacetCatComparisonTable
### Title: Create a merged gtsummary table by faceting comparisons across
###   multiple categorical variables
### Aliases: MakeFacetCatComparisonTable

### ** Examples





cleanEx()
nameEx("MakePairwiseHeatmap")
### * MakePairwiseHeatmap

flush(stderr()); flush(stdout())

### Name: MakePairwiseHeatmap
### Title: Make a pairwise referent heatmap
### Aliases: MakePairwiseHeatmap

### ** Examples





cleanEx()
nameEx("MakeTable1")
### * MakeTable1

flush(stderr()); flush(stdout())

### Name: MakeTable1
### Title: Create Summary Table using gtsummary
### Aliases: MakeTable1

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels so the table is publication-ready
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A Table 1 across 16 mixed-type variables
vars <- c(
  "Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin",
  "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin", "Apolipoprotein_A1",
  "Apolipoprotein_B", "C_Reactive_Protein", "Cortisol", "Cystatin_C",
  "Ferritin", "Insulin", "Leptin"
)

MakeTable1(Labelled, variables = vars)




cleanEx()
nameEx("MakeUnivariateRegressionTable")
### * MakeUnivariateRegressionTable

flush(stderr()); flush(stdout())

### Name: MakeUnivariateRegressionTable
### Title: Univariate Regression Table
### Aliases: MakeUnivariateRegressionTable UnivariateRegressionTable

### ** Examples





cleanEx()
nameEx("MergeCodebooks")
### * MergeCodebooks

flush(stderr()); flush(stdout())

### Name: MergeCodebooks
### Title: Merge multiple codebooks using harmonization rules
### Aliases: MergeCodebooks

### ** Examples

## Not run: 
##D 
##D MergedCB <- MergeCodebooks(
##D 
##D   codebooks = list(
##D     Study1 = cb1,
##D     Study2 = cb2
##D   ),
##D 
##D   Rules = MergeRules
##D )
##D 
## End(Not run)




cleanEx()
nameEx("Merge_ByClosestTime")
### * Merge_ByClosestTime

flush(stderr()); flush(stdout())

### Name: Merge_ByClosestTime
### Title: Merge Two Data Frames by Closest Time
### Aliases: Merge_ByClosestTime

### ** Examples

# Clinic visits with blood pressure
visits <- data.frame(
  id         = c("A", "B", "C", "D"),
  visit_date = as.Date(c("2024-03-01", "2024-05-20", "2024-02-10",
                         "2024-04-01")),
  sbp        = c(120, 135, 128, 142)
)

# Lab draws (multiple per participant, on different dates)
labs <- data.frame(
  id         = c("A", "A", "B", "B", "C", "D"),
  lab_date   = as.Date(c("2024-01-05", "2024-03-10", "2024-01-20",
                         "2024-06-01", "2024-02-15", "2023-08-15")),
  creatinine = c(0.9, 1.1, 0.8, 1.0, 1.2, 1.4)
)

ShowTable <- function(x, caption = NULL) {
  htmltools::browsable(htmltools::HTML(as.character(
    kableExtra::kable_styling(
      knitr::kable(x, format = "html", caption = caption),
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
  )))
}

# The two tables do not line up
ShowTable(visits, "Clinic visits")
ShowTable(labs, "Lab draws")

# For each visit, attach the lab drawn closest in time (within participant)
res <- Merge_ByClosestTime(
  visits, labs,
  TimeVar1 = "visit_date",
  TimeVar2 = "lab_date",
  keys = "id",
  is_date = TRUE
)

# One row per visit, with the nearest lab attached
ShowTable(res$merged_dataframe, "Visits with nearest lab")

# The gap for each match, which decides whether it is usable
ShowTable(
  data.frame(
    id = res$merged_dataframe$id,
    visit_date = res$merged_dataframe$visit_date,
    DaysToNearestLab = as.numeric(res$time_differences)
  ),
  "Time gap between each visit and its matched lab"
)




cleanEx()
nameEx("MultivariableRegressionTable")
### * MultivariableRegressionTable

flush(stderr()); flush(stdout())

### Name: MultivariableRegressionTable
### Title: Multivariable regression table
### Aliases: MultivariableRegressionTable

### ** Examples




cleanEx()
nameEx("Plot2GroupStats")
### * Plot2GroupStats

flush(stderr()); flush(stdout())

### Name: Plot2GroupStats
### Title: Plot & Summarize Group Stats via MakeComparisonTable (BH q from
###   p; SHAPE by p; COLOR by Category (vector or data frame); stable point
###   size; palette via paletteer)
### Aliases: Plot2GroupStats

### ** Examples




cleanEx()
nameEx("PlotAnovaRelationshipsMatrix")
### * PlotAnovaRelationshipsMatrix

flush(stderr()); flush(stdout())

### Name: PlotAnovaRelationshipsMatrix
### Title: Plot ANOVA Relationships Matrix
### Aliases: PlotAnovaRelationshipsMatrix

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotAnovaRelationshipsMatrix(
  Labelled,
  CatVars = c("Diagnosis", "sex", "Genotype"),
  ContVars = c("age", "ACE_CD143_Angiotensin_Converti",
               "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
               "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
               "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin",
               "Apolipoprotein_A1")
)

# Raw p-value associations
result$Unadjusted$plot

# FDR-adjusted associations
result$FDRCorrected$plot

# The same matrix using Kruskal-Wallis instead of ANOVA
result_NonParametric <- PlotAnovaRelationshipsMatrix(
  Labelled,
  CatVars = c("Diagnosis", "sex", "Genotype"),
  ContVars = c("age", "ACE_CD143_Angiotensin_Converti",
               "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
               "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
               "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin",
               "Apolipoprotein_A1"),
  Parametric = FALSE
)

result_NonParametric$Unadjusted$plot
result_NonParametric$FDRCorrected$plot

# Where the two tests disagree
cols_Key <- c("CategoricalVariable", "ContinuousVariable", "p")
df_Compare <- merge(
  result$Unadjusted$PvalTable[, cols_Key],
  result_NonParametric$Unadjusted$PvalTable[, cols_Key],
  by = c("CategoricalVariable", "ContinuousVariable"),
  suffixes = c("_ANOVA", "_KruskalWallis")
)
df_Compare$AgreesAt05 <-
  (df_Compare$p_ANOVA < 0.05) == (df_Compare$p_KruskalWallis < 0.05)

htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(
    dplyr::mutate(
      df_Compare,
      dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 3))
    ),
    height = "320px", full_width = TRUE
  )
)))

# Covariate adjustment (parametric only)
PlotAnovaRelationshipsMatrix(
  Labelled,
  CatVars = c("Diagnosis", "sex"),
  ContVars = c("AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  covariates = "age"
)$FDRCorrected$plot



cleanEx()
nameEx("PlotAssociations")
### * PlotAssociations

flush(stderr()); flush(stdout())

### Name: PlotAssociations
### Title: Plot Associations
### Aliases: PlotAssociations

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Two categorical variables (grouped bar chart)
PlotAssociations(Labelled, "Diagnosis", "Genotype")

# Two continuous variables (scatter plot with correlation)
PlotAssociations(Labelled, "age", "AXL")

# One continuous and one categorical variable (box/violin plot)
PlotAssociations(Labelled, "Diagnosis", "AXL")



cleanEx()
nameEx("PlotBlandAltman")
### * PlotBlandAltman

flush(stderr()); flush(stdout())

### Name: PlotBlandAltman
### Title: Plot Bland-Altman Agreement Plot
### Aliases: PlotBlandAltman

### ** Examples

# Two devices measuring the same quantity; device B carries a 2-unit bias
set.seed(101)
n <- 80
truth <- rnorm(n, mean = 100, sd = 15)
method_data <- data.frame(
  SampleID = paste0("S", 1:n),
  DeviceA  = truth + rnorm(n, 0, 3),
  DeviceB  = truth + 2 + rnorm(n, 0, 3)
)

result <- PlotBlandAltman(method_data, "DeviceA", "DeviceB")

# Agreement plot: mean difference (bias) and 95% limits of agreement
result$plot

# Bias and limits of agreement
result$stats$mean.diffs
result$stats$lines

# Correlation for the same pair, which the offset does not affect
cor(method_data$DeviceA, method_data$DeviceB)

# A device with no bias but poor precision
method_data$DeviceC <- truth + rnorm(n, 0, 12)
noisy <- PlotBlandAltman(method_data, "DeviceA", "DeviceC")
noisy$plot
noisy$stats$mean.diffs
noisy$stats$lines
cor(method_data$DeviceA, method_data$DeviceC)



cleanEx()
nameEx("PlotCatInteractionEffectsMatrix")
### * PlotCatInteractionEffectsMatrix

flush(stderr()); flush(stdout())

### Name: PlotCatInteractionEffectsMatrix
### Title: Plot Categorical Interaction Effects Matrix
### Aliases: PlotCatInteractionEffectsMatrix

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotCatInteractionEffectsMatrix(
  Labelled,
  predictor_vars = c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                     "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                     "Apolipoprotein_A2", "Apolipoprotein_B",
                     "Apolipoprotein_CI", "Apolipoprotein_CIII",
                     "Apolipoprotein_D", "Apolipoprotein_E"),
  outcome_vars = c("age", "ACE_CD143_Angiotensin_Converti",
                   "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
                   "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
                   "Alpha_1_Microglobulin"),
  interVar = "Diagnosis"
)

# Raw p-value interaction matrix
result$p

# FDR-adjusted interaction matrix
result$p_FDR



cleanEx()
nameEx("PlotCategoricalDistributions")
### * PlotCategoricalDistributions

flush(stderr()); flush(stdout())

### Name: PlotCategoricalDistributions
### Title: Plot categorical distributions
### Aliases: PlotCategoricalDistributions

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

PlotCategoricalDistributions(
  Labelled,
  variables = c("Diagnosis", "Genotype")
)



cleanEx()
nameEx("PlotChiSqCovar")
### * PlotChiSqCovar

flush(stderr()); flush(stdout())

### Name: PlotChiSqCovar
### Title: Plot Chi-Square Tests for Categorical Associations (optionally
###   stratified by covariates)
### Aliases: PlotChiSqCovar

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Derive a few more categorical variables so the matrix has off-diagonal
# structure to read; self-associations are dropped.
Labelled$APOE4 <- ifelse(
  grepl("E4", as.character(Labelled$Genotype)), "Carrier", "Non-carrier")
Labelled$AgeGroup <- cut(
  Labelled$age, breaks = c(-Inf, 65, 80, Inf),
  labels = c("<65", "65-79", "80+"))
Labelled$TauTertile <- cut(
  Labelled$tau,
  breaks = stats::quantile(Labelled$tau, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE),
  labels = c("Low", "Middle", "High"), include.lowest = TRUE)

result <- PlotChiSqCovar(
  Labelled,
  predictor_vars = c("Diagnosis", "sex", "APOE4"),
  outcome_vars = c("Genotype", "AgeGroup", "TauTertile")
)

# Raw p-value associations
result$p

# FDR-adjusted associations
result$p_FDR



cleanEx()
nameEx("PlotClusterBoxplot")
### * PlotClusterBoxplot

flush(stderr()); flush(stdout())

### Name: PlotClusterBoxplot
### Title: Plot cluster boxplots by variable
### Aliases: PlotClusterBoxplot

### ** Examples

PlotClusterBoxplot(
  data = mtcars,
  ClusterVar = "cyl",
  variables = c("mpg", "disp", "hp"),
  Scale = TRUE,
  ReferenceLines = "z"
)

PlotClusterBoxplot(
  data = mtcars,
  ClusterVar = "cyl",
  variables = c("mpg", "disp", "hp"),
  ClusterLabel = "n"
)




cleanEx()
nameEx("PlotClusterCentreHeatmap")
### * PlotClusterCentreHeatmap

flush(stderr()); flush(stdout())

### Name: PlotClusterCentreHeatmap
### Title: Plot cluster centre profiles as a heatmap
### Aliases: PlotClusterCentreHeatmap PlotClusterCentreProfile

### ** Examples




cleanEx()
nameEx("PlotClusterComposition")
### * PlotClusterComposition

flush(stderr()); flush(stdout())

### Name: PlotClusterComposition
### Title: Plot categorical composition by cluster
### Aliases: PlotClusterComposition

### ** Examples




cleanEx()
nameEx("PlotClusterDiagnostic")
### * PlotClusterDiagnostic

flush(stderr()); flush(stdout())

### Name: PlotClusterDiagnostic
### Title: Plot a per-cluster diagnostic value
### Aliases: PlotClusterDiagnostic

### ** Examples




cleanEx()
nameEx("PlotClusterFitReview")
### * PlotClusterFitReview

flush(stderr()); flush(stdout())

### Name: PlotClusterFitReview
### Title: Plot cluster fit-review metrics
### Aliases: PlotClusterFitReview

### ** Examples




cleanEx()
nameEx("PlotClusterMap")
### * PlotClusterMap

flush(stderr()); flush(stdout())

### Name: PlotClusterMap
### Title: Plot a two-dimensional cluster review map
### Aliases: PlotClusterMap PlotClusterAssignment

### ** Examples




cleanEx()
nameEx("PlotClusterProfiles")
### * PlotClusterProfiles

flush(stderr()); flush(stdout())

### Name: PlotClusterProfiles
### Title: Plot labelled numeric profiles by cluster
### Aliases: PlotClusterProfiles

### ** Examples




cleanEx()
nameEx("PlotClusterSilhouette")
### * PlotClusterSilhouette

flush(stderr()); flush(stdout())

### Name: PlotClusterSilhouette
### Title: Plot a per-participant silhouette profile
### Aliases: PlotClusterSilhouette

### ** Examples




cleanEx()
nameEx("PlotContinuousDistributions")
### * PlotContinuousDistributions

flush(stderr()); flush(stdout())

### Name: PlotContinuousDistributions
### Title: Plot Continuous Distributions
### Aliases: PlotContinuousDistributions

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Eight variables across three columns show wrapped labels and multi-row facets.
PlotContinuousDistributions(
  data = Labelled,
  variables = c("AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Ferritin",
                "Gamma_Interferon_induced_Monokin", "MMP7", "tau", "p_tau"),
  ncol = 3
)

# Grouped rain-clouds use the Diagnosis fill to compare distributions.
PlotContinuousDistributions(
  data = Labelled,
  variables = c("Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10", "TRAIL_R3"),
  Fill = "Diagnosis",
  ncol = 3
)



cleanEx()
nameEx("PlotCorrelationsHeatmap")
### * PlotCorrelationsHeatmap

flush(stderr()); flush(stdout())

### Name: PlotCorrelationsHeatmap
### Title: Plot correlations heatmap
### Aliases: PlotCorrelationsHeatmap

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Square heatmap: the same 10 variables on both axes
vars <- c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin",
          "Alpha_2_Macroglobulin", "Apolipoprotein_A1", "Apolipoprotein_B",
          "C_Reactive_Protein", "Cortisol", "Insulin")

square <- PlotCorrelationsHeatmap(
  Labelled,
  predictor_vars = vars,
  outcome_vars = vars
)
square$Unadjusted$plot
square$FDRCorrected$plot

# Rectangular heatmap: different variables on x and y
rectangular <- PlotCorrelationsHeatmap(
  Labelled,
  predictor_vars = c("age", "AXL", "Adiponectin", "Cortisol", "Insulin"),
  outcome_vars = c("Apolipoprotein_A1", "Apolipoprotein_B",
                   "C_Reactive_Protein", "Ferritin", "Leptin")
)
rectangular$Unadjusted$plot



cleanEx()
nameEx("PlotDatasetComparison")
### * PlotDatasetComparison

flush(stderr()); flush(stdout())

### Name: PlotDatasetComparison
### Title: Plot dataset comparison diagnostics
### Aliases: PlotDatasetComparison

### ** Examples

data(SampleData)

# Build two versions of a keyed dataset that differ in records, variables,
# and values, so every diagnostic panel has something to show.
old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)

new_data <- old_data[-(1:8), ]
new_data <- rbind(
  new_data,
  transform(old_data[1:3, ], id = max(old_data$id) + 1:3)
)
new_data$Cohort <- "Wave2"
new_data$Genotype <- NULL
new_data$age[1:25] <- new_data$age[1:25] + 1
new_data$Cortisol[1:15] <- new_data$Cortisol[1:15] * 1.2

comparison <- CompareDatasets(old_data, new_data, keys = "id")

diagnostics <- PlotDatasetComparison(
  comparison,
  Plot = "All",
  interactive = FALSE
)

# Check status, summary metrics, structural, value-change, and top-change views
diagnostics$Checks
diagnostics$SummaryMetrics
diagnostics$StructureChanges
diagnostics$VariableChanges
diagnostics$TopChangedVariables



cleanEx()
nameEx("PlotDirectionalHeatmaps")
### * PlotDirectionalHeatmaps

flush(stderr()); flush(stdout())

### Name: PlotDirectionalHeatmaps
### Title: Create directional heatmaps across continuous & binary variables
### Aliases: PlotDirectionalHeatmaps

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels so binary variables are detected and
# axis labels are readable
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A mix of binary categorical (Diagnosis, sex) and continuous variables
result <- PlotDirectionalHeatmaps(
  Labelled,
  variables = c("Diagnosis", "sex", "age", "AXL", "Adiponectin",
                "Alpha_1_Antitrypsin", "C_Reactive_Protein", "Cortisol",
                "Insulin", "Leptin")
)

# Raw p-value directional heatmap
result$Unadjusted$plot

# FDR-adjusted directional heatmap
result$FDRCorrected$plot

# How binary variables were coded (which level counts as the positive one)
result$BinaryMapping



cleanEx()
nameEx("PlotForestFromTable")
### * PlotForestFromTable

flush(stderr()); flush(stdout())

### Name: PlotForestFromTable
### Title: Create a Forest Plot from Univariate Regression Tables
### Aliases: PlotForestFromTable plotForestFromTable

### ** Examples




cleanEx()
nameEx("PlotInteractionEffectsContinuous")
### * PlotInteractionEffectsContinuous

flush(stderr()); flush(stdout())

### Name: PlotInteractionEffectsContinuous
### Title: Plot Single Interaction Effect
### Aliases: PlotInteractionEffectsContinuous

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axis titles and legend
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Categorical moderator: one regression line per Diagnosis group
PlotInteractionEffectsContinuous(
  Labelled,
  interVar = "Diagnosis",
  outcome_var = "AXL",
  predictor_var = "age"
)

# Continuous moderator: lines at mean and +/- 1 SD of the moderator
PlotInteractionEffectsContinuous(
  Labelled,
  interVar = "Adiponectin",
  outcome_var = "AXL",
  predictor_var = "age"
)



cleanEx()
nameEx("PlotInteractionEffectsMatrix")
### * PlotInteractionEffectsMatrix

flush(stderr()); flush(stdout())

### Name: PlotInteractionEffectsMatrix
### Title: Plot Interaction Effects Matrix
### Aliases: PlotInteractionEffectsMatrix

### ** Examples





cleanEx()
nameEx("PlotMergeValidation")
### * PlotMergeValidation

flush(stderr()); flush(stdout())

### Name: PlotMergeValidation
### Title: Plot merge validation diagnostics
### Aliases: PlotMergeValidation

### ** Examples

set.seed(1)

# `site` comes from both sources and disagrees, leaving a site.x/site.y pair
left <- data.frame(
  id = 1:50,
  site = sample(c("A", "B"), 50, replace = TRUE),
  x = rnorm(50)
)
right <- data.frame(
  id = c(1:45, 101:105),
  site = sample(c("A", "B"), 50, replace = TRUE),
  y = rnorm(50)
)
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

diagnostics <- PlotMergeValidation(
  validation,
  Plot = "All",
  interactive = FALSE
)

# Merge-check status, key coverage, join audit, agreement, and conflicts
diagnostics$Checks
diagnostics$Coverage
diagnostics$JoinAudit
diagnostics$Agreement
diagnostics$Conflicts



cleanEx()
nameEx("PlotMiningMatrix")
### * PlotMiningMatrix

flush(stderr()); flush(stdout())

### Name: PlotMiningMatrix
### Title: PlotMiningMatrix
### Aliases: PlotMiningMatrix

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A mining matrix over 11 mixed-type variables (categorical + continuous)
result <- PlotMiningMatrix(
  Labelled,
  outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
  predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                     "Apolipoprotein_B", "C_Reactive_Protein",
                     "Cortisol", "Insulin")
)

# The relationship plot (point shape/size encode significance from raw p)
result$Unadjusted$plot

# The p-value table carries both unadjusted (p) and FDR-adjusted (p_adj)
# p-values so results can be inspected with and without FDR correction
result$Unadjusted$PvalTable[, c("XVar", "YVar", "p", "p_adj", "Test")]

# Per-outcome FDR correction instead of matrix-wide
result_perout <- PlotMiningMatrix(
  Labelled,
  outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
  predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                     "Apolipoprotein_B", "C_Reactive_Protein",
                     "Cortisol", "Insulin"),
  fdr_scope = "per_outcome"
)

# An interactive version of the matrix
if (requireNamespace("plotly", quietly = TRUE)) {
  plotly::ggplotly(result$Unadjusted$plot)
}



cleanEx()
nameEx("PlotMissingData")
### * PlotMissingData

flush(stderr()); flush(stdout())

### Name: PlotMissingData
### Title: Plot Missing Data
### Aliases: PlotMissingData

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# The revalued data includes missingness defined in the codebook
vars <- c("age", "AXL", "Angiotensinogen", "BMP_6", "IL_6",
          "Fetuin_A", "NT_proBNP", "ENA_78")

PlotMissingData(
  Labelled,
  variables = vars,
  HoverVars = "Diagnosis"
)




cleanEx()
nameEx("PlotNumInteractionEffectsMatrix")
### * PlotNumInteractionEffectsMatrix

flush(stderr()); flush(stdout())

### Name: PlotNumInteractionEffectsMatrix
### Title: Plot Numerical Interaction Effects Matrix
### Aliases: PlotNumInteractionEffectsMatrix

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotNumInteractionEffectsMatrix(
  Labelled,
  predictor_vars = c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                     "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                     "Apolipoprotein_A2", "Apolipoprotein_B",
                     "Apolipoprotein_CI", "Apolipoprotein_CIII",
                     "Apolipoprotein_D", "Apolipoprotein_E"),
  outcome_vars = c("ACE_CD143_Angiotensin_Converti",
                   "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
                   "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
                   "Alpha_1_Microglobulin"),
  interVar = "age"
)

# Raw p-value interaction matrix
result$p

# FDR-adjusted interaction matrix
result$p_FDR



cleanEx()
nameEx("PlotPValueComparisons")
### * PlotPValueComparisons

flush(stderr()); flush(stdout())

### Name: PlotPValueComparisons
### Title: Plot P-Value Comparisons
### Aliases: PlotPValueComparisons

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

PlotPValueComparisons(
  Labelled,
  group_var = "Diagnosis",
  variables = c(
    "age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Cortisol",
    "Ferritin", "GRO_alpha", "MMP10", "MMP7", "NT_proBNP", "PAI_1",
    "TRAIL_R3", "VEGF", "Ab_42", "p_tau", "tau"
  )
)



cleanEx()
nameEx("PlotPartialRegressionScatter")
### * PlotPartialRegressionScatter

flush(stderr()); flush(stdout())

### Name: PlotPartialRegressionScatter
### Title: Partial Regression Plot
### Aliases: PlotPartialRegressionScatter

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

result <- PlotPartialRegressionScatter(
  Labelled,
  IndepVar = "age",
  DepVar = "AXL",
  covariates = "Adiponectin"
)

result$plot



cleanEx()
nameEx("PlotPathway_KT")
### * PlotPathway_KT

flush(stderr()); flush(stdout())

### Name: PlotPathway_KT
### Title: Plot the kynurenine-tryptophan pathway
### Aliases: PlotPathway_KT CreatePathwayPlot_KT

### ** Examples

# A results table keyed by kynurenine-pathway metabolite. Real workflows
# build this with calculate_pathway_results(); here it is entered directly.
results <- data.frame(
  Metabolite = c(
    "Tryptophan", "Serotonin", "N-Formylkynurenine", "Kynurenine",
    "Kynurenic Acid", "3-Hydroxykynurenine", "Anthranilic Acid",
    "Xanthurenic Acid", "3-Hydroxyanthranilic acid", "Quinolinic Acid"
  ),
  correlation = c(0.30, -0.20, 0.50, 0.10, -0.40, 0.60, 0.20, -0.10, 0.30, 0.45),
  p_value = c(0.01, 0.20, 0.03, 0.50, 0.04, 0.001, 0.30, 0.60, 0.02, 0.008),
  p_adj = c(0.05, 0.40, 0.10, 0.70, 0.10, 0.01, 0.50, 0.80, 0.08, 0.03)
)

# Basic usage with raw p-values
PlotPathway_KT(results, "Kynurenine pathway")

# Use FDR-adjusted p-values for significance
PlotPathway_KT(results, "Kynurenine pathway", use_fdr = TRUE)



cleanEx()
nameEx("PlotPhiHeatmap")
### * PlotPhiHeatmap

flush(stderr()); flush(stdout())

### Name: PlotPhiHeatmap
### Title: Plot Phi Correlations Between Binary Variables
### Aliases: PlotPhiHeatmap

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# CatVars must be binary (exactly two unique non-NA values). Derive a few
# more binary indicators so the matrix has off-diagonal structure to read;
# self-associations on the diagonal are masked out.
Labelled$APOE4 <- ifelse(
  grepl("E4", as.character(Labelled$Genotype)), "Carrier", "Non-carrier")
Labelled$HighTau <- ifelse(
  Labelled$tau > stats::median(Labelled$tau, na.rm = TRUE), "High", "Low")
Labelled$LowAbeta <- ifelse(
  Labelled$Ab_42 < stats::median(Labelled$Ab_42, na.rm = TRUE), "Low", "High")

result <- PlotPhiHeatmap(
  Labelled,
  CatVars = c("Diagnosis", "sex", "APOE4", "HighTau", "LowAbeta")
)

# Raw p-value phi heatmap
result$Unadjusted$plot

# FDR-adjusted phi heatmap
result$FDRCorrected$plot



cleanEx()
nameEx("PlotPointCorrelationsHeatmap")
### * PlotPointCorrelationsHeatmap

flush(stderr()); flush(stdout())

### Name: PlotPointCorrelationsHeatmap
### Title: Plot Point-Biserial Correlations Between Binary and Continuous
###   Variables
### Aliases: PlotPointCorrelationsHeatmap

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# CatVars must be binary (exactly two unique non-NA values). These markers
# separate on Diagnosis (Ab_42, tau, p_tau) or on sex (Leptin, PAI_1,
# NT_proBNP), so both rows of the heatmap carry signal in both directions.
result <- PlotPointCorrelationsHeatmap(
  Labelled,
  CatVars = c("Diagnosis", "sex"),
  ContVars = c("Ab_42", "tau", "p_tau", "Leptin", "PAI_1", "NT_proBNP")
)

# Raw p-value point-biserial heatmap
result$Unadjusted$plot

# FDR-adjusted point-biserial heatmap
result$FDRCorrected$plot



cleanEx()
nameEx("PlotSpiderChart")
### * PlotSpiderChart

flush(stderr()); flush(stdout())

### Name: PlotSpiderChart
### Title: Plot a spider chart across continuous and binary variables
### Aliases: PlotSpiderChart

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

vars_biomarkers <- c(
  "Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10", "MMP7", "TRAIL_R3"
)
categories_biomarkers <- c(
  "Neurodegeneration", "Neurodegeneration", "Neurodegeneration",
  "Inflammation", "Inflammation", "Matrix remodeling", "TNF signaling"
)

# Input order preserves the supplied clinical/domain sequence.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "input"
)

# Discrimination puts the largest between-group differences next to each other.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "discrimination"
)

# Hierarchical order groups biomarkers with similar Diagnosis profiles.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "hierarchical"
)

# Greedy order places consecutive spokes with maximally different profiles.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "greedy"
)

# Category/discrimination order keeps domains together, then ranks the
# biomarkers within each domain by their between-group difference.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "category_discrimination",
  VariableCategories = categories_biomarkers
)

# Interactive (plotly) version of the category-aware chart.
PlotSpiderChart(
  data = Labelled,
  variables = vars_biomarkers,
  group_var = "Diagnosis",
  VariableOrder = "category_discrimination",
  VariableCategories = categories_biomarkers,
  interactive = TRUE
)



cleanEx()
nameEx("PlotSplitViolin")
### * PlotSplitViolin

flush(stderr()); flush(stdout())

### Name: PlotSplitViolin
### Title: Split violin with aligned half-boxplots, significance label,
###   sample sizes, and label-aware title
### Aliases: PlotSplitViolin

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Ab_42 has a clear Diagnosis-group difference in the bundled teaching data.
PlotSplitViolin(Labelled, Var = "Ab_42", group_var = "Diagnosis")



cleanEx()
nameEx("PlotSwimmerTransitions")
### * PlotSwimmerTransitions

flush(stderr()); flush(stdout())

### Name: PlotSwimmerTransitions
### Title: Plot swimmer-style transitions for a binary condition over
###   repeated visits
### Aliases: PlotSwimmerTransitions

### ** Examples




cleanEx()
nameEx("PlotTimeDistribution")
### * PlotTimeDistribution

flush(stderr()); flush(stdout())

### Name: PlotTimeDistribution
### Title: Plot Time Distribution
### Aliases: PlotTimeDistribution

### ** Examples

set.seed(1)
df <- data.frame(
  Date = as.Date("2024-01-01") + sample(0:364, 200, replace = TRUE)
)

PlotTimeDistribution(df, DateVariable = "Date")



cleanEx()
nameEx("PlotTimeSwimmer")
### * PlotTimeSwimmer

flush(stderr()); flush(stdout())

### Name: PlotTimeSwimmer
### Title: Plot longitudinal swimmer timelines
### Aliases: PlotTimeSwimmer

### ** Examples

df_swimmer <- tibble::tibble(
  ID = c("P01", "P01", "P01", "P02", "P02", "P03", "P03", "P03", "P04", "P04"),
  Day = c(0, 90, 365, 0, 270, 0, 180, 540, 0, 120),
  State = c("Stable", "Flare", "Recovered", "Stable", "Stable",
            "Flare", "Flare", "Recovered", "Stable", "Withdrawn"),
  Event = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE),
  EventType = c(NA, "Flare", NA, NA, NA, "Flare", NA, NA, NA, "Withdrawal")
)

# State paths, ordered by duration
PlotTimeSwimmer(
  data = df_swimmer,
  id_var = "ID",
  Time = "Day",
  State = "State",
  Event = "Event",
  EventType = "EventType",
  Format = "state_path",
  SortBy = "duration",
  TimeUnit = "months"
)

# Visit points, ordered by last observed follow-up
PlotTimeSwimmer(
  data = df_swimmer,
  id_var = "ID",
  Time = "Day",
  State = "State",
  Format = "visit_points",
  SortBy = "last_time",
  TimeUnit = "months"
)

# Event rugs, centered on each participant's first flare
PlotTimeSwimmer(
  data = subset(df_swimmer, ID %in% c("P01", "P03")),
  id_var = "ID",
  Time = "Day",
  State = "State",
  Event = "Event",
  EventType = "EventType",
  TimeScale = "from_event",
  EventReference = TRUE,
  Format = "event_rug",
  SortBy = "state",
  TimeUnit = "months"
)




cleanEx()
nameEx("PlotVolcanoEffects")
### * PlotVolcanoEffects

flush(stderr()); flush(stdout())

### Name: PlotVolcanoEffects
### Title: Plot volcano-style association effects
### Aliases: PlotVolcanoEffects

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable point labels
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

predictors <- c(
  "ACE_CD143_Angiotensin_Converti", "ACTH_Adrenocorticotropic_Hormon",
  "Adiponectin", "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
  "Alpha_2_Macroglobulin", "Apolipoprotein_A1", "Apolipoprotein_B",
  "B_Lymphocyte_Chemoattractant_BL", "C_Reactive_Protein", "Cortisol",
  "Eotaxin_3", "Ferritin", "Fibrinogen", "GRO_alpha", "IGF_BP_2", "MIF",
  "MMP10", "MMP7", "NT_proBNP", "PAI_1", "Resistin", "TRAIL_R3", "VEGF",
  "Ab_42", "p_tau", "tau"
)

# Continuous outcome, adjusted for age
cont <- PlotVolcanoEffects(
  data = Labelled,
  predictor_vars = predictors,
  outcome_var = "AXL",
  covariates = "age",
  OutcomeType = "continuous",
  LabelMode = "top_n",
  TopN = 3
)
# Raw and FDR-adjusted continuous-outcome volcano plots
cont$RawPPlot
cont$FDRPlot

# Categorical outcome (Diagnosis), Cohen's d effect metric
cat_res <- PlotVolcanoEffects(
  data = Labelled,
  predictor_vars = predictors,
  outcome_var = "Diagnosis",
  OutcomeType = "categorical",
  EffectMetric = "cohens_d",
  LabelMode = "top_n",
  TopN = 3
)
# Raw and FDR-adjusted categorical-outcome volcano plots
cat_res$RawPPlot
cat_res$FDRPlot



cleanEx()
nameEx("PlotZScore")
### * PlotZScore

flush(stderr()); flush(stdout())

### Name: PlotZScore
### Title: Plot Z-score group differences with statistical significance
### Aliases: PlotZScore CreateZScorePlot

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axis labels
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Many variables, grouped into categories
vars <- c("AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
          "Apolipoprotein_A1", "Apolipoprotein_B", "C_Reactive_Protein",
          "Cortisol", "Cystatin_C", "Ferritin", "Insulin", "Leptin")
cats <- c(rep("Metabolic", 4), rep("Lipids", 2),
          rep("Inflammation", 3), rep("Endocrine", 3))

p <- PlotZScore(
  Labelled,
  TargetVar = "Diagnosis",
  variables = vars,
  VariableCategories = cats
)
p

# An interactive version via plotly
if (requireNamespace("plotly", quietly = TRUE)) {
  plotly::ggplotly(p)
}



cleanEx()
nameEx("PrepNumericData")
### * PrepNumericData

flush(stderr()); flush(stdout())

### Name: PrepNumericData
### Title: Prepare numeric data safely for analysis
### Aliases: PrepNumericData

### ** Examples

df <- data.frame(
  x = c("1", "2", "3"),
  y = c(1, 2, Inf),
  z = factor(c("4", "5", "6"))
)

PrepNumericData(df)




cleanEx()
nameEx("ProjectCluster")
### * ProjectCluster

flush(stderr()); flush(stdout())

### Name: ProjectCluster
### Title: Project cases through a fitted clustering model
### Aliases: ProjectCluster Project_SOMClust ProjectSOMCluster

### ** Examples




cleanEx()
nameEx("ProjectRCI")
### * ProjectRCI

flush(stderr()); flush(stdout())

### Name: ProjectRCI
### Title: Project a trained RCI object onto new data
### Aliases: ProjectRCI

### ** Examples

## Not run: 
##D # NOTE: This example is not run because ProjectRCI() currently errors with
##D # "object 'VariableTable' not found" - VariableTable is used but never
##D # defined in the function body (tracked bug, to be fixed separately).
##D rci_data <- data.frame(
##D   id = rep(1:30, each = 2),
##D   visit = rep(c("Baseline", "Followup"), 30),
##D   Score = round(rnorm(60, mean = 50, sd = 10), 1)
##D )
##D 
##D rci <- CreateRCIObject(
##D   data = rci_data,
##D   variables = "Score",
##D   DataFormat = "long",
##D   id_var = "id",
##D   VisitColumn = "visit",
##D   BaselineVisit = "Baseline"
##D )
##D 
##D projected <- ProjectRCI(data = rci_data, Object = rci, id_var = "id")
##D projected$Plots$Spaghetti$Score
## End(Not run)



cleanEx()
nameEx("ReadSciData")
### * ReadSciData

flush(stderr()); flush(stdout())

### Name: ReadSciData
### Title: Read a scientific data file with optional inspection
### Aliases: ReadSciData

### ** Examples





cleanEx()
nameEx("ReplaceMissingCode")
### * ReplaceMissingCode

flush(stderr()); flush(stdout())

### Name: ReplaceMissingCode
### Title: Replace Missing Codes with NA
### Aliases: ReplaceMissingCode

### ** Examples

# `age` uses three different codes to record three different reasons
df <- data.frame(
  id     = 1:6,
  age    = c(34, 999, 52, -7, 41, -8),
  score  = c(10, -9, -9, 15, 20, 12),
  status = c("Active", "Unknown", "Active", "Withdrawn", "Unknown", "Active")
)

# Several markers for one variable go in one cell
codebook <- data.frame(
  Variable    = c("age", "score", "status"),
  MissingCode = c("999, -7, -8", "-9", "Unknown")
)

# Before: the sentinels are averaged in as if they were ages
df
mean(df$age)

# After: every listed code becomes NA
cleaned <- ReplaceMissingCode(df, codebook)
cleaned
mean(cleaned$age, na.rm = TRUE)
colSums(is.na(cleaned))

# Blank codes are skipped
ReplaceMissingCode(
  df,
  data.frame(Variable = c("id", "age"), MissingCode = c(NA, "999, -7, -8"))
)




cleanEx()
nameEx("ReplaceMissingLabels")
### * ReplaceMissingLabels

flush(stderr()); flush(stdout())

### Name: ReplaceMissingLabels
### Title: Replace Missing Labels in Dataframe Columns
### Aliases: ReplaceMissingLabels

### ** Examples

# A data frame where only some columns carry a variable label
df <- data.frame(
  age    = c(52, 61, 77),
  bmi    = c(24.1, 29.7, 22.0),
  smoker = c(0, 1, 0)
)
labelled::var_label(df$age) <- "Age (years)"

# Before: bmi and smoker have no label
sapply(df, function(x) sjlabelled::get_label(x))

# After: unlabelled columns are labelled with their column name
filled <- ReplaceMissingLabels(df)
sapply(filled, function(x) sjlabelled::get_label(x))




cleanEx()
nameEx("RevalueData")
### * RevalueData

flush(stderr()); flush(stdout())

### Name: RevalueData
### Title: Revalue Data
### Aliases: RevalueData

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Before: no labels, and sex is stored as 0/1
sjlabelled::get_label(SampleData$age)   # NULL
class(SampleData$sex)                    # "integer"

# Revalue using the codebook
revalued <- RevalueData(SampleData, SampleVariableTypes)
Labelled <- revalued$RevaluedData

# After: labels attached, sex is a labelled factor
sjlabelled::get_label(Labelled$age)     # "Age"
levels(Labelled$sex)                     # "Female" "Male"

# Recoded variables, and codebook entries not found in the data
revalued$recodedvars
revalued$not_in_data





cleanEx()
nameEx("SampleData")
### * SampleData

flush(stderr()); flush(stdout())

### Name: SampleData
### Title: SampleData for practicing SciDataReportR functions
### Aliases: SampleData
### Keywords: datasets

### ** Examples

data(SampleData)

# 333 participants, 131 columns, two diagnosis groups.
dim(SampleData)
table(SampleData$Diagnosis)

# The first columns are demographics; the rest are biomarkers.
names(SampleData)[1:10]

# As shipped, it is unlabelled and `sex` is still a bare numeric code.
str(SampleData[, c("Diagnosis", "age", "sex", "Genotype", "AXL")])
sjlabelled::get_label(SampleData$age)




cleanEx()
nameEx("SampleVariableTypes")
### * SampleVariableTypes

flush(stderr()); flush(stdout())

### Name: SampleVariableTypes
### Title: Example Dataset: SampleVariableTypes
### Aliases: SampleVariableTypes
### Keywords: datasets

### ** Examples

data(SampleVariableTypes)

dim(SampleVariableTypes)
table(SampleVariableTypes$Type)




cleanEx()
nameEx("SciDataPalette")
### * SciDataPalette

flush(stderr()); flush(stdout())

### Name: SciDataPalette
### Title: SciDataReportR qualitative color palette
### Aliases: SciDataPalette

### ** Examples

SciDataPalette(3)
SciDataPalette(8)
SciDataPalette()




cleanEx()
nameEx("SimulatedPhenotypeData")
### * SimulatedPhenotypeData

flush(stderr()); flush(stdout())

### Name: SimulatedPhenotypeData
### Title: Neutral simulated clustering and phenotyping benchmark
### Aliases: SimulatedPhenotypeData SimulatedPhenotypeVariableTypes

### ** Examples

data(SimulatedPhenotypeData)

# 480 participants, with the truth clusters balanced across both cohorts
dim(SimulatedPhenotypeData)
table(SimulatedPhenotypeData$Cohort, SimulatedPhenotypeData$TruthCluster)





cleanEx()
nameEx("SummarizeTransitions")
### * SummarizeTransitions

flush(stderr()); flush(stdout())

### Name: SummarizeTransitions
### Title: Summarize participant transitions for a binary longitudinal
###   condition
### Aliases: SummarizeTransitions

### ** Examples

toy_df <- tibble::tibble(
  ParticipantID = rep(paste0("P", 1:4), each = 4),
  VisitOrder = rep(1:4, times = 4),
  VisitDate = rep(seq.Date(as.Date("2024-01-01"), by = "month", length.out = 4), times = 4),
  MetSBinary = c(
    0, 0, 1, 1,
    1, 1, 0, 0,
    0, 0, 0, 0,
    TRUE, TRUE, TRUE, TRUE
  )
)

transitions <- SummarizeTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary,
  date_var = VisitDate,
  x_axis_type = "time_from_baseline",
  time_from_baseline_unit = "months"
)

transitions$condition_summary




cleanEx()
nameEx("ValidateMerge")
### * ValidateMerge

flush(stderr()); flush(stdout())

### Name: ValidateMerge
### Title: Validate a merge between two source data frames and a merged
###   result
### Aliases: ValidateMerge

### ** Examples

left <- data.frame(
  record_id = c(1, 1, 2, 2),
  visit_type = c(1, 2, 1, 2),
  age = c(40, 40, 55, 55)
)

right <- data.frame(
  record_id = c(1, 2),
  imaging_score = c(0.4, 0.8)
)

merged <- dplyr::left_join(left, right, by = "record_id")

validation <- ValidateMerge(
  LeftData = left,
  RightData = right,
  MergedData = merged,
  keys = "record_id",
  expected_relationship = "many-to-one"
)

validation$Summary



cleanEx()
nameEx("calculate_pathway_results")
### * calculate_pathway_results

flush(stderr()); flush(stdout())

### Name: calculate_pathway_results
### Title: Calculate Pathway Results for Metabolite Comparisons
### Aliases: calculate_pathway_results

### ** Examples

## Not run: 
##D results <- calculate_pathway_results(
##D   data = my_data,
##D   comparison_var = "poorSleep",
##D   covariates = c("Age", "BMI"),
##D   metabolites = c("Tryptophan", "Kynurenine"),
##D   comparison_type = "binary"
##D )
## End(Not run)



cleanEx()
nameEx("createBinaryMapping")
### * createBinaryMapping

flush(stderr()); flush(stdout())

### Name: createBinaryMapping
### Title: Create a Mapping Table for Binary Variables
### Aliases: createBinaryMapping

### ** Examples





cleanEx()
nameEx("createFacetLabels")
### * createFacetLabels

flush(stderr()); flush(stdout())

### Name: createFacetLabels
### Title: Create facet labels for ggplot2 based on variable labels in a
###   data frame
### Aliases: createFacetLabels

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
vars_Show <- c("AXL", "Adiponectin", "tau", "p_tau")

# Column name on the first line, attached label on the second
labels_Facet <- createFacetLabels(Labelled[vars_Show])
labels_Facet





cleanEx()
nameEx("geom_starcaption")
### * geom_starcaption

flush(stderr()); flush(stdout())

### Name: geom_starcaption
### Title: Add a Caption Explaining Star Annotations
### Aliases: geom_starcaption

### ** Examples

data(SampleData)
data(SampleVariableTypes)

df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Add the caption to a correlation heatmap whose tiles already show stars.
heatmap <- PlotCorrelationsHeatmap(
  data = df_Labelled,
  predictor_vars = c("Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10"),
  outcome_vars = c("MMP7", "TRAIL_R3", "Ferritin", "Fibrinogen", "MIF")
)
heatmap$Unadjusted$plot + geom_starcaption()



cleanEx()
nameEx("getBinaryVars")
### * getBinaryVars

flush(stderr()); flush(stdout())

### Name: getBinaryVars
### Title: Identify Binary Variables
### Aliases: getBinaryVars

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Two-level factors, which are the ones that can be modelled as 0/1
vars_Binary <- getBinaryVars(Labelled)
vars_Binary

# `Revalued = FALSE` looks for any column with two distinct values instead
# of two factor levels, for frames that have not been through RevalueData().
getBinaryVars(SampleData, Revalued = FALSE)

# Which level each one is scored against
createBinaryMapping(Labelled, vars_Binary)




cleanEx()
nameEx("getCatVars")
### * getCatVars

flush(stderr()); flush(stdout())

### Name: getCatVars
### Title: Get Categorical Variables
### Aliases: getCatVars

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Every factor in the frame
getCatVars(Labelled)

# `Ordinal = FALSE` drops ordered factors, which is what you want when the
# ordered variables are going to be analyzed on their numeric scale instead.
getCatVars(Labelled, Ordinal = FALSE)

# Only meaningful after RevalueData(): in the raw extract `sex` is still a
# bare 0/1 numeric column and is not detected as categorical.
getCatVars(SampleData)




cleanEx()
nameEx("getNumVars")
### * getNumVars

flush(stderr()); flush(stdout())

### Name: getNumVars
### Title: Get Numeric Variables
### Aliases: getNumVars

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# 128 numeric columns: age plus the biomarker panel
vars_Numeric <- getNumVars(Labelled)
length(vars_Numeric)
utils::head(vars_Numeric)

# Ordinal variables are excluded by default; include them when they should
# be modelled on their numeric scale.
length(getNumVars(Labelled, Ordinal = TRUE))

# The point of deriving the set: it feeds straight into an analysis without
# a hand-typed vector that can fall out of date.
MakeTable1(Labelled, variables = utils::head(vars_Numeric, 5))




cleanEx()
nameEx("merge_detail")
### * merge_detail

flush(stderr()); flush(stdout())

### Name: merge_detail
### Title: Print a plain-text detail report for one safe_merge result
### Aliases: merge_detail

### ** Examples

baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
merge_detail(m)




cleanEx()
nameEx("merge_summary_table")
### * merge_summary_table

flush(stderr()); flush(stdout())

### Name: merge_summary_table
### Title: Combine safe_merge logs into a single summary table
### Aliases: merge_summary_table

### ** Examples

baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
vitals <- data.frame(id = 1:4, sbp = c(120, 135, 118, 141))
m1 <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
m2 <- safe_merge(m1$data, vitals, by = "id", name = "+ vitals")
merge_summary_table(list(m1$log, m2$log))




cleanEx()
nameEx("plotPCA")
### * plotPCA

flush(stderr()); flush(stdout())

### Name: plotPCA
### Title: Plot PCA scores
### Aliases: plotPCA

### ** Examples





cleanEx()
nameEx("plotSigAssociations")
### * plotSigAssociations

flush(stderr()); flush(stdout())

### Name: plotSigAssociations
### Title: Plot Significant Associations
### Aliases: plotSigAssociations

### ** Examples




cleanEx()
nameEx("plotSigCorrelations")
### * plotSigCorrelations

flush(stderr()); flush(stdout())

### Name: plotSigCorrelations
### Title: Plot Significant Correlations
### Aliases: plotSigCorrelations

### ** Examples




cleanEx()
nameEx("removeString")
### * removeString

flush(stderr()); flush(stdout())

### Name: removeString
### Title: Remove Strings from a Vector
### Aliases: removeString

### ** Examples

data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Every numeric column, including ones that are not really measurements
vars_Numeric <- getNumVars(Labelled)
length(vars_Numeric)

# Drop the ones that should not be screened as biomarkers
vars_Biomarkers <- removeString(vars_Numeric, c("age", "tau", "p_tau"))
length(vars_Biomarkers)
utils::head(vars_Biomarkers, 4)

# Names that are not present are ignored, so an over-broad exclusion list
# is harmless
removeString(c("age", "sex", "AXL"), c("sex", "not_a_column"))

# Unlike setdiff(), duplicates in the original are kept
removeString(c("a", "a", "b", "c"), "b")
setdiff(c("a", "a", "b", "c"), "b")




cleanEx()
nameEx("safe_merge")
### * safe_merge

flush(stderr()); flush(stdout())

### Name: safe_merge
### Title: Safely merge two data frames with relationship-aware validation
### Aliases: safe_merge

### ** Examples

left <- data.frame(
  record_id = c(1, 1, 2, 2),
  visit_type = c(1, 2, 1, 2),
  age = c(40, 40, 55, 55)
)

right <- data.frame(
  record_id = c(1, 2),
  imaging_score = c(0.4, 0.8)
)

m <- safe_merge(
  df_before = left,
  df_add = right,
  by = "record_id",
  name = "Example imaging merge",
  expected_relationship = "many-to-one"
)

m$log




cleanEx()
nameEx("scale_color_SciData")
### * scale_color_SciData

flush(stderr()); flush(stdout())

### Name: scale_color_SciData
### Title: SciDataReportR discrete color scale
### Aliases: scale_color_SciData

### ** Examples

ggplot2::ggplot(
  iris,
  ggplot2::aes(
    x = Sepal.Length,
    y = Sepal.Width,
    color = Species
  )
) +
  ggplot2::geom_point() +
  scale_color_SciData()




cleanEx()
nameEx("scale_color_pvalue")
### * scale_color_pvalue

flush(stderr()); flush(stdout())

### Name: scale_color_pvalue
### Title: Apply an evidence-aware p-value color scale
### Aliases: scale_color_pvalue

### ** Examples





cleanEx()
nameEx("scale_fill_SciData")
### * scale_fill_SciData

flush(stderr()); flush(stdout())

### Name: scale_fill_SciData
### Title: SciDataReportR discrete fill scale
### Aliases: scale_fill_SciData

### ** Examples

ggplot2::ggplot(
  iris,
  ggplot2::aes(
    x = Species,
    y = Sepal.Length,
    fill = Species
  )
) +
  ggplot2::geom_boxplot() +
  scale_fill_SciData()




cleanEx()
nameEx("scale_fill_pvalue")
### * scale_fill_pvalue

flush(stderr()); flush(stdout())

### Name: scale_fill_pvalue
### Title: Apply an evidence-aware p-value fill scale
### Aliases: scale_fill_pvalue

### ** Examples





cleanEx()
nameEx("windsorize")
### * windsorize

flush(stderr()); flush(stdout())

### Name: windsorize
### Title: Winsorize a numeric vector using SD or IQR thresholds
### Aliases: windsorize

### ** Examples

x <- c(rnorm(100), 10, 15, -12)

# SD-based winsorization
windsorize(x, method = "sd", sdlim = 2.5)

# IQR-based winsorization
windsorize(x, method = "iqr", iqrlim = 1.5)

# Compare the distribution before and after winsorization. Both panels are
# drawn on the raw data's x range, because free scales would rescale the
# winsorized panel to its own narrower range and hide the very thing being
# demonstrated.
set.seed(42)
x <- c(rnorm(200, mean = 10, sd = 2), 30, 32, -8, -10)
compare_df <- data.frame(
  raw        = x,
  winsorized = windsorize(x, method = "iqr", iqrlim = 1.5)
)

df_Compare <- tidyr::pivot_longer(
  compare_df,
  cols = dplyr::everything(),
  names_to = "Version",
  values_to = "Value"
)
df_Compare$Version <- factor(
  df_Compare$Version,
  levels = c("raw", "winsorized"),
  labels = c("Raw", "Winsorized")
)

ggplot2::ggplot(df_Compare, ggplot2::aes(x = Value)) +
  ggplot2::geom_histogram(bins = 40, na.rm = TRUE) +
  ggplot2::facet_wrap(~ Version, ncol = 1) +
  ggplot2::coord_cartesian(xlim = range(compare_df$raw, na.rm = TRUE)) +
  ggplot2::labs(
    title = "Winsorization pulls outliers to the limits",
    subtitle = "Both panels share the raw data's x range",
    x = "Value", y = "Count"
  ) +
  ggplot2::theme_bw()

# The four extreme values are gone from the tails and have reappeared as
# taller bars at the winsorized limits; nothing was dropped.
range(compare_df$raw)
range(compare_df$winsorized)
length(compare_df$raw) == length(compare_df$winsorized)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
