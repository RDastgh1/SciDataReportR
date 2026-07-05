pkgname <- "SciDataReportR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SciDataReportR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SciDataReportR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AddToCodebook")
### * AddToCodebook

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: AddToCodebook
### Title: Add a new variable to a codebook
### Aliases: AddToCodebook

### ** Examples

# Create an empty codebook
codebook <- data.frame(Variable = character(0), Label = character(0),
                       Type = character(0), Category = character(0),
                       Recode = character(0), Code = character(0),
                       Exclude = logical(0), Notes = character(0))

# Add a new variable to the codebook
codebook <- AddToCodebook(codebook, "Age", "Age of participants", "numeric", "Demographics")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("AddToCodebook", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ApplyFDRCorrection")
### * ApplyFDRCorrection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ApplyFDRCorrection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ApplyNormativeTScores")
### * ApplyNormativeTScores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ApplyNormativeTScores
### Title: Apply a normative T-score model to new data
### Aliases: ApplyNormativeTScores

### ** Examples

df <- tibble::tibble(
  Group = c(
    rep("Reference", 8),
    rep("Clinical", 2)
  ),
  Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
  Education = factor(c(
    "College", "College", "Graduate", "Graduate",
    "College", "Graduate", "College", "Graduate",
    "College", "Graduate"
  )),
  Sex = factor(c(
    "F", "M", "F", "M", "F", "M", "F", "M", "F", "M"
  )),
  Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
  TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49) * 1000
)

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

scored_df <- ApplyNormativeTScores(
  data = df,
  normative_obj = norm_obj
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ApplyNormativeTScores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("AssemblePlots")
### * AssemblePlots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: AssemblePlots
### Title: Assemble ggplot objects into a unified multi-panel figure
### Aliases: AssemblePlots

### ** Examples

library(ggplot2)

p1 <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point()

p2 <- ggplot(mtcars, aes(hp, wt)) +
  geom_point()

AssemblePlots(list(p1, p2))

AssemblePlots(
  list(MPG = p1, Horsepower = p2),
  UseNamesAsTitles = TRUE,
  LegendPosition = "top"
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("AssemblePlots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CodebookMergeApp")
### * CodebookMergeApp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CodebookMergeApp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ConvertOrdinalToNumeric")
### * ConvertOrdinalToNumeric

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ConvertOrdinalToNumeric
### Title: Convert ordinal variables to numeric
### Aliases: ConvertOrdinalToNumeric

### ** Examples

# An ordered factor with numeric levels, and one with non-numeric levels
df <- data.frame(
  id     = 1:5,
  likert = factor(c("1", "2", "3", "2", "1"),
                  levels = c("1", "2", "3"), ordered = TRUE),
  grade  = factor(c("A", "B", "A", "C", "B"),
                  levels = c("A", "B", "C"), ordered = TRUE)
)

out <- ConvertOrdinalToNumeric(df)

# likert becomes numeric; grade stays an ordered factor (levels are not numeric)
sapply(out, class)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ConvertOrdinalToNumeric", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateMCAObject")
### * CreateMCAObject

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateMCAObject
### Title: Create a reusable MCA object and visualizations
### Aliases: CreateMCAObject

### ** Examples

data(SampleData)

mca <- CreateMCAObject(
  SampleData,
  VarsToReduce = c("Diagnosis", "Genotype")
)

# Display the scree plot from the returned object
mca$p_scree



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateMCAObject", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateMCATable")
### * CreateMCATable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateMCATable
### Title: Create MCA table and visualization
### Aliases: CreateMCATable

### ** Examples

data(SampleData)

mca <- CreateMCATable(SampleData, VarsToReduce = c("Diagnosis", "Genotype"))
mca$p_scree



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateMCATable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateNormativeTScoreModel")
### * CreateNormativeTScoreModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateNormativeTScoreModel
### Title: Create normative T-scores from a regression model
### Aliases: CreateNormativeTScoreModel

### ** Examples

df <- tibble::tibble(
  Group = c(
    rep("Reference", 8),
    rep("Clinical", 2)
  ),
  Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
  Education = factor(c(
    "College", "College", "Graduate", "Graduate",
    "College", "Graduate", "College", "Graduate",
    "College", "Graduate"
  )),
  Sex = factor(c(
    "F", "M", "F", "M", "F", "M", "F", "M", "F", "M"
  )),
  Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
  TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49) * 1000
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

# Display the fitted T-score plot from the returned object
out$plots$tscore



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateNormativeTScoreModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateNormativeTScores")
### * CreateNormativeTScores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateNormativeTScores
### Title: Create normative T-scores from a regression model
### Aliases: CreateNormativeTScores

### ** Examples

df <- tibble::tibble(
  Group = c(rep("Reference", 8), rep("Clinical", 2)),
  Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
  Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
  TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49)
)

out <- CreateNormativeTScores(
  data = df,
  test_var = "TrailsA",
  count_var = "Visit",
  covariates = "Age",
  reference_var = "Group",
  reference_value = "Reference",
  return_plots = TRUE
)
out$plots$tscore



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateNormativeTScores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreatePCAObject")
### * CreatePCAObject

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreatePCAObject
### Title: Create a reusable PCA object and visualizations
### Aliases: CreatePCAObject

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreatePCAObject", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreatePCATable")
### * CreatePCATable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreatePCATable
### Title: Create PCA table and visualization
### Aliases: CreatePCATable

### ** Examples

PCA <- CreatePCATable(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

# Display the scree plot from the returned object
PCA$p_scree




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreatePCATable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreatePathwayPlot_KT")
### * CreatePathwayPlot_KT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreatePathwayPlot_KT
### Title: Create Kynurenine-Tryptophan Pathway Plot
### Aliases: CreatePathwayPlot_KT

### ** Examples

results <- data.frame(
  Metabolite = c("Tryptophan", "Kynurenine", "Quinolinic Acid"),
  correlation = c(0.3, 0.1, 0.45),
  p_value = c(0.01, 0.5, 0.008),
  p_adj = c(0.05, 0.7, 0.03)
)

CreatePathwayPlot_KT(results, "Kynurenine pathway")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreatePathwayPlot_KT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateRCIObject")
### * CreateRCIObject

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateRCIObject
### Title: Create a Reliable Change Index (RCI) object
### Aliases: CreateRCIObject

### ** Examples

set.seed(1)
rci_data <- data.frame(
  id = rep(1:30, each = 2),
  visit = rep(c("Baseline", "Followup"), 30),
  Score = round(rnorm(60, mean = 50, sd = 10), 1)
)

rci <- CreateRCIObject(
  data = rci_data,
  variables = "Score",
  DataFormat = "long",
  id_var = "id",
  VisitColumn = "visit",
  BaselineVisit = "Baseline"
)

# Display a spaghetti plot from the returned object
rci$Plots$Spaghetti$Score



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateRCIObject", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateSOMClusterModel")
### * CreateSOMClusterModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateSOMClusterModel
### Title: SOM + latent profile clustering pipeline (with AHP and distance
###   baselines)
### Aliases: CreateSOMClusterModel

### ** Examples

## Not run: 
##D # NOTE: This example is kept in \dontrun{} because the function currently
##D # errors with "could not find function 'get_data'": tidyLPA::get_data()
##D # fails to dispatch under this call path (tracked bug, possibly a tidyLPA
##D # version pin issue). The reduced settings below (k_range = 2:4, models = 1)
##D # are otherwise fast enough to run once the bug is resolved.
##D data(SampleData)
##D 
##D model <- CreateSOMClusterModel(
##D   data = SampleData,
##D   variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
##D   method = "exploratory",
##D   k_range = 2:4,
##D   models = 1
##D )
##D model$plots$cluster_fit_summary_plot
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateSOMClusterModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateVariableTypesTemplate")
### * CreateVariableTypesTemplate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateVariableTypesTemplate
### Title: Create a Template for Variable Types
### Aliases: CreateVariableTypesTemplate

### ** Examples

df <- data.frame(
  num = c(1.1, 2.2),
  int = c(1L, 2L),
  fact = factor(c("A", "B")),
  char = c("a", "b"),
  date = as.Date(c("2021-01-01", "2021-01-02"))
)
CreateVariableTypesTemplate(df)
CreateVariableTypesTemplate(df, file.path(tempdir(), "variable_types.csv"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateVariableTypesTemplate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CreateZScorePlot")
### * CreateZScorePlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CreateZScorePlot
### Title: Create a Z-score plot with statistical significance
### Aliases: CreateZScorePlot

### ** Examples

data(SampleData)

CreateZScorePlot(
  SampleData,
  TargetVar = "Diagnosis",
  variables = c("age", "AXL", "Adiponectin")
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CreateZScorePlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ExploreDatasetComparison")
### * ExploreDatasetComparison

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ExploreDatasetComparison
### Title: Explore dataset comparison results interactively
### Aliases: ExploreDatasetComparison

### ** Examples

data(SampleData)

old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data
new_data$age[1:5] <- new_data$age[1:5] + 1

comparison <- CompareDatasets(old_data, new_data, keys = "id")

# Produce the interactive comparison dashboard
ExploreDatasetComparison(comparison)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ExploreDatasetComparison", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ExploreMergeValidation")
### * ExploreMergeValidation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ExploreMergeValidation
### Title: Explore merge validation results interactively
### Aliases: ExploreMergeValidation

### ** Examples

set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Produce the interactive merge-validation dashboard
ExploreMergeValidation(validation)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ExploreMergeValidation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ExtractPCAComponentSummary")
### * ExtractPCAComponentSummary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ExtractPCAComponentSummary
### Title: Extract PCA component summaries
### Aliases: ExtractPCAComponentSummary

### ** Examples

## No test: 
if (requireNamespace("gt", quietly = TRUE)) {
  pca_obj <- CreatePCAObject(
    data = mtcars,
    VarsToReduce = colnames(mtcars)
  )

  summary_obj <- ExtractPCAComponentSummary(pca_obj)

  summary_obj$LongTable
  summary_obj$FormattedSummaryTable
  summary_obj$FormattedSummaryTableLines
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ExtractPCAComponentSummary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FormattedDataDictionary")
### * FormattedDataDictionary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FormattedDataDictionary
### Title: Create a formatted data dictionary table
### Aliases: FormattedDataDictionary

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach variable labels and factor levels first so the dictionary shows
# human-readable labels and correct variable types.
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Formatted data dictionary for a subset of variables
FormattedDataDictionary(
  Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin")]
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FormattedDataDictionary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("IQROutliers")
### * IQROutliers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("IQROutliers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("InspectCategoricalSummary")
### * InspectCategoricalSummary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("InspectCategoricalSummary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MakeComparisonTable")
### * MakeComparisonTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MakeComparisonTable
### Title: Make comparison table with covariate adjustment, effect sizes,
###   and pairwise contrasts
### Aliases: MakeComparisonTable

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels so the table shows readable output
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Compare variables across Diagnosis groups with effect sizes and
# pairwise contrasts
MakeComparisonTable(
  data = Labelled,
  group_var = "Diagnosis",
  variables = c("age", "sex", "AXL", "Adiponectin"),
  AddEffectSize = TRUE,
  AddPairwise = TRUE
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MakeComparisonTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MakeDataDictionary")
### * MakeDataDictionary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MakeDataDictionary
### Title: Create a data dictionary for a data frame
### Aliases: MakeDataDictionary

### ** Examples

df <- tibble::tibble(
  group = factor(c("A", "B", "A")),
  status = c("yes", "no", "yes")
)
if (requireNamespace("codebook", quietly = TRUE)) {
  MakeDataDictionary(df)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MakeDataDictionary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MakeTable1")
### * MakeTable1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MakeTable1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MergeCodebooks")
### * MergeCodebooks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MergeCodebooks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Merge_ByClosestTime")
### * Merge_ByClosestTime

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Merge_ByClosestTime
### Title: Merge Two Data Frames by Closest Time
### Aliases: Merge_ByClosestTime

### ** Examples

# Clinic visits with blood pressure
visits <- data.frame(
  id         = c("A", "B", "C"),
  visit_date = as.Date(c("2024-03-01", "2024-05-20", "2024-02-10")),
  sbp        = c(120, 135, 128)
)

# Lab draws (multiple per participant, on different dates)
labs <- data.frame(
  id         = c("A", "A", "B", "B", "C"),
  lab_date   = as.Date(c("2024-01-05", "2024-03-10", "2024-01-20",
                         "2024-06-01", "2024-02-15")),
  creatinine = c(0.9, 1.1, 0.8, 1.0, 1.2)
)

# For each visit, attach the lab drawn closest in time (within participant)
res <- Merge_ByClosestTime(
  visits, labs,
  TimeVar1 = "visit_date",
  TimeVar2 = "lab_date",
  keys = "id",
  is_date = TRUE
)

res$merged_dataframe
res$time_differences




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Merge_ByClosestTime", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MultivariableRegressionTable")
### * MultivariableRegressionTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MultivariableRegressionTable
### Title: Multivariable regression table
### Aliases: MultivariableRegressionTable

### ** Examples

## No test: 
data(SampleData)

result <- MultivariableRegressionTable(
  SampleData,
  outcome_vars = "AXL",
  predictor_vars = c("Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin"),
  covariates = "age"
)

# Display the regression coefficient matrix plot
result$Plots$RegressionMatrix
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MultivariableRegressionTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Pipeline_SOMClust")
### * Pipeline_SOMClust

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Pipeline_SOMClust
### Title: SOM + latent profile clustering pipeline (with AHP and distance
###   baselines)
### Aliases: Pipeline_SOMClust

### ** Examples

## Not run: 
##D # NOTE: Not run - see CreateSOMClusterModel() for the tracked get_data() bug.
##D data(SampleData)
##D 
##D Pipeline_SOMClust(
##D   data = SampleData,
##D   variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
##D   method = "exploratory",
##D   k_range = 2:4,
##D   models = 1
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Pipeline_SOMClust", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Plot2GroupStats")
### * Plot2GroupStats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Plot2GroupStats
### Title: Plot & Summarize Group Stats via MakeComparisonTable (BH q from
###   p; SHAPE by p; COLOR by Category (vector or data frame); stable point
###   size; palette via paletteer)
### Aliases: Plot2GroupStats

### ** Examples

## No test: 
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable output
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Compare 30 analytes between Diagnosis groups
vars <- c(
  "age", "ACE_CD143_Angiotensin_Converti", "ACTH_Adrenocorticotropic_Hormon",
  "AXL", "Adiponectin", "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
  "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
  "Angiotensinogen", "Apolipoprotein_A_IV", "Apolipoprotein_A1",
  "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
  "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E",
  "Apolipoprotein_H", "B_Lymphocyte_Chemoattractant_BL", "BMP_6",
  "Beta_2_Microglobulin", "Betacellulin", "C_Reactive_Protein", "CD40",
  "CD5L", "Calbindin", "Calcitonin", "CgA"
)

result <- Plot2GroupStats(
  Labelled,
  variables = vars,
  group_var = "Diagnosis",
  impClust = "Impaired",
  normalClust = "Control"
)

# Display the group-comparison plot
result$plot
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Plot2GroupStats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotAnovaRelationshipsMatrix")
### * PlotAnovaRelationshipsMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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

result$Unadjusted$plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotAnovaRelationshipsMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotAssociations")
### * PlotAssociations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotAssociations
### Title: Plot Associations
### Aliases: PlotAssociations

### ** Examples

data(SampleData)

# Two categorical variables (grouped bar chart)
PlotAssociations(SampleData, "Diagnosis", "Genotype")

# Two continuous variables (scatter plot with correlation)
PlotAssociations(SampleData, "age", "AXL")

# One continuous and one categorical variable (box/violin plot)
PlotAssociations(SampleData, "Diagnosis", "AXL")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotAssociations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotBlandAltman")
### * PlotBlandAltman

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotBlandAltman
### Title: Plot Bland-Altman Agreement Plot
### Aliases: PlotBlandAltman

### ** Examples

# Bland-Altman compares two measurements of the SAME quantity on the same
# scale. Here two devices measure the same underlying value, with device B
# carrying a small constant bias plus noise.
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

# Underlying statistics
result$stats$mean.diffs



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotBlandAltman", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCatInteractionEffectsMatrix")
### * PlotCatInteractionEffectsMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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

result$p



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCatInteractionEffectsMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCategoricalDistributions")
### * PlotCategoricalDistributions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotCategoricalDistributions
### Title: Plot categorical distributions
### Aliases: PlotCategoricalDistributions

### ** Examples

data(SampleData)

PlotCategoricalDistributions(
  SampleData,
  variables = c("Diagnosis", "Genotype")
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCategoricalDistributions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotChiSqCovar")
### * PlotChiSqCovar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotChiSqCovar
### Title: Plot Chi-Square Tests for Categorical Associations (optionally
###   stratified by covariates)
### Aliases: PlotChiSqCovar

### ** Examples

data(SampleData)

result <- PlotChiSqCovar(
  SampleData,
  predictor_vars = c("Diagnosis", "Genotype"),
  outcome_vars = c("Diagnosis", "Genotype")
)

result$p



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotChiSqCovar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotClusterBoxplot")
### * PlotClusterBoxplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotClusterBoxplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotContinuousDistributions")
### * PlotContinuousDistributions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotContinuousDistributions
### Title: Plot Continuous Distributions
### Aliases: PlotContinuousDistributions

### ** Examples

data(SampleData)

PlotContinuousDistributions(
  SampleData,
  variables = c("age", "AXL", "Adiponectin")
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotContinuousDistributions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCorrelationsHeatmap")
### * PlotCorrelationsHeatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCorrelationsHeatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotDatasetComparison")
### * PlotDatasetComparison

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotDatasetComparison
### Title: Plot dataset comparison diagnostics
### Aliases: PlotDatasetComparison

### ** Examples

data(SampleData)

# Build two versions of a keyed dataset to compare
old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data
new_data$age[1:5] <- new_data$age[1:5] + 1

comparison <- CompareDatasets(old_data, new_data, keys = "id")

# Display a single diagnostic plot
PlotDatasetComparison(comparison, Plot = "Checks")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotDatasetComparison", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotDirectionalHeatmaps")
### * PlotDirectionalHeatmaps

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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

result$Unadjusted$plot

# How binary variables were coded (which level counts as the positive one)
result$BinaryMapping



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotDirectionalHeatmaps", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotForestFromTable")
### * PlotForestFromTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotForestFromTable
### Title: Create a Forest Plot from Univariate Regression Tables
### Aliases: PlotForestFromTable plotForestFromTable

### ** Examples

## No test: 
data(SampleData)

# Build univariate regression tables to plot
urt <- MakeUnivariateRegressionTable(
  data = SampleData,
  outcome_vars = "AXL",
  predictor_vars = c("age", "Adiponectin", "Alpha_1_Antitrypsin")
)

PlotForestFromTable(urt)

# Or plot a filtered subset of the results
PlotForestFromTable(subset(urt$Results, Predictor != "age"))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotForestFromTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotInteractionEffectsContinuous")
### * PlotInteractionEffectsContinuous

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotInteractionEffectsContinuous", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotInteractionEffectsMatrix")
### * PlotInteractionEffectsMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotInteractionEffectsMatrix
### Title: Plot Interaction Effects Matrix
### Aliases: PlotInteractionEffectsMatrix

### ** Examples

## No test: 
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

outcomes <- c("age", "ACE_CD143_Angiotensin_Converti",
              "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
              "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
              "Alpha_1_Microglobulin")
predictors <- c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
                "Apolipoprotein_A_IV", "Apolipoprotein_A1",
                "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
                "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E")

# With a categorical interaction variable (Diagnosis)
results <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "Diagnosis",
  outcome_vars = outcomes,
  predictor_vars = predictors
)

# With a continuous interaction variable (age)
results_cont <- PlotInteractionEffectsMatrix(
  data = Labelled,
  interVar = "age",
  outcome_vars = setdiff(outcomes, "age"),
  predictor_vars = predictors
)

# Display the plot (colored tiles mark significant interactions)
results$Unadjusted$plot
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotInteractionEffectsMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotMergeValidation")
### * PlotMergeValidation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotMergeValidation
### Title: Plot merge validation diagnostics
### Aliases: PlotMergeValidation

### ** Examples

set.seed(1)
left  <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")

validation <- ValidateMerge(left, right, merged, keys = "id")

# Display a single diagnostic plot
PlotMergeValidation(validation, Plot = "Checks")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotMergeValidation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotMiningMatrix")
### * PlotMiningMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotMiningMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotMissingData")
### * PlotMissingData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotMissingData
### Title: Plot Missing Data
### Aliases: PlotMissingData

### ** Examples

data(SampleData)

# SampleData has real missingness in several assays
vars <- c("age", "AXL", "Angiotensinogen", "BMP_6", "IL_6",
          "Fetuin_A", "NT_proBNP", "ENA_78")

PlotMissingData(
  SampleData,
  variables = vars,
  HoverVars = "Diagnosis"
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotMissingData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotNumInteractionEffectsMatrix")
### * PlotNumInteractionEffectsMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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

result$p



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotNumInteractionEffectsMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPValueComparisons")
### * PlotPValueComparisons

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPValueComparisons
### Title: Plot P-Value Comparisons
### Aliases: PlotPValueComparisons

### ** Examples

data(SampleData)

PlotPValueComparisons(
  SampleData,
  group_var = "Diagnosis",
  variables = c("age", "AXL", "Adiponectin")
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPValueComparisons", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPartialRegressionScatter")
### * PlotPartialRegressionScatter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPartialRegressionScatter
### Title: Partial Regression Plot
### Aliases: PlotPartialRegressionScatter

### ** Examples

data(SampleData)

result <- PlotPartialRegressionScatter(
  SampleData,
  IndepVar = "age",
  DepVar = "AXL",
  covariates = "Adiponectin"
)

result$plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPartialRegressionScatter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPathway_KT")
### * PlotPathway_KT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPathway_KT
### Title: Plot the kynurenine-tryptophan pathway
### Aliases: PlotPathway_KT

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPathway_KT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPhiHeatmap")
### * PlotPhiHeatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPhiHeatmap
### Title: Plot Phi Correlations Between Binary Variables
### Aliases: PlotPhiHeatmap

### ** Examples

data(SampleData)

# CatVars must be binary (exactly two unique non-NA values)
result <- PlotPhiHeatmap(SampleData, CatVars = c("Diagnosis", "sex"))

result$Unadjusted$plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPhiHeatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPointCorrelationsHeatmap")
### * PlotPointCorrelationsHeatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPointCorrelationsHeatmap
### Title: Plot Point-Biserial Correlations Between Binary and Continuous
###   Variables
### Aliases: PlotPointCorrelationsHeatmap

### ** Examples

data(SampleData)

# CatVars must be binary (exactly two unique non-NA values)
result <- PlotPointCorrelationsHeatmap(
  SampleData,
  CatVars = c("Diagnosis", "sex"),
  ContVars = c("age", "AXL", "Adiponectin")
)

result$Unadjusted$plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPointCorrelationsHeatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSpiderChart")
### * PlotSpiderChart

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSpiderChart
### Title: Plot a spider chart across continuous and binary variables
### Aliases: PlotSpiderChart

### ** Examples

data(SampleData)

# Static spider chart
PlotSpiderChart(
  SampleData,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  group_var = "Diagnosis"
)

# Interactive (plotly) version of the same chart
PlotSpiderChart(
  SampleData,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  group_var = "Diagnosis",
  interactive = TRUE
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSpiderChart", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSplitViolin")
### * PlotSplitViolin

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSplitViolin
### Title: Split violin with aligned half-boxplots, significance label,
###   sample sizes, and label-aware title
### Aliases: PlotSplitViolin

### ** Examples

data(SampleData)

PlotSplitViolin(SampleData, Var = "AXL", group_var = "Diagnosis")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSplitViolin", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSwimmerTransitions")
### * PlotSwimmerTransitions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSwimmerTransitions
### Title: Plot swimmer-style transitions for a binary condition over
###   repeated visits
### Aliases: PlotSwimmerTransitions

### ** Examples

## No test: 
# Build a swimmer-shaped dataset from the survival::colon data: each subject
# contributes a variable number of follow-up visits and a binary condition
# (recurrence) that can develop over time.
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(2024)

  subjects <- survival::colon |>
    dplyr::filter(etype == 1) |>
    dplyr::distinct(id, .keep_all = TRUE) |>
    dplyr::slice_sample(n = 40) |>
    dplyr::transmute(
      SubjectID = paste0("P", id),
      n_visits  = pmin(pmax(round(time / 400) + 2, 2), 6),
      onset     = purrr::map_int(n_visits, ~ sample(0:.x, 1))
    )

  swimmer_df <- subjects |>
    dplyr::mutate(Visit = purrr::map(n_visits, seq_len)) |>
    tidyr::unnest(Visit) |>
    dplyr::group_by(SubjectID) |>
    dplyr::mutate(
      VisitDate  = as.Date("2020-01-01") +
        cumsum(c(0, sample(60:200, dplyr::n() - 1, TRUE))),
      Recurrence = as.integer(onset > 0 & Visit >= onset)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(SubjectID, Visit, VisitDate, Recurrence)

  # Aligned by visit number, ordered by when the condition first develops
  PlotSwimmerTransitions(
    data = swimmer_df,
    id_var = SubjectID,
    time_var = Visit,
    status_var = Recurrence,
    order_participants_by = "first_transition"
  )

  # Positioned by elapsed time from each participant's baseline visit
  PlotSwimmerTransitions(
    data = swimmer_df,
    id_var = SubjectID,
    time_var = Visit,
    status_var = Recurrence,
    date_var = VisitDate,
    x_axis_type = "time_from_baseline",
    time_from_baseline_unit = "months"
  )
}
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSwimmerTransitions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotTimeDistribution")
### * PlotTimeDistribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotTimeDistribution
### Title: Plot Time Distribution
### Aliases: PlotTimeDistribution

### ** Examples

set.seed(1)
df <- data.frame(
  Date = as.Date("2024-01-01") + sample(0:364, 200, replace = TRUE)
)

PlotTimeDistribution(df, DateVariable = "Date")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotTimeDistribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotTimeSwimmer")
### * PlotTimeSwimmer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotTimeSwimmer
### Title: Plot longitudinal swimmer timelines
### Aliases: PlotTimeSwimmer

### ** Examples

df <- tibble::tibble(
  ID = c(1,1,1,2,2,2),
  Visit = c(0,6,12,0,6,12),
  Cluster = c("A","A","B","B","B","C")
)

PlotTimeSwimmer(
  data = df,
  id_var = "ID",
  Time = "Visit",
  State = "Cluster"
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotTimeSwimmer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotVolcanoEffects")
### * PlotVolcanoEffects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotVolcanoEffects
### Title: Plot volcano-style association effects
### Aliases: PlotVolcanoEffects

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable point labels
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

predictors <- c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                "Apolipoprotein_A1", "Apolipoprotein_B", "C_Reactive_Protein",
                "Cortisol", "Insulin", "Leptin")

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
cont$RawPPlot

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
cat_res$RawPPlot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotVolcanoEffects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotZScore")
### * PlotZScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotZScore
### Title: Plot Z-score group differences with statistical significance
### Aliases: PlotZScore

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotZScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PrepNumericData")
### * PrepNumericData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PrepNumericData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ProjectRCI")
### * ProjectRCI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ProjectRCI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ProjectSOMCluster")
### * ProjectSOMCluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ProjectSOMCluster
### Title: Project new data onto an existing SOM clinical phenotype space
### Aliases: ProjectSOMCluster

### ** Examples

## Not run: 
##D # NOTE: Not run - projection requires a trained SOM model, and
##D # CreateSOMClusterModel() currently errors on the tracked get_data() bug
##D # (see CreateSOMClusterModel() for details).
##D data(SampleData)
##D 
##D model <- CreateSOMClusterModel(
##D   data = SampleData,
##D   variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
##D   method = "finalize",
##D   final_k = 3,
##D   final_model = 1
##D )
##D 
##D projected <- ProjectSOMCluster(object = model, new_df = SampleData)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ProjectSOMCluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Project_SOMClust")
### * Project_SOMClust

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Project_SOMClust
### Title: Project new data onto an existing SOM clinical phenotype space
### Aliases: Project_SOMClust

### ** Examples

## Not run: 
##D # NOTE: Not run - see ProjectSOMCluster() and CreateSOMClusterModel() for the
##D # tracked get_data() bug that blocks the SOM workflow.
##D Project_SOMClust(object = model, new_df = SampleData)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Project_SOMClust", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReplaceMissingCode")
### * ReplaceMissingCode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReplaceMissingCode
### Title: Replace Missing Codes with NA
### Aliases: ReplaceMissingCode

### ** Examples

# A data frame that uses sentinel values to encode missingness
df <- data.frame(
  id    = 1:6,
  age   = c(34, 999, 52, 999, 41, 29),
  score = c(10, -9, -9, 15, 20, 12)
)

# A codebook mapping each variable to its missing code(s)
codebook <- data.frame(
  Variable    = c("age", "score"),
  MissingCode = c("999", "-9")
)

# Before: sentinel codes still present
df

# After: sentinel codes replaced with NA
ReplaceMissingCode(df, codebook)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReplaceMissingCode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReplaceMissingLabels")
### * ReplaceMissingLabels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReplaceMissingLabels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RevalueData")
### * RevalueData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RevalueData
### Title: Revalue Data
### Aliases: RevalueData

### ** Examples

data(SampleData)
data(SampleVariableTypes)

# Before: no variable labels and coded factors are still raw
# (e.g. sex is stored as 0/1 with no label)
sjlabelled::get_label(SampleData$age)   # NULL
class(SampleData$sex)                    # "integer"

# Revalue using the codebook: attach labels and recode coded factors
revalued <- RevalueData(SampleData, SampleVariableTypes)
Labelled <- revalued$RevaluedData

# After: labels are attached and sex is a labelled factor
sjlabelled::get_label(Labelled$age)     # "Age"
levels(Labelled$sex)                     # "Female" "Male"

# Variables that were recoded, and any codebook entries not found in the data
revalued$recodedvars
revalued$not_in_data




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RevalueData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SampleVariableTypes")
### * SampleVariableTypes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SampleVariableTypes
### Title: Example Dataset: SampleVariableTypes
### Aliases: SampleVariableTypes
### Keywords: datasets

### ** Examples

data(SampleVariableTypes)
head(SampleVariableTypes)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SampleVariableTypes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SummarizeTransitions")
### * SummarizeTransitions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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

SummarizeTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary,
  date_var = VisitDate,
  x_axis_type = "time_from_baseline",
  time_from_baseline_unit = "months"
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SummarizeTransitions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calculate_pathway_results")
### * calculate_pathway_results

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_pathway_results", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("geom_starcaption")
### * geom_starcaption

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: geom_starcaption
### Title: Add a Caption Explaining Star Annotations
### Aliases: geom_starcaption

### ** Examples

library(ggplot2)

# Compose the caption onto any ggplot with `+`
ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  geom_starcaption()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("geom_starcaption", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("merge_detail")
### * merge_detail

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: merge_detail
### Title: Print a plain-text detail report for one safe_merge result
### Aliases: merge_detail

### ** Examples

baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
merge_detail(m)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("merge_detail", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("merge_summary_table")
### * merge_summary_table

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("merge_summary_table", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotPCA")
### * plotPCA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotPCA
### Title: Plot PCA scores
### Aliases: plotPCA

### ** Examples

## No test: 
# Build a PCA object to plot
PCAObj <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

# Default 3D scatter of the first three components
plotPCA(PCAObj)

# 2D scatter colored by a variable in the data
plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  Mode = "2D",
  Var = "cyl",
  ColorType = "factor"
)

# Add hover information
plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  Mode = "2D",
  HoverVars = c("mpg", "hp")
)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotPCA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotSigAssociations")
### * plotSigAssociations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSigAssociations
### Title: Plot Significant Associations
### Aliases: plotSigAssociations

### ** Examples

## No test: 
# Build an ANOVA relationships matrix, then plot the significant pairs
av <- PlotAnovaRelationshipsMatrix(
  mtcars,
  CatVars = c("cyl", "gear"),
  ContVars = c("mpg", "wt", "hp")
)

plots <- plotSigAssociations(mtcars, av)

# Display the first significant-association plot
plots[[1]]
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSigAssociations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotSigCorrelations")
### * plotSigCorrelations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSigCorrelations
### Title: Plot Significant Correlations
### Aliases: plotSigCorrelations

### ** Examples

## No test: 
# Build a correlation heatmap, then plot the significant pairs
ch <- PlotCorrelationsHeatmap(
  mtcars,
  predictor_vars = c("mpg", "wt", "hp"),
  outcome_vars = c("mpg", "wt", "hp")
)

plots <- plotSigCorrelations(mtcars, ch)

# Display the first significant-correlation scatterplot
plots[[1]]
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSigCorrelations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("windsorize")
### * windsorize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: windsorize
### Title: Winsorize a numeric vector using SD or IQR thresholds
### Aliases: windsorize

### ** Examples

x <- c(rnorm(100), 10, 15, -12)

# SD-based winsorization
windsorize(x, method = "sd", sdlim = 2.5)

# IQR-based winsorization
windsorize(x, method = "iqr", iqrlim = 1.5)

# Compare the distribution before and after winsorization
set.seed(42)
x <- c(rnorm(200, mean = 10, sd = 2), 30, 32, -8, -10)
compare_df <- data.frame(
  raw        = x,
  winsorized = windsorize(x, method = "iqr", iqrlim = 1.5)
)

PlotContinuousDistributions(compare_df, variables = c("raw", "winsorized"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("windsorize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
