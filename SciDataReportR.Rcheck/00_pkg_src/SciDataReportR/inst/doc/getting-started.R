## -----------------------------------------------------------------------------
library(SciDataReportR)

data("SampleData")
data("SampleVariableTypes")


## -----------------------------------------------------------------------------
library(readr)

clinical_data <- read_csv("path/to/data.csv")
variable_types <- read_csv("path/to/variable_types.csv")


## -----------------------------------------------------------------------------
variable_types <- CreateVariableTypesTemplate(SampleData)


## -----------------------------------------------------------------------------
data_dictionary <- MakeDataDictionary(SampleData)
FormattedDataDictionary(SampleData)


## -----------------------------------------------------------------------------
revalued <- RevalueData(SampleData, variable_types)
analysis_data <- revalued$RevaluedData


## -----------------------------------------------------------------------------
PlotMissingData(analysis_data, Relabel = TRUE)


## -----------------------------------------------------------------------------
MakeTable1(analysis_data, TreatOrdinalAs = "Continuous")


## -----------------------------------------------------------------------------
continuous_vars <- getNumVars(analysis_data, Ordinal = FALSE)
categorical_vars <- getCatVars(analysis_data)

CreateSummaryTable(analysis_data, continuous_vars, Relabel = TRUE)
PlotContinuousDistributions(analysis_data, continuous_vars[1:12], ncol = 3)
PlotCategoricalDistributions(analysis_data, categorical_vars)


## -----------------------------------------------------------------------------
PlotAssociations(analysis_data, "age", "Adiponectin")
PlotAssociations(analysis_data, "Diagnosis", "Ab_42")


## -----------------------------------------------------------------------------
correlation_result <- PlotCorrelationsHeatmap(
  analysis_data,
  xVars = continuous_vars[1:5],
  yVars = continuous_vars[20:40],
  method = "pearson",
  covars = NULL,
  Relabel = TRUE,
  Ordinal = FALSE
)

correlation_result$Unadjusted$plot


## -----------------------------------------------------------------------------
add_r_and_stars(correlation_result) + geom_starcaption()


## -----------------------------------------------------------------------------
MakeComparisonTable(
  DataFrame = analysis_data,
  CompVariable = "Diagnosis",
  Variables = c("age", "tau", "p_tau"),
  AddEffectSize = TRUE
)


## -----------------------------------------------------------------------------
lab_measures <- continuous_vars[10:60]

PlotZScore(
  analysis_data,
  TargetVar = "Diagnosis",
  Variables = lab_measures,
  sort = FALSE
)


## -----------------------------------------------------------------------------
pca_object <- CreatePCAObject(
  analysis_data,
  lab_measures,
  minThresh = 0.75,
  scale = TRUE,
  center = TRUE
)

pca_object$p_scree
pca_object$Lollipop


## -----------------------------------------------------------------------------
projected_scores <- ProjectPCA(new_data, pca_object)


## -----------------------------------------------------------------------------
sessionInfo()

