## -----------------------------------------------------------------------------
#| label: LoadPackages
library(SciDataReportR)


## -----------------------------------------------------------------------------
#| label: LoadData
data("SampleData")
data("SampleVariableTypes")

df_Raw <- SampleData
VariableTypes <- SampleVariableTypes


## -----------------------------------------------------------------------------
#| label: ProjectDataExample
#| eval: false
# df_Clinical <- readr::read_csv(here::here("data", "clinical_data.csv"))
# VariableTypes_Clinical <- readr::read_csv(
#   here::here("data", "variable_types.csv")
# )


## -----------------------------------------------------------------------------
#| label: CreateVariableMetadata
VariableTypes_Template <- CreateVariableTypesTemplate(df_Raw)


## -----------------------------------------------------------------------------
#| label: ReviewDataDictionary
df_DataDictionary <- MakeDataDictionary(df_Raw)
FormattedDataDictionary(df_Raw)


## -----------------------------------------------------------------------------
#| label: RevalueSampleData
df_Revalued <- RevalueData(df_Raw, VariableTypes)$RevaluedData


## -----------------------------------------------------------------------------
#| label: PlotMissingness
PlotMissingData(df_Revalued, Relabel = TRUE)


## -----------------------------------------------------------------------------
#| label: CreateTableOne
MakeTable1(df_Revalued, TreatOrdinalAs = "Continuous")


## -----------------------------------------------------------------------------
#| label: ProfileDistributions
vars_Continuous <- getNumVars(df_Revalued, Ordinal = FALSE)
vars_Categorical <- getCatVars(df_Revalued)

CreateSummaryTable(df_Revalued, vars_Continuous, Relabel = TRUE)
PlotContinuousDistributions(df_Revalued, vars_Continuous[1:12], ncol = 3)
PlotCategoricalDistributions(df_Revalued, vars_Categorical)


## -----------------------------------------------------------------------------
#| label: PlotFocusedAssociations
PlotAssociations(df_Revalued, "age", "Adiponectin")
PlotAssociations(df_Revalued, "Diagnosis", "Ab_42")


## -----------------------------------------------------------------------------
#| label: PlotCorrelationHeatmap
correlation_result <- PlotCorrelationsHeatmap(
  df_Revalued,
  xVars = vars_Continuous[1:5],
  yVars = vars_Continuous[20:40],
  method = "pearson",
  covars = NULL,
  Relabel = TRUE,
  Ordinal = FALSE
)

correlation_result$Unadjusted$plot


## -----------------------------------------------------------------------------
#| label: AnnotateCorrelationHeatmap
add_r_and_stars(correlation_result) + geom_starcaption()


## -----------------------------------------------------------------------------
#| label: CompareDiagnosisGroups
MakeComparisonTable(
  DataFrame = df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c("age", "tau", "p_tau"),
  AddEffectSize = TRUE
)


## -----------------------------------------------------------------------------
#| label: PlotGroupZScores
vars_LabMeasures <- vars_Continuous[10:60]

PlotZScore(
  df_Revalued,
  TargetVar = "Diagnosis",
  Variables = vars_LabMeasures,
  sort = FALSE
)


## -----------------------------------------------------------------------------
#| label: CreatePCAObject
pca_object <- CreatePCAObject(
  df_Revalued,
  vars_LabMeasures,
  minThresh = 0.75,
  scale = TRUE,
  center = TRUE
)

pca_object$p_scree
pca_object$Lollipop


## -----------------------------------------------------------------------------
#| label: ProjectPCAExample
#| eval: false
# projected_scores <- ProjectPCA(new_data, pca_object)


## -----------------------------------------------------------------------------
#| label: Reproducibility
# save.image(here::here("results", "getting-started.RData"))
print(sessionInfo())

