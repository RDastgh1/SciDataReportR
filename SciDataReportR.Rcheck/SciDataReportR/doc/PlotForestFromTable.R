## -----------------------------------------------------------------------------
#| label: LoadPackages
library(SciDataReportR)
library(dplyr)


## -----------------------------------------------------------------------------
#| label: LoadData
data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData


## -----------------------------------------------------------------------------
#| label: FitOutcomeScreen
vars_Outcomes <- c(
  "Calbindin",
  "Ferritin",
  "MMP7",
  "Sortilin"
)

vars_Predictors <- c(
  "Diagnosis",
  "age"
)

Reg_Obj_Outcomes <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors,
  Standardize = TRUE
)


## -----------------------------------------------------------------------------
#| label: PlotDefaultForest
PlotForestFromTable(Reg_Obj_Outcomes)


## -----------------------------------------------------------------------------
#| label: PlotFilteredResults
Reg_Obj_Outcomes$Results %>%
  filter(Predictor == "Diagnosis") %>%
  PlotForestFromTable()


## -----------------------------------------------------------------------------
#| label: PlotFlippedForest
PlotForestFromTable(
  Reg_Obj_Outcomes,
  Flip = TRUE
)


## -----------------------------------------------------------------------------
#| label: FitPredictorScreen
Reg_Obj_Predictors <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = "Sortilin",
  predictor_vars = c("Diagnosis", "age", "Calbindin", "Ferritin", "MMP7"),
  Standardize = TRUE
)


## -----------------------------------------------------------------------------
#| label: PlotManyPredictors
PlotForestFromTable(Reg_Obj_Predictors)


## -----------------------------------------------------------------------------
#| label: InspectPlotLabels
Reg_Obj_Outcomes$Metadata$Outcomes


## -----------------------------------------------------------------------------
#| label: FitLogisticScreen
Reg_Obj_Logistic <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = "Diagnosis",
  predictor_vars = c("age", "Calbindin", "Ferritin", "MMP7"),
  Standardize = TRUE
)

Reg_Obj_Logistic$FormattedTable


## -----------------------------------------------------------------------------
#| label: InspectLogisticMetadata
Reg_Obj_Logistic$Metadata$Outcomes


## -----------------------------------------------------------------------------
#| label: PlotLogisticForest
PlotForestFromTable(Reg_Obj_Logistic)


## -----------------------------------------------------------------------------
#| label: Reproducibility
# save.image("PlotForestFromTable_workspace.RData")
print(sessionInfo())

