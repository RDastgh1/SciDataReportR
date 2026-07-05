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
#| label: DefineLinearVariables
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


## -----------------------------------------------------------------------------
#| label: FitLinearScreen
Reg_Obj_Un <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors,
  Standardize = TRUE
)


## -----------------------------------------------------------------------------
#| label: DisplayFormattedTable
Reg_Obj_Un$FormattedTable


## -----------------------------------------------------------------------------
#| label: DisplayLargeTable
Reg_Obj_Un$LargeTable


## -----------------------------------------------------------------------------
#| label: DisplayResults
Reg_Obj_Un$Results


## -----------------------------------------------------------------------------
#| label: FilterResults
Reg_Obj_Un$Results %>%
  filter(Significant) %>%
  arrange(PValue)


## -----------------------------------------------------------------------------
#| label: FitCovariateScreen
Reg_Obj_Un_Covar <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = vars_Outcomes,
  predictor_vars = vars_Predictors,
  covariates = "sex",
  Standardize = TRUE
)

Reg_Obj_Un_Covar$FormattedTable


## -----------------------------------------------------------------------------
#| label: FitManyPredictors
Reg_Obj_ManyPredictors <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = "Sortilin",
  predictor_vars = c("Diagnosis", "age", "Calbindin", "Ferritin", "MMP7"),
  Standardize = TRUE
)

Reg_Obj_ManyPredictors$FormattedTable


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
#| label: InspectModelObjects
names(Reg_Obj_Un$ModelSummaries)
names(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin)

summary(Reg_Obj_ManyPredictors$ModelSummaries$Sortilin$Ferritin)


## -----------------------------------------------------------------------------
#| label: Reproducibility
# save.image("MakeUnivariateRegressionTable_workspace.RData")
print(sessionInfo())

