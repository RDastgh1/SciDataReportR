## -----------------------------------------------------------------------------
library(SciDataReportR)
library(dplyr)


## -----------------------------------------------------------------------------
data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData


## -----------------------------------------------------------------------------
OutcomeVars <- c(
  "Diagnosis",
  "sex"
)

PredictorVars <- c(
  "Genotype",
  "Calbindin",
  "Ferritin",
  "MMP7",
  "Calbindin",
  "Sortilin",
  "Osteopontin"
)


## -----------------------------------------------------------------------------
UniObj <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars
)


## -----------------------------------------------------------------------------
UniObj$FormattedTable


## -----------------------------------------------------------------------------
UniObj$LargeTable


## -----------------------------------------------------------------------------
UniObj_Covar <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars,
  covariates = "age"
)


## -----------------------------------------------------------------------------
UniObj_Covar$FormattedTable


## -----------------------------------------------------------------------------
UniObj_Std <- MakeUnivariateRegressionTable(
  data = df_Revalued,
  outcome_vars = OutcomeVars,
  predictor_vars = PredictorVars,
  Standardize = TRUE
)


## -----------------------------------------------------------------------------
UniObj_Std$FormattedTable


## -----------------------------------------------------------------------------
names(
  UniObj$ModelSummaries
)


## -----------------------------------------------------------------------------
summary(
  UniObj$ModelSummaries$Diagnosis$Ferritin
)


## -----------------------------------------------------------------------------
ForestPlot <- PlotForestFromTable(
  UniObj
)

ForestPlot


## -----------------------------------------------------------------------------
PlotForestFromTable(
  UniObj_Covar
)


## -----------------------------------------------------------------------------
sessionInfo()

