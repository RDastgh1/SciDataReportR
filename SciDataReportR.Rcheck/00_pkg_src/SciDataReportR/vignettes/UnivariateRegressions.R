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
  "Calbindin",
  "Ferritin"
)

PredictorVars <- c(
  "Diagnosis",
  "age",
  "sex",
  "Genotype"
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
  covariates = "MMP7"
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
  UniObj$ModelSummaries$Calbindin$Diagnosis
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

