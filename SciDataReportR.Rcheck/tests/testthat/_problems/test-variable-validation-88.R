# Extracted from test-variable-validation.R:88

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("rstatix")
data("SampleData", package = "SciDataReportR")
data("SampleVariableTypes", package = "SciDataReportR")
df_Test <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
Calls <- list(
    PlotCorrelationsHeatmap = function() PlotCorrelationsHeatmap(
      df_Test, predictor_vars = "age", outcome_vars = "AXL", covariates = "NOPE"),
    PlotChiSqCovar = function() PlotChiSqCovar(
      df_Test, predictor_vars = "sex", outcome_vars = "Diagnosis", covariates = "NOPE"),
    PlotAnovaRelationshipsMatrix = function() PlotAnovaRelationshipsMatrix(
      df_Test, CatVars = "sex", ContVars = "age", covariates = "NOPE"),
    PlotPointCorrelationsHeatmap = function() PlotPointCorrelationsHeatmap(
      df_Test, CatVars = "sex", ContVars = "age", covariates = "NOPE"),
    PlotMiningMatrix = function() PlotMiningMatrix(
      df_Test, outcome_vars = "age", predictor_vars = "AXL", covariates = "NOPE"),
    PlotNumInteractionEffectsMatrix = function() PlotNumInteractionEffectsMatrix(
      df_Test, predictor_vars = "AXL", outcome_vars = "age",
      interVar = "Cortisol", covariates = "NOPE")
  )
