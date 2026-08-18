# Extracted from test-variable-validation.R:70

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
      df_Test, predictor_vars = c("age", "NOPE"), outcome_vars = "age"),
    PlotAnovaRelationshipsMatrix = function() PlotAnovaRelationshipsMatrix(
      df_Test, CatVars = c("sex", "NOPE"), ContVars = "age"),
    PlotChiSqCovar = function() PlotChiSqCovar(
      df_Test, predictor_vars = c("sex", "NOPE"), outcome_vars = "sex"),
    PlotPointCorrelationsHeatmap = function() PlotPointCorrelationsHeatmap(
      df_Test, CatVars = c("sex", "NOPE"), ContVars = "age"),
    PlotPhiHeatmap = function() PlotPhiHeatmap(
      df_Test, CatVars = c("sex", "NOPE")),
    PlotMiningMatrix = function() PlotMiningMatrix(
      df_Test, outcome_vars = c("age", "NOPE"), predictor_vars = "AXL"),
    PlotDirectionalHeatmaps = function() PlotDirectionalHeatmaps(
      df_Test, variables = c("age", "NOPE")),
    PlotNumInteractionEffectsMatrix = function() PlotNumInteractionEffectsMatrix(
      df_Test, predictor_vars = c("AXL", "NOPE"), outcome_vars = "age",
      interVar = "Cortisol"),
    PlotCatInteractionEffectsMatrix = function() PlotCatInteractionEffectsMatrix(
      df_Test, predictor_vars = c("AXL", "NOPE"), outcome_vars = "age",
      interVar = "sex")
  )
for (nm in names(Calls)) {
    expect_error(Calls[[nm]](), "Variables not found in `data`: NOPE",
                 fixed = TRUE, label = nm)
  }
