# Extracted from test-assemble-plots.R:5

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SampleData", package = "SciDataReportR")
data("SampleVariableTypes", package = "SciDataReportR")
df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
