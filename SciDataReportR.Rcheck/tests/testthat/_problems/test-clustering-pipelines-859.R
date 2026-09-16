# Extracted from test-clustering-pipelines.R:859

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("dbscan")
data("SimulatedPhenotypeData")
df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
