# Extracted from test-clustering-pipelines.R:837

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SimulatedPhenotypeData")
df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
