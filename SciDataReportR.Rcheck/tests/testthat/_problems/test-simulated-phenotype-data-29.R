# Extracted from test-simulated-phenotype-data.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SimulatedPhenotypeData", package = "SciDataReportR")
df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
