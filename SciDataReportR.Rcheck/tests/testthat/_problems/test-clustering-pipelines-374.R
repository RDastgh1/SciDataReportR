# Extracted from test-clustering-pipelines.R:374

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
data("SimulatedPhenotypeData")
df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
