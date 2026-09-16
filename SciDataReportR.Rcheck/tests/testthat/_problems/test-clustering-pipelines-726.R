# Extracted from test-clustering-pipelines.R:726

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:80, ]
