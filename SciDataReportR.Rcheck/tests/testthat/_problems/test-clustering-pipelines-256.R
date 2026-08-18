# Extracted from test-clustering-pipelines.R:256

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("dbscan")
skip_if_not_installed("mclust")
data("SimulatedPhenotypeData")
df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:160, ]
