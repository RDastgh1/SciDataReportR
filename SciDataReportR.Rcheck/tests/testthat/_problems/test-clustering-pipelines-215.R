# Extracted from test-clustering-pipelines.R:215

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
skip_if_not_installed("cluster")
data("SimulatedPhenotypeData")
df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:120, ]
