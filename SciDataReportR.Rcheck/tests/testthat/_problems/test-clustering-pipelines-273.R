# Extracted from test-clustering-pipelines.R:273

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
skip_if_not_installed("cluster")
skip_if_not_installed("dbscan")
skip_if_not_installed("poLCA")
skip_if_not_installed("FactoMineR")
data("SimulatedPhenotypeData")
df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
