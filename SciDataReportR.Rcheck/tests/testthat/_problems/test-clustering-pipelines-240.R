# Extracted from test-clustering-pipelines.R:240

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
metrics_same <- .ClusterPartitionMetrics(c(1, 1, 2, 2), c(2, 2, 1, 1))
