# Extracted from test-clustering-pipelines.R:426

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(
    CreateClusterModel_KMeans(
      data.frame(x = 1:5, y = 6:10), c("x", "y"),
      stability_resamples = -1),
    "non-negative integer")
