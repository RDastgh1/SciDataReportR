# Extracted from test-clustering-pipelines.R:442

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("dbscan")
set.seed(290)
df_Test <- data.frame(x = rnorm(20), y = rnorm(20))
model <- CreateClusterModel_HDBSCAN(
    df_Test, c("x", "y"), method = "finalize", final_minPts = 20,
    final_cluster_selection_epsilon = 0, stability_resamples = 1)
