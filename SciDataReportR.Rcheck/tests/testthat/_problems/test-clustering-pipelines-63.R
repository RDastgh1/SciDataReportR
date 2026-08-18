# Extracted from test-clustering-pipelines.R:63

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
set.seed(112)
df_Test <- data.frame(x = c(rnorm(20), rnorm(20, 4)),
    y = c(rnorm(20), rnorm(20, 4)))
model <- CreateClusterModel_MClust(
    df_Test, c("x", "y"), method = "finalize", final_k = 2,
    final_model = 1, ZScoreType = "None")
