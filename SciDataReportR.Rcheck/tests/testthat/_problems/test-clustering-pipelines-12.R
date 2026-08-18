# Extracted from test-clustering-pipelines.R:12

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
set.seed(101)
df_Test <- data.frame(
    x = c(rnorm(20), rnorm(20, 4), NA_real_),
    y = c(rnorm(20), rnorm(20, 4), NA_real_)
  )
mclust_model <- CreateClusterModel_MClust(
    data = df_Test, variables = c("x", "y"), method = "finalize",
    final_k = 2, final_model = 1
  )
