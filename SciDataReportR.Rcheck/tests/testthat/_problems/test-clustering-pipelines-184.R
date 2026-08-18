# Extracted from test-clustering-pipelines.R:184

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
skip_if_not_installed("dbscan")
set.seed(102)
df_Numeric <- data.frame(a = c(rnorm(20), rnorm(20, 3), NA), b = c(rnorm(20), rnorm(20, 3), NA), c = c(rnorm(20), rnorm(20, 3), NA))
pca_model <- CreateClusterModel_PCA_KMeans(df_Numeric, c("a", "b", "c"), method = "finalize", final_k = 2)
