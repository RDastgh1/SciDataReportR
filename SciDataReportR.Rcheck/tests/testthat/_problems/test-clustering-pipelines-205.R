# Extracted from test-clustering-pipelines.R:205

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("poLCA")
skip_if_not_installed("FactoMineR")
skip_if_not_installed("mclust")
set.seed(103)
df_Categorical <- data.frame(
    symptom = factor(sample(c("none", "mild", "severe"), 60, replace = TRUE)),
    diagnosis = factor(sample(c("control", "case"), 60, replace = TRUE))
  )
lca_model <- CreateClusterModel_LatentClass(df_Categorical, c("symptom", "diagnosis"), final_k = 2, nrep = 2)
