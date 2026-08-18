# Extracted from test-clustering-pipelines.R:168

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("mclust")
skip_if_not_installed("cluster")
data("SimulatedPhenotypeData")
df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:120, ]
first <- CreateKMeansClusterModel(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
second <- CreateKMeansClusterModel(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
expect_equal(first$Stability$replicates, second$Stability$replicates)
expect_equal(sum(first$ModelInfo$fit_table$Recommended), 1)
expect_true(all(c("ReproducibilityScore", "ahp_index") %in% names(first$ModelInfo$fit_table)))
expect_false("bootstrap_agreement" %in% names(first$Stability$plots))
