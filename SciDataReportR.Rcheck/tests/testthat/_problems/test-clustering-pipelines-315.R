# Extracted from test-clustering-pipelines.R:315

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("kohonen")
skip_if_not_installed("aweSOM")
suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
skip_if_not_installed("mclust")
data("SimulatedPhenotypeData")
df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:100, ]
model <- CreateSOMMclustClusterModel(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    final_model = 1, som_xdim = 4, som_ydim = 4,
    stability_resamples = 1, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = NULL)
