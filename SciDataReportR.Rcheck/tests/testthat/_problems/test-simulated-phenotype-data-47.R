# Extracted from test-simulated-phenotype-data.R:47

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("poLCA")
skip_if_not_installed("FactoMineR")
skip_if_not_installed("mclust")
data("SimulatedPhenotypeData", package = "SciDataReportR")
df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
df_Projection <- subset(SimulatedPhenotypeData, Cohort == "Projection")
vars_Categorical <- paste0("CatVar", 1:3)
model_LCA <- CreateLatentClassClusterModel(
    df_Training, vars_Categorical, method = "finalize", final_k = 4,
    nrep = 5)
