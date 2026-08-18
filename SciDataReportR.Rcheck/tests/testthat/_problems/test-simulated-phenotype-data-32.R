# Extracted from test-simulated-phenotype-data.R:32

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SimulatedPhenotypeData", package = "SciDataReportR")
df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
model <- CreatePCAKMeansClusterModel(
    df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4,
    nstart = 5)
