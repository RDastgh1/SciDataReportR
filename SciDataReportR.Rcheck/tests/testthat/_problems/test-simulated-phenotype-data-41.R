# Extracted from test-simulated-phenotype-data.R:41

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
