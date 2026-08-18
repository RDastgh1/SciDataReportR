# Extracted from test-clustering-pipelines.R:403

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
