# Extracted from test-simulated-phenotype-data.R:4

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SimulatedPhenotypeData", package = "SciDataReportR")
data("SimulatedPhenotypeVariableTypes", package = "SciDataReportR")
expect_equal(nrow(SimulatedPhenotypeData), 480)
