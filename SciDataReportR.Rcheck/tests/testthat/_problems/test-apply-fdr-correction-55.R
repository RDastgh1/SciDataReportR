# Extracted from test-apply-fdr-correction.R:55

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
pm <- matrix(c(0.005, 0.04, 0.03, 0.01), nrow = 2)
res <- ApplyFDRCorrection(pm, fdr_scope = "per_predictor")
