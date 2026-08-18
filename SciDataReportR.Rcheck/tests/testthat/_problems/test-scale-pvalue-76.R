# Extracted from test-scale-pvalue.R:76

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
p_values <- c(1, 0.7, 0.1, 0.075, 0.05, 0.02, 0.01, 0.001, 1e-5, 1e-8)
warped <- SciDataReportR:::.warp_pvalues(p_values)
