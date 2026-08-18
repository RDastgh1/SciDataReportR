# Extracted from test-scale-pvalue.R:89

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
positions <- SciDataReportR:::.warp_pvalues(c(1, 0.5, 0.05))
