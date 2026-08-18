# Extracted from test-scidata-colors.R:188

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
values <- .ClusterPalette(factor(paste0("C", seq_len(40))))
