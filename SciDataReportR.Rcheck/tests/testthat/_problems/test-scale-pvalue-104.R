# Extracted from test-scale-pvalue.R:104

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
p_values <- c(1, 0.05, 0.01, 0.001, 1e-8)
inferno <- scale_color_pvalue()
