# Extracted from test-variable-validation.R:133

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
Formals <- names(formals(PlotDirectionalHeatmaps))
expect_equal(Formals[1:3], c("data", "variables", "Relabel"))
expect_equal(tail(Formals, 3), c("Data", "xVars", "yVars"))
