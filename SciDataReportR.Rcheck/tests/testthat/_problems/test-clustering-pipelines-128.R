# Extracted from test-clustering-pipelines.R:128

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
DefaultRange <- function(fun, argument) eval(formals(fun)[[argument]])
