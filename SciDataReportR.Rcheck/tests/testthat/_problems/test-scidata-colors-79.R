# Extracted from test-scidata-colors.R:79

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
plain <- .SciDataNamedValues(c("a", "b"))
