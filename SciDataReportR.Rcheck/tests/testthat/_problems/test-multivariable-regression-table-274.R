# Extracted from test-multivariable-regression-table.R:274

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
vals <- c(rnorm(100, sd = 0.3), 600)
lims <- SciDataReportR:::ScidrRobustFillLimits(vals)
