# Extracted from test-scidata-colors.R:178

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
relabelled <- factor(c(
    "1 (n = 40, 33%)", "2 (n = 40, 33%)", "Noise (n = 12, 10%)"
  ))
values <- .ClusterPalette(relabelled)
