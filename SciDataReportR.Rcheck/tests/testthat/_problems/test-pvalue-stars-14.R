# Extracted from test-pvalue-stars.R:14

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
Boundaries <- c(0.0001, 0.001, 0.01, 0.05)
expect_equal(ScidrPValueStars(Boundaries), c("****", "***", "**", "*"))
expect_equal(
    ScidrPValueStars(Boundaries, tiers = ScidrStarTiersThree,
                     ns_label = "", na_label = ""),
    c("***", "***", "**", "*")
  )
