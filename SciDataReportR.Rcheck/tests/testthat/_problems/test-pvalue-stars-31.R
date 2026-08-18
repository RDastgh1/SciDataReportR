# Extracted from test-pvalue-stars.R:31

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
Probs <- c(0.2, NA)
expect_equal(ScidrPValueStars(Probs), c("ns", NA_character_))
expect_equal(ScidrPValueStars(Probs, ns_label = "", na_label = ""), c("", ""))
