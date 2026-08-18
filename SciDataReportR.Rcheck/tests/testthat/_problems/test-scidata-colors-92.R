# Extracted from test-scidata-colors.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_equal(
    unname(.SciDataContrastText(c("#0B1F5E", "#641B68"))),
    c("white", "white")
  )
