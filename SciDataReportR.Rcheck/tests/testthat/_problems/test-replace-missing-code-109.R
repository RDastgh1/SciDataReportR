# Extracted from test-replace-missing-code.R:109

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(age = c(34, 999))
expect_error(ReplaceMissingCode(df, "not a data frame"), "data frame")
