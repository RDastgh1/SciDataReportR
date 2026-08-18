# Extracted from test-replace-missing-code.R:102

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(age = c(34, 999, 52))
expect_silent(
    out <- ReplaceMissingCode(df, data.frame(Variable = "age", MissingCode = "Unknown"))
  )
