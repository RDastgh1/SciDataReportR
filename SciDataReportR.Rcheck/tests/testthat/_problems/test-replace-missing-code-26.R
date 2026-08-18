# Extracted from test-replace-missing-code.R:26

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(age = c(34, 999, -7, -8, 41))
expected <- c(34, NA, NA, NA, 41)
for (spec in c("999, -7, -8", "999,-7,-8", "999; -7; -8", "999;-7;-8", " 999 , -7 ,-8 ")) {
    out <- ReplaceMissingCode(df, data.frame(Variable = "age", MissingCode = spec))
    expect_equal(out$age, expected, info = spec)
  }
