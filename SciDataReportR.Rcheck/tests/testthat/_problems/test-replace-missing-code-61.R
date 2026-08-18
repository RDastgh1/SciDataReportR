# Extracted from test-replace-missing-code.R:61

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(id = 1:3, age = c(34, 999, 52))
codebook <- data.frame(
    Variable = c("age", "not_selected"),
    MissingCode = c("999", "-9")
  )
expect_warning(out <- ReplaceMissingCode(df, codebook), "not_selected")
