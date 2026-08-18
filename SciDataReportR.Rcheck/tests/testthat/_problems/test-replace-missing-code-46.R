# Extracted from test-replace-missing-code.R:46

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(id = 1:3, age = c(34, 999, 52))
codebook <- data.frame(
    Variable = c("id", "age"),
    MissingCode = c(NA_character_, NA_character_)
  )
expect_silent(out <- ReplaceMissingCode(df, codebook))
