# Extracted from test-replace-missing-code.R:34

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(age = c(34, 999, -7, 41))
codebook <- data.frame(Variable = "age", MissingCode = I(list(c(999, -7))))
out <- ReplaceMissingCode(df, codebook)
