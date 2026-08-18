# Extracted from test-variable-validation.R:26

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(a = 1:5, b = 6:10)
expect_equal(ScidrValidateVariable(df_Test, "a", "interVar"), "a")
