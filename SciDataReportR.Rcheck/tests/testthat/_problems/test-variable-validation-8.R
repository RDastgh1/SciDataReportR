# Extracted from test-variable-validation.R:8

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(a = 1:5, b = 6:10)
expect_error(
    ScidrValidateVariables(df_Test, c("a", "NOPE"), "predictor_vars"),
    "Variables not found in `data`: NOPE \\(supplied to `predictor_vars`\\)",
    fixed = FALSE
  )
