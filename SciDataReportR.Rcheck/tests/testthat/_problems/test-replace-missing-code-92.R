# Extracted from test-replace-missing-code.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(
    grp = factor(c("A", "Refused", "B"), levels = c("A", "B", "C", "Refused"))
  )
out <- ReplaceMissingCode(df, data.frame(Variable = "grp", MissingCode = "Refused"))
expect_equal(levels(out$grp), c("A", "B", "C"))
