# Extracted from test-multivariable-regression-table.R:327

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(508)
df <- data.frame(
    cohort = factor(sample(c("Control", "Long Covid"), 80, replace = TRUE)),
    x1 = rnorm(80),
    x2 = rnorm(80)
  )
df$x1[sample(seq_len(80), 20)] <- NA
df$x2[sample(seq_len(80), 20)] <- NA
res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "cohort",
    PredictorVars = c("x1", "x2"),
    MissingDataStrategy = "complete_cases"
  )
