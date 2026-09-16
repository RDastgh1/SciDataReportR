# Extracted from test-multivariable-regression-table.R:300

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(407)
df <- data.frame(
    cohort = factor(sample(c("Control", "Long Covid"), 90, replace = TRUE)),
    x1 = rnorm(90),
    x2 = rnorm(90),
    x_sparse = rnorm(90)
  )
df$x1[sample(seq_len(90), 12)] <- NA
df$x2[sample(seq_len(90), 10)] <- NA
df$x_sparse[sample(seq_len(90), 50)] <- NA
res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "cohort",
    PredictorVars = c("x1", "x2", "x_sparse")
  )
