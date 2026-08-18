# Extracted from test-make-pairwise-heatmap.R:10

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 5)),
    y = rnorm(10)
  )
expect_error(
    MakePairwiseHeatmap(df, group_var = "Group", variables = "y"),
    "Referent"
  )
