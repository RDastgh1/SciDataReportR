# Extracted from test-make-pairwise-heatmap.R:156

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("emmeans")
set.seed(101)
df <- data.frame(
    Group = factor(rep(c("Control", "A", "B"), each = 10), levels = c("Control", "A", "B")),
    y1 = c(rnorm(10, 0), rnorm(10, 1.5), rnorm(10, -1.5)),
    y2 = c(rnorm(10, 0), rnorm(10, 1.0), rnorm(10, -1.0)),
    y3 = c(rnorm(10, 0), rnorm(10, 0.5), rnorm(10, -0.5))
  )
res_matrix <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2", "y3"),
    Referent = "Control",
    adjust_scope = "matrix",
    p_adjust_method = "holm"
  )
