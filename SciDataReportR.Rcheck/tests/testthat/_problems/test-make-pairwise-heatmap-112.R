# Extracted from test-make-pairwise-heatmap.R:112

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("emmeans")
skip_if_not_installed("sandwich")
df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 8), levels = c("Control", "A")),
    y = c(1:8, 6:13)
  )
res <- MakePairwiseHeatmap(
    data = df,
    group_var = "Group",
    variables = "y",
    Referent = "Control",
    Parametric = FALSE,
    p_adjust_method = "none"
  )
