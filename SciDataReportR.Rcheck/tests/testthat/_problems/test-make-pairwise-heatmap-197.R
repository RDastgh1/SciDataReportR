# Extracted from test-make-pairwise-heatmap.R:197

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("emmeans")
df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 15), levels = c("Control", "A")),
    y1 = c(seq(-0.2, 0.2, length.out = 15), seq(4.8, 5.2, length.out = 15)),
    y2 = c(seq(-0.2, 0.2, length.out = 15), seq(-5.2, -4.8, length.out = 15))
  )
res <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2"),
    Referent = "Control",
    star_p = "adjusted",
    adjusted_outline = TRUE
  )
