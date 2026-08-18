# Extracted from test-make-pairwise-heatmap.R:39

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("emmeans")
df <- data.frame(
    Group = factor(
      rep(c("Control", "A", "B"), each = 12),
      levels = c("Control", "A", "B")
    ),
    y1 = c(seq(-1, 1, length.out = 12),
           seq(1.5, 3.5, length.out = 12),
           seq(-3.5, -1.5, length.out = 12)),
    y2 = c(seq(0, 2, length.out = 12),
           seq(2, 4, length.out = 12),
           seq(4, 6, length.out = 12))
  )
df$y1 <- sjlabelled::set_label(df$y1, "Marker one")
res <- MakePairwiseHeatmap(
    data = df,
    group_var = "Group",
    variables = c("y1", "y2"),
    Referent = "Control",
    Parametric = TRUE,
    adjust_scope = "matrix",
    p_adjust_method = "fdr"
  )
