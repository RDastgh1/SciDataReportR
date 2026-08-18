# Extracted from test-plot-missing-data.R:12

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    id = 1:4,
    MarkerOne = c(1, NA, 3, 4),
    MarkerTwo = c(NA, 2, 3, 4)
  )
p <- PlotMissingData(df_Test, variables = c("MarkerOne", "MarkerTwo"), Relabel = FALSE)
expect_s3_class(p, "ggplot")
expect_equal(p$scales$get_scales("x")$name, "Observations")
expect_equal(sort(unique(p$data$.x_plot)), 1:4)
