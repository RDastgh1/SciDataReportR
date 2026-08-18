# Extracted from test-plot-missing-data.R:36

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(Group = c("A", "B"), Marker = c(1, NA))
expect_error(PlotMissingData(df_Test, x_var = "Unknown"), "x_var was not found")
