# Extracted from test-label-ordinal-contract.R:87

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("gtsummary")
df_Test <- data.frame(score = c(1, 2, 3))
sjlabelled::set_label(df_Test$score) <- "Clinical score"
labelled <- MakeTable1(df_Test, variables = "score", Relabel = TRUE)
