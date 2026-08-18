# Extracted from test-label-ordinal-contract.R:28

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(severity = ordered(
    c("Mild", "Severe", NA), levels = c("Mild", "Moderate", "Severe")
  ))
sjlabelled::set_label(df_Test$severity) <- "Clinical severity"
out <- ConvertOrdinalToNumeric(df_Test)
expect_equal(as.numeric(out$severity), c(1, 3, NA_real_))
expect_identical(attr(out$severity, "scidr_ordinal_score_source"), "rank")
