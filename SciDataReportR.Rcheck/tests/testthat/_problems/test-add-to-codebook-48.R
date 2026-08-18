# Extracted from test-add-to-codebook.R:48

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_codebook <- data.frame(Variable = "age", stringsAsFactors = FALSE)
expect_error(AddToCodebook(df_codebook, "sex", "a", Domain = c("A", "B")), "single atomic")
