# Extracted from test-add-to-codebook.R:41

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(AddToCodebook(list(), "age"), "data frame")
expect_error(AddToCodebook(data.frame(x = 1), "age"), "Variable")
df_codebook <- data.frame(Variable = "age", stringsAsFactors = FALSE)
expect_error(AddToCodebook(df_codebook, NA_character_), "non-missing")
expect_error(AddToCodebook(df_codebook, " "), "non-empty")
