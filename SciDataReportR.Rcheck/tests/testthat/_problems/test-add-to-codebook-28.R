# Extracted from test-add-to-codebook.R:28

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_codebook <- data.frame(Variable = "sex", stringsAsFactors = FALSE)
expect_warning(
    actual <- AddToCodebook(df_codebook, "age", Domain = "Clinical"),
    "not an existing codebook column"
  )
