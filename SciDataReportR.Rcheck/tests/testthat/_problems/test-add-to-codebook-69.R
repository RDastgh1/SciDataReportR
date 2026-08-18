# Extracted from test-add-to-codebook.R:69

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_codebook <- data.frame(
    Variable = c("sex", "site"),
    Type = c("Categorical", "Categorical"),
    Recode = c(1, 1),
    Exclude = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
expect_warning(
    AddToCodebook(
      df_codebook,
      "age",
      VariableType = "Double"
    ),
    "`Type` has not previously appeared"
  )
