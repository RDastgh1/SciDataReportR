# Extracted from test-add-to-codebook.R:77

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
expect_warning(
    AddToCodebook(df_codebook, "age_recode", VariableRecode = 4),
    "`Recode` has not previously appeared"
  )
expect_warning(
    actual <- AddToCodebook(df_codebook, "age_exclude", VariableExclude = 3),
    "`Exclude` has not previously appeared"
  )
