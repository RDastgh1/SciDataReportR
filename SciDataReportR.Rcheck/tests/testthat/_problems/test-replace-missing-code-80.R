# Extracted from test-replace-missing-code.R:80

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_chr <- data.frame(
    grp = c("A", "Unknown", "B", "Unknown"),
    stringsAsFactors = FALSE
  )
out_chr <- ReplaceMissingCode(
    df_chr, data.frame(Variable = "grp", MissingCode = "Unknown")
  )
expect_equal(out_chr$grp, c("A", NA, "B", NA))
df_fct <- data.frame(grp = factor(c("A", "Refused", "B", "Refused", "C")))
out_fct <- ReplaceMissingCode(
    df_fct, data.frame(Variable = "grp", MissingCode = "Refused")
  )
expect_equal(as.character(out_fct$grp), c("A", NA, "B", NA, "C"))
