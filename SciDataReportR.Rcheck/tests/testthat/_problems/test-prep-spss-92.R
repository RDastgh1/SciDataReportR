# Extracted from test-prep-spss.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
df <- data.frame(
    group = factor(c("A", "B")),
    value = c(1, 2)
  )
res <- PrepSPSS(df, quiet = TRUE)
expect_equal(nrow(res$label_map), 0)
expect_named(
    res$label_map,
    c("spss_name", "label_type", "value", "original_label", "spss_label")
  )
