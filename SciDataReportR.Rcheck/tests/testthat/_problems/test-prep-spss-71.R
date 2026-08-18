# Extracted from test-prep-spss.R:71

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
long_name <- paste(rep("question", 40), collapse = " ")
df <- data.frame(x = 1:3)
names(df) <- long_name
res <- PrepSPSS(df, quiet = TRUE)
new_label <- attr(res$data[[1]], "label", exact = TRUE)
expect_true(nchar(new_label, type = "bytes") <= 256)
expect_true(endsWith(new_label, "..."))
