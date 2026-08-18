# Extracted from test-prep-spss.R:34

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
shared_prefix <- strrep("x", 150)
lv1 <- paste0(shared_prefix, " ending one")
lv2 <- paste0(shared_prefix, " ending two")
df <- data.frame(group = factor(c(lv1, lv2)))
res <- PrepSPSS(df, quiet = TRUE)
new_levels <- levels(res$data$group)
expect_true(all(nchar(new_levels, type = "bytes") <= 120))
