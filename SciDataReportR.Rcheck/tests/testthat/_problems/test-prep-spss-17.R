# Extracted from test-prep-spss.R:17

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
long_a <- paste(rep("alpha", 40), collapse = " ")
long_b <- paste(rep("bravo", 40), collapse = " ")
df <- data.frame(group = factor(c(long_a, long_b, "short")))
res <- PrepSPSS(df, quiet = TRUE)
new_levels <- levels(res$data$group)
expect_true(all(nchar(new_levels, type = "bytes") <= 120))
expect_length(unique(new_levels), 3)
expect_true("short" %in% new_levels)
value_rows <- res$label_map[res$label_map$label_type == "value_label", ]
expect_equal(nrow(value_rows), 2)
