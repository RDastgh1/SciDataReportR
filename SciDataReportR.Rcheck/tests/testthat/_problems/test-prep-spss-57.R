# Extracted from test-prep-spss.R:57

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
long_label <- strrep("a", 130)
x <- structure(
    c(1, 2, 1),
    labels = stats::setNames(c(1, 2), c(long_label, "ok"))
  )
df <- data.frame(score = 1:3)
df$score <- x
res <- PrepSPSS(df, quiet = TRUE)
new_labels <- names(attr(res$data$score, "labels", exact = TRUE))
expect_true(all(nchar(new_labels, type = "bytes") <= 120))
expect_true("ok" %in% new_labels)
expect_equal(sum(res$label_map$label_type == "value_label"), 1)
expect_equal(res$label_map$original_label[[1]], long_label)
