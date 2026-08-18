# Extracted from test-prep-spss.R:109

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("janitor")
skip_if_not_installed("haven")
long_level <- paste(rep("chronic condition", 15), collapse = " ")
df <- data.frame(
    diagnosis = factor(c(long_level, "none")),
    value = c(1, 2)
  )
sav_path <- tempfile(fileext = ".sav")
on.exit(unlink(sav_path), add = TRUE)
expect_no_error(PrepSPSS(df, path = sav_path, quiet = TRUE))
expect_true(file.exists(sav_path))
