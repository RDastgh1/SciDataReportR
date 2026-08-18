# Extracted from test-scale-pvalue.R:36

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
label_function <- function(x) paste0("p=", x)
scale <- scale_color_pvalue(
    direction = 1,
    limits = c(1e-6, 0.5),
    breaks = c(0.5, 0.05, 0.001, 1e-6),
    labels = label_function
  )
