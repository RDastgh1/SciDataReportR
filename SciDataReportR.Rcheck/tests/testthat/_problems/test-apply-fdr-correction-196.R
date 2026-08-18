# Extracted from test-apply-fdr-correction.R:196

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
MakeSymmetricPmat <- function(p) {
  n <- (1 + sqrt(1 + 8 * length(p))) / 2
  m <- matrix(NA_real_, n, n, dimnames = list(paste0("v", seq_len(n)), paste0("v", seq_len(n))))
  m[lower.tri(m)] <- p
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

# test -------------------------------------------------------------------------
rect <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 2)
expect_error(ApplyFDRCorrection(rect, symmetric = TRUE), "not a symmetric")
