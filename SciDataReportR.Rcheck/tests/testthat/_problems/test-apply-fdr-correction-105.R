# Extracted from test-apply-fdr-correction.R:105

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
p <- c(0.001, 0.020, 0.040)
m <- MakeSymmetricPmat(p)
res <- ApplyFDRCorrection(m)
expect_equal(res[lower.tri(res)], stats::p.adjust(p, "fdr"))
expect_equal(
    ApplyFDRCorrection(m, method = "bonferroni")[lower.tri(m)],
    stats::p.adjust(p, "bonferroni")
  )
