# Extracted from test-apply-fdr-correction.R:178

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
diag(m) <- 0
res <- ApplyFDRCorrection(m, fdr_scope = "per_outcome", outcome_margin = 2)
expect_true(all(is.na(diag(res))))
for (j in seq_len(ncol(m))) {
    col_p <- m[, j]
    col_p[j] <- NA
    expect_equal(res[, j], ApplyFDRCorrection(col_p))
  }
