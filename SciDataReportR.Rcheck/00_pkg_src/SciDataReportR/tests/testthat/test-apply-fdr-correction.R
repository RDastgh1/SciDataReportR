test_that("ApplyFDRCorrection matches hand-computed BH on a vector", {
  # BH by hand for p = c(0.005, 0.04, 0.03, 0.01), n = 4:
  #   sorted:   0.005 (r1) -> 0.02,  0.01 (r2) -> 0.02,
  #             0.03  (r3) -> 0.04,  0.04 (r4) -> 0.04
  p <- c(0.005, 0.04, 0.03, 0.01)
  expect_equal(ApplyFDRCorrection(p), c(0.02, 0.04, 0.04, 0.02))
})

test_that("ApplyFDRCorrection matrix scope adjusts across all cells, keeps NA", {
  pm <- matrix(c(0.005, 0.04, 0.03, 0.01, NA, 0.002),
               nrow = 2,
               dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3")))
  res <- ApplyFDRCorrection(pm, fdr_scope = "matrix")

  # BH by hand over the 5 finite values c(0.005, 0.04, 0.03, 0.01, 0.002):
  #   0.002 -> 0.010, 0.005 -> 0.0125, 0.01 -> 0.05/3, 0.03 -> 0.0375, 0.04 -> 0.04
  expected <- matrix(c(0.0125, 0.04, 0.0375, 0.05 / 3, NA, 0.010),
                     nrow = 2, dimnames = dimnames(pm))
  expect_equal(res, expected)
  expect_identical(dimnames(res), dimnames(pm))
})

test_that("ApplyFDRCorrection per_outcome adjusts within columns (outcome_margin = 2)", {
  pm <- matrix(c(0.005, 0.04, 0.03, 0.01, NA, 0.002),
               nrow = 2,
               dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3")))
  res <- ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 2)

  # BH by hand per column:
  #   out1 c(0.005, 0.04) -> c(0.01, 0.04)
  #   out2 c(0.03, 0.01)  -> c(0.03, 0.02)
  #   out3 c(NA, 0.002)   -> c(NA, 0.002)
  expected <- matrix(c(0.01, 0.04, 0.03, 0.02, NA, 0.002),
                     nrow = 2, dimnames = dimnames(pm))
  expect_equal(res, expected)

  # rows as outcomes must equal adjusting each row separately
  res_rows <- ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 1)
  expect_equal(res_rows[1, ], ApplyFDRCorrection(pm[1, ]))
  expect_equal(res_rows[2, ], ApplyFDRCorrection(pm[2, ]))
})

test_that("ApplyFDRCorrection vector input with outcome_ids groups correctly", {
  p <- c(0.01, 0.04, 0.02, 0.03)
  ids <- c("y1", "y1", "y2", "y2")
  # BH by hand: y1 c(0.01, 0.04) -> c(0.02, 0.04); y2 c(0.02, 0.03) -> c(0.03, 0.03)
  expect_equal(
    ApplyFDRCorrection(p, fdr_scope = "per_outcome", outcome_ids = ids),
    c(0.02, 0.04, 0.03, 0.03)
  )
})

test_that("ApplyFDRCorrection supports per-predictor correction", {
  pm <- matrix(c(0.005, 0.04, 0.03, 0.01), nrow = 2)
  res <- ApplyFDRCorrection(pm, fdr_scope = "per_predictor")
  expect_equal(res[1, ], ApplyFDRCorrection(pm[1, ]))
  expect_equal(res[2, ], ApplyFDRCorrection(pm[2, ]))
  expect_error(ApplyFDRCorrection(c(0.01, 0.02), fdr_scope = "per_predictor"), "predictor_ids")
})

test_that("ApplyFDRCorrection matrix scope reproduces stats::p.adjust exactly", {
  set.seed(7)
  pm <- matrix(runif(24), nrow = 4)
  pm[2, 3] <- NA
  res <- ApplyFDRCorrection(pm)
  pvec <- as.vector(pm)
  expected <- rep(NA_real_, length(pvec))
  ok <- is.finite(pvec)
  expected[ok] <- stats::p.adjust(pvec[ok], method = "fdr")
  expect_equal(as.vector(res), expected)
})

test_that("ApplyFDRCorrection supports data frames and other methods", {
  pdf_in <- data.frame(out1 = c(0.01, 0.02), out2 = c(0.03, 0.04),
                       row.names = c("pred1", "pred2"))
  res <- ApplyFDRCorrection(pdf_in, fdr_scope = "per_outcome", outcome_margin = 2)
  expect_s3_class(res, "data.frame")
  expect_equal(res$out1, stats::p.adjust(pdf_in$out1, "fdr"))
  expect_equal(res$out2, stats::p.adjust(pdf_in$out2, "fdr"))

  p <- c(0.01, 0.02, 0.03)
  expect_equal(ApplyFDRCorrection(p, method = "bonferroni"),
               stats::p.adjust(p, "bonferroni"))
})

MakeSymmetricPmat <- function(p) {
  n <- (1 + sqrt(1 + 8 * length(p))) / 2
  m <- matrix(NA_real_, n, n, dimnames = list(paste0("v", seq_len(n)), paste0("v", seq_len(n))))
  m[lower.tri(m)] <- p
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

test_that("ApplyFDRCorrection corrects each pair of a symmetric matrix once", {
  p <- c(0.001, 0.020, 0.040)
  m <- MakeSymmetricPmat(p)

  res <- ApplyFDRCorrection(m)

  # The family is the three unique pairs, not the six filled cells.
  expect_equal(res[lower.tri(res)], stats::p.adjust(p, "fdr"))
  expect_equal(
    ApplyFDRCorrection(m, method = "bonferroni")[lower.tri(m)],
    stats::p.adjust(p, "bonferroni")
  )
  # Result stays symmetric, diagonal is excluded from the family.
  expect_equal(res, t(res))
  expect_true(all(is.na(diag(res))))
  expect_identical(dimnames(res), dimnames(m))
})

test_that("ApplyFDRCorrection symmetric handling halves the Bonferroni family", {
  p <- c(0.001, 0.020, 0.040)
  m <- MakeSymmetricPmat(p)

  doubled <- ApplyFDRCorrection(m, method = "bonferroni", symmetric = FALSE)
  once <- ApplyFDRCorrection(m, method = "bonferroni", symmetric = TRUE)

  expect_equal(doubled[lower.tri(m)], pmin(1, 2 * once[lower.tri(m)]))
})

test_that("ApplyFDRCorrection keeps a diagonal of self-comparisons out of the family", {
  p <- c(0.001, 0.020, 0.040)
  m <- MakeSymmetricPmat(p)
  diag(m) <- 0

  excluded <- ApplyFDRCorrection(m)
  included <- ApplyFDRCorrection(m, include_diagonal = TRUE)

  expect_true(all(is.na(diag(excluded))))
  expect_equal(excluded[lower.tri(m)], stats::p.adjust(p, "fdr"))

  # Folding perfect self-correlations into the family pulls the off-diagonal
  # values down, which is exactly what excluding the diagonal prevents.
  expect_true(all(included[lower.tri(m)] <= excluded[lower.tri(m)]))
  expect_true(any(included[lower.tri(m)] < excluded[lower.tri(m)]))
  expect_equal(diag(included), rep(0, 3), ignore_attr = TRUE)
})

test_that("ApplyFDRCorrection only treats genuinely symmetric matrices specially", {
  p <- c(0.001, 0.020, 0.040)
  m <- MakeSymmetricPmat(p)

  # Different names down the two margins means these are not the same tests.
  renamed <- m
  colnames(renamed) <- c("a", "b", "c")
  expect_equal(ApplyFDRCorrection(renamed),
               ApplyFDRCorrection(renamed, symmetric = FALSE))

  # One mismatched cell is enough to fall back to whole-matrix correction.
  nearly <- m
  nearly[1, 2] <- 0.5
  expect_equal(ApplyFDRCorrection(nearly),
               ApplyFDRCorrection(nearly, symmetric = FALSE))

  # So is a missing value that only appears in one triangle.
  half_missing <- m
  half_missing[1, 2] <- NA
  expect_equal(ApplyFDRCorrection(half_missing),
               ApplyFDRCorrection(half_missing, symmetric = FALSE))

  rect <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 2)
  expect_equal(ApplyFDRCorrection(rect),
               ApplyFDRCorrection(rect, symmetric = FALSE))
})

test_that("ApplyFDRCorrection per-group scope drops the diagonal of symmetric input", {
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
})

test_that("ApplyFDRCorrection symmetric preserves data frame input", {
  p <- c(0.001, 0.020, 0.040)
  pdf_in <- as.data.frame(MakeSymmetricPmat(p))

  res <- ApplyFDRCorrection(pdf_in)

  expect_s3_class(res, "data.frame")
  expect_identical(names(res), names(pdf_in))
  expect_equal(unname(unlist(res[lower.tri(as.matrix(res))])),
               unname(stats::p.adjust(p, "fdr")))
})

test_that("ApplyFDRCorrection validates symmetric arguments", {
  rect <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 2)
  expect_error(ApplyFDRCorrection(rect, symmetric = TRUE), "not a symmetric")
  expect_error(ApplyFDRCorrection(rect, symmetric = "yes"), "symmetric must be")
  expect_error(ApplyFDRCorrection(rect, include_diagonal = NA), "include_diagonal")
})

test_that("ApplyFDRCorrection validates its inputs", {
  expect_error(ApplyFDRCorrection(c(0.1, 0.2), fdr_scope = "per_outcome"),
               "outcome_ids")
  expect_error(ApplyFDRCorrection(c(0.1, 0.2), fdr_scope = "per_outcome",
                                  outcome_ids = "y1"),
               "same length")
  expect_error(ApplyFDRCorrection(matrix(0.5, 2, 2), outcome_margin = 3),
               "outcome_margin")
  expect_error(ApplyFDRCorrection(matrix("a", 2, 2)), "numeric")
})
