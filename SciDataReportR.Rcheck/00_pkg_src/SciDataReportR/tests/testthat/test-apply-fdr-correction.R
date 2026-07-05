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
