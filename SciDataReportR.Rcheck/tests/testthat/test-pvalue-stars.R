# The package used to carry six separate star implementations that disagreed on
# whether the cut points were inclusive. A p of exactly 0.001 earned "***" in
# MultivariableRegressionTable() and PlotInteractionEffectsMatrix() but only
# "**" in PlotSplitViolin(), while every caption stated the strict "<" form.

test_that("star thresholds are inclusive upper bounds", {
  Boundaries <- c(0.0001, 0.001, 0.01, 0.05)

  expect_equal(ScidrPValueStars(Boundaries), c("****", "***", "**", "*"))
  expect_equal(
    ScidrPValueStars(Boundaries, tiers = ScidrStarTiersThree,
                     ns_label = "", na_label = ""),
    c("***", "***", "**", "*")
  )
})

test_that("star scale matches rstatix::add_significance", {
  skip_if_not_installed("rstatix")
  Probs <- c(0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 0.06, 0.5)

  expect_equal(
    ScidrPValueStars(Probs),
    rstatix::add_significance(data.frame(p = Probs), "p")$p.signif
  )
})

test_that("ns and NA labels are configurable and NA never becomes a star", {
  Probs <- c(0.2, NA)

  expect_equal(ScidrPValueStars(Probs), c("ns", NA_character_))
  expect_equal(ScidrPValueStars(Probs, ns_label = "", na_label = ""), c("", ""))
  expect_equal(ScidrPValueStars(NA_real_, na_label = "n/a"), "n/a")
})

test_that("caption text states the same thresholds the code applies", {
  Caption <- ScidrStarCaptionText()

  expect_match(Caption, "p <= 0.05", fixed = TRUE)
  expect_match(Caption, "p <= 0.01", fixed = TRUE)
  expect_match(Caption, "p <= 0.001", fixed = TRUE)
  # ASCII only, for CRAN
  expect_false(any(grepl("[^\x01-\x7F]", Caption)))
})
