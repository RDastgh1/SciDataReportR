make_messy_validation <- function() {
  baseline <- data.frame(
    id = 1:5,
    age = c(50, 61, 45, 58, 63),
    sex = c("F", "M", "F", "F", "M")
  )
  labs <- data.frame(
    id = c(1, 2, 2, 6),
    sex = c("F", "M", "M", "F"),
    glucose = c(90, 110, 115, 100)
  )
  merged <- dplyr::left_join(baseline, labs, by = "id")

  ValidateMerge(baseline, labs, merged, Keys = "id")
}

test_that("ExploreMergeValidation works without a Detail argument", {
  skip_if_not_installed("reactable")
  skip_if_not_installed("htmltools")

  v <- make_messy_validation()

  dashboard <- ExploreMergeValidation(v)
  expect_s3_class(dashboard, "shiny.tag.list")

  html <- as.character(htmltools::renderTags(dashboard)$html)

  # Fingerprint/summary cards are gone; context folded into checks table
  expect_no_match(html, "sdr-fingerprint")
  expect_no_match(html, "sdr-card-grid")

  # Compact default: accordion sections present but collapsed
  expect_match(html, "Coverage explorer \\(4 unmatched\\)")
  expect_match(html, "Duplicate-variable conflicts \\(1 variable\\)")
  expect_no_match(html, "<details[^>]*open")
})

test_that("ExploreMergeValidation Detail = 'Full' expands accordions", {
  skip_if_not_installed("reactable")
  skip_if_not_installed("htmltools")

  v <- make_messy_validation()

  dashboard <- ExploreMergeValidation(v, Detail = "Full")
  html <- as.character(htmltools::renderTags(dashboard)$html)

  expect_match(html, "<details[^>]*open")
})

test_that("ExploreMergeValidation omits empty sections", {
  skip_if_not_installed("reactable")
  skip_if_not_installed("htmltools")

  baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
  labs <- data.frame(id = 1:4, glucose = c(90, 110, 100, 95))
  merged <- dplyr::left_join(baseline, labs, by = "id")

  v <- ValidateMerge(baseline, labs, merged, Keys = "id")

  dashboard <- ExploreMergeValidation(v)
  html <- as.character(htmltools::renderTags(dashboard)$html)

  expect_match(html, "Validation checks")
  expect_no_match(html, "Coverage explorer")
  expect_no_match(html, "Duplicate-variable conflicts")
  expect_no_match(html, "Suggested actions")
})
