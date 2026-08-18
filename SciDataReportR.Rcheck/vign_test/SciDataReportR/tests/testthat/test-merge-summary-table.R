test_that("merge_summary_table binds safe_merge logs", {
  baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
  labs_clean <- data.frame(id = 1:4, glucose = c(90, 110, 100, 95))
  labs_duplicated <- data.frame(id = c(1, 2, 2), platelet = c(200, 180, 185))

  m1 <- safe_merge(baseline, labs_clean, by = "id", name = "Clean merge")
  m2 <- safe_merge(baseline, labs_duplicated, by = "id", name = "Duplicated merge")

  combined <- merge_summary_table(list(m1$log, m2$log))

  expect_s3_class(combined, "tbl_df")
  expect_equal(nrow(combined), 2)
  expect_equal(combined$Merge, c("Clean merge", "Duplicated merge"))
  expect_equal(combined$Status, c("PASS", "FAIL"))
})

test_that("merge_summary_table filters to flagged merges", {
  baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
  labs_clean <- data.frame(id = 1:4, glucose = c(90, 110, 100, 95))
  labs_duplicated <- data.frame(id = c(1, 2, 2), platelet = c(200, 180, 185))

  m1 <- safe_merge(baseline, labs_clean, by = "id", name = "Clean merge")
  m2 <- safe_merge(baseline, labs_duplicated, by = "id", name = "Duplicated merge")

  flagged <- merge_summary_table(list(m1$log, m2$log), flagged_only = TRUE)

  expect_equal(nrow(flagged), 1)
  expect_equal(flagged$Merge, "Duplicated merge")

  # Full safe_merge results are also accepted
  flagged_from_results <- merge_summary_table(list(m1, m2), flagged_only = TRUE)
  expect_identical(flagged, flagged_from_results)
})

test_that("merge_summary_table rejects invalid input", {
  expect_error(merge_summary_table(list()), "non-empty list")
  expect_error(merge_summary_table(list(1)), "log tibble")
})
