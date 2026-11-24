test_that("variables narrows the queue without refitting", {
  set.seed(5)
  df_Example <- data.frame(Age = c(rnorm(39, 50, 5), 500), Weight = c(rnorm(39, 70, 6), 400),
    Height = rnorm(40, 170, 8))
  result <- DetectDataAnomalies(df_Example, c("Age", "Weight", "Height"), methods = "conditional")
  full <- CreateAnomalyReviewQueue(result)
  subset <- CreateAnomalyReviewQueue(result, variables = "Age")

  kept <- subset$Evidence$Variable
  expect_true(all(is.na(kept) | kept == "Age"))
  expect_true(all(subset$Issues$IssueID %in% full$Issues$IssueID))
  # Issue identities are positions, so narrowing must not renumber a finding.
  shared <- intersect(subset$Issues$IssueID, full$Issues$IssueID)
  expect_equal(subset$Issues$Observed[match(shared, subset$Issues$IssueID)],
    full$Issues$Observed[match(shared, full$Issues$IssueID)])
  # Coverage still describes the whole screen.
  expect_equal(subset$Summary, full$Summary)
  expect_equal(subset$Metadata$Variables, "Age")
  expect_equal(subset$Issues$EvidenceCount,
    as.integer(table(factor(subset$Evidence$IssueID, levels = subset$Issues$IssueID))[subset$Issues$IssueID]))
  expect_silent(ScidrReviewValidate(subset))
})

test_that("variables rejects columns that were not screened", {
  df_Example <- data.frame(Visit = as.Date("2020-01-01") + seq_len(40))
  df_Example$Visit[1] <- as.Date("2099-01-01")
  result <- DetectDataAnomalies(df_Example, "Visit", methods = "dates")
  expect_error(CreateAnomalyReviewQueue(result, variables = "Nonexistent"), "not screened")
  expect_error(CreateAnomalyReviewQueue(result, variables = c("Visit", "Visit")), "unique column names")
})

test_that("the explorer honours variables and says that it did", {
  skip_if_not_installed("DT")
  set.seed(6)
  df_Example <- data.frame(Age = c(rnorm(39, 50, 5), 500), Weight = c(rnorm(39, 70, 6), 400))
  result <- DetectDataAnomalies(df_Example, c("Age", "Weight"), methods = "conditional")
  html <- paste(as.character(ExploreDataAnomalies(result, variables = "Age", IncludePlots = FALSE)), collapse = "")
  expect_true(grepl("selected variable\\(s\\) only", html))
})
