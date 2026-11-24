test_that("automatic selection supports unrevalued frames and protects source classes", {
  df_Raw <- data.frame(ID = sprintf("%03d", 1:40), Site = rep(c("A", "B"),20),
    Value = seq_len(40), Category = rep(c("yes", "no"),20), Logical = rep(c(TRUE,FALSE),20),
    Date = as.Date("2020-01-01") + seq_len(40))
  df_Raw$Extra <- I(as.list(seq_len(40)))
  before <- df_Raw
  result <- DetectDataAnomalies(df_Raw, id_var = "ID", stratify_by = "Site", methods = "dates")
  expect_setequal(result$Metadata$Variables, c("Value","Category","Logical","Date"))
  expect_true("unsupported_auto_excluded" %in% result$Diagnostics$Issue)
  expect_true("missing_label" %in% result$Diagnostics$Issue)
  expect_identical(df_Raw, before)
  expect_identical(result$ReviewData$Date, df_Raw$Date)
  expect_error(DetectDataAnomalies(df_Raw, "Extra", methods="dates"), "Unsupported")
})

test_that("native mixed backends accept plain unrevalued data with default selection", {
  skip_if_not_installed("isotree")
  set.seed(22)
  df_Raw <- data.frame(Value = rnorm(60), Category = rep(c("yes", "no"),30), Logical = rep(c(TRUE,FALSE),30))
  before <- df_Raw
  result <- DetectDataAnomalies(df_Raw, methods="isolation", method_args=list(isolation=list(ntrees=20,nthreads=1)))
  expect_true(any(result$Results$Status == "evaluated"))
  expect_identical(df_Raw, before)
})

test_that("print projection includes the full queue and score-only checks are unclassified", {
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  before <- queue
  html <- htmltools::renderTags(ExploreDataAnomalies(queue, IncludePlots=FALSE))$html
  expect_match(html, "Complete review queue")
  expect_match(html, "scidr-static-table")
  expect_match(html, "Source row 2")
  expect_match(html, "Technical scrolling tables are omitted")
  expect_identical(queue, before)
  result <- ReviewQueueFixture()
  result$Results$Flag <- NA
  html <- htmltools::renderTags(ExploreDataAnomalies(result, IncludePlots=FALSE))$html
  expect_match(html, "No cutoff supplied; scores only")
  expect_match(html, "Not classified")
})
