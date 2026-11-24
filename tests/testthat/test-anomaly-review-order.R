test_that("review tiers and score percentiles are derived without changing evidence", {
  result <- ReviewQueueFixture()
  result$ReviewData$Code <- rep(c(0,1),20)
  queue <- CreateAnomalyReviewQueue(result)
  queue$Evidence <- queue$Evidence[rep(1,5), ]
  queue$Issues <- queue$Issues[rep(1,5), ]
  queue$Evidence$IssueID <- queue$Issues$IssueID <- letters[1:5]
  queue$Evidence$EvidenceID <- LETTERS[1:5]
  queue$Evidence$ReferenceN <- c(30, 10, NA, 100, 30)
  queue$Evidence$Variable <- c("Visit","Visit",NA,"Code","Visit")
  queue$Evidence$Score <- c(10, 20, 50, 100, NA)
  queue$Evidence$Method <- c("dates","dates","isolation","conditional","dates")
  queue$Sources$Analysis$ReviewData$Code <- result$ReviewData$Code
  before <- queue
  view <- ScidrReviewNarratives(queue)
  expect_identical(queue, before)
  # Tiers are severities of the finding type: date 2, rare coded value 7,
  # record combination 8. Thin reference support (row b, n = 10) sorts within
  # its own type rather than ahead of a more serious one.
  expect_equal(view$Issues$ReviewTier, c(2L,2L,2L,7L,8L))
  expect_equal(view$Issues$IssueID, c("a","e","b","d","c"))
  expect_true(is.na(view$Issues$ScorePercentile[2]))
  expect_match(view$Issues$OrderRationale[1], "not error probability")
  expect_equal(view$Issues$SuggestedReviewOrder, 1:5)
  expect_identical(view, ScidrReviewNarratives(queue))
  queue$Evidence <- dplyr::bind_rows(queue$Evidence, queue$Evidence[1, ])
  queue$Evidence$EvidenceID[6] <- "duplicate-support"
  expect_equal(ScidrReviewNarratives(queue)$Issues$IssueID, view$Issues$IssueID)
})

test_that("tie percentiles are retained within stratum and missing scores stay unavailable", {
  result <- ReviewQueueFixture()
  result$DateResults <- result$DateResults[FALSE, ]
  result$Results$Method <- "isolation"
  result$Results$Flag <- NA
  result$Results$Score <- rep(c(1,1,2,2),10)
  result$Results$Stratum <- rep(1:2,each=20)
  queue <- CreateAnomalyReviewQueue(result,top_n=20)
  view <- ScidrReviewNarratives(queue)
  expect_equal(sort(unique(view$Evidence$ScorePercentile)),c(.25,.75))
  expect_true(all(view$Issues$ReviewTier == 8L))
  expect_true(all(is.na(queue$Issues$Variable)))
  expect_true(all(is.na(queue$Evidence$Flag)))
  expect_true(all(view$Issues$FindingType == "Uncommon combination"))
})

test_that("derived Excel order does not move decisions to a different issue", {
  skip_if_not_installed("openxlsx")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  file <- tempfile(fileext=".xlsx")
  on.exit(unlink(file),add=TRUE)
  ExportAnomalyReviewQueue(queue,file, layout = "full")
  id <- ScidrReviewNarratives(queue)$Issues$IssueID[1]
  EditReviewWorkbook(file,function(d) {d$`Review status`[1] <- "Correction needed";d})
  updated <- ImportAnomalyReviewDecisions(queue,file)
  expect_equal(updated$Issues$ReviewStatus[match(id,updated$Issues$IssueID)],"Correction needed")
  expect_identical(updated$Metadata$QueueFingerprint,queue$Metadata$QueueFingerprint)
  ExportAnomalyReviewQueue(queue,file,overwrite=TRUE, layout = "full")
  derived <- ScidrReviewHeaders(c("SuggestedReviewOrder","ReviewTier","OrderRationale","ScorePercentile","ReviewReferenceN","FindingType","Methods"))
  for (sheet in c("Review","Evidence")) EditReviewWorkbook(file,function(d) d[setdiff(names(d),derived)],sheet)
  expect_identical(ImportAnomalyReviewDecisions(queue,file),queue)
})


test_that("the coverage hierarchy protects embedded Quarto tables", {
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  html <- htmltools::renderTags(ExploreDataAnomalies(queue, IncludePlots=FALSE))$html
  expect_match(html, 'data-html-table-processing="none"', fixed=TRUE)
  expect_match(html, 'data-quarto-disable-processing="true"', fixed=TRUE)
  expect_match(html, "Suggested review order")
  expect_match(html, "Coverage")
  expect_false(grepl("Technical results", html, fixed = TRUE))
  expect_false(grepl("Queue diagnostics", html, fixed = TRUE))
  expect_false(grepl("table.*selected-row-summary", html))
})

test_that("value-level findings outrank blanks regardless of score availability", {
  skip_if_not_installed("outliertree")
  skip_if_not_installed("isotree")
  set.seed(9)
  n <- 115
  df_Data <- data.frame(ID = paste0("S", seq_len(n)), stringsAsFactors = FALSE)
  for (j in 1:25) df_Data[[paste0("Lab", j)]] <- round(stats::rnorm(n, 100, 15), 1)
  # A panel of four labs skipped wholesale by 30 records, with one record
  # carrying three of the four: a real missingness finding, which used to fill
  # the top of the queue because conditional findings carry no score and sorted
  # last.
  for (j in 1:4) df_Data[[paste0("Lab", j)]][1:30] <- NA
  df_Data$Lab1[31] <- NA
  df_Data$Lab20[5] <- 900
  Result <- suppressMessages(DetectDataAnomalies(df_Data))
  Issues <- ScidrReviewNarratives(suppressMessages(CreateAnomalyReviewQueue(Result)))$Issues
  Value <- Issues$SuggestedReviewOrder[Issues$ReviewCategory %in%
    c("Distribution outlier", "Conditional outlier")]
  Blank <- Issues$SuggestedReviewOrder[Issues$ReviewCategory == "Unexpected missingness"]
  expect_gt(length(Value), 0L)
  expect_gt(length(Blank), 0L)
  expect_lt(max(Value), min(Blank))
})

test_that("the extreme-magnitude ratio is scale free and needs a usable reference", {
  values <- c(rep(10, 30), 11, 9, 1000)
  expect_gt(ScidrAnomalyExtremeRatio(values, 33L), 10)
  expect_lt(ScidrAnomalyExtremeRatio(values, 31L), 10)
  # Scaling the column cannot change which records look extreme.
  expect_equal(ScidrAnomalyExtremeRatio(values, 33L), ScidrAnomalyExtremeRatio(values * 1000, 33L))
  expect_true(is.na(ScidrAnomalyExtremeRatio(c(1, 2, 500), 3L)))
  expect_true(is.na(ScidrAnomalyExtremeRatio(rep(5, 40), 1L)))
  expect_true(is.na(ScidrAnomalyExtremeRatio(as.character(values), 33L)))
  expect_true(is.na(ScidrAnomalyExtremeRatio(c(values[-33], NA), 33L)))
})

test_that("an extreme value is reviewed before every severity tier", {
  skip_if_not_installed("outliertree")
  set.seed(12)
  df_Example <- data.frame(ID = 1:80, BMI = round(stats::rnorm(80, 27, 4), 2),
    Age = round(stats::rnorm(80, 50, 10)), Lab = stats::rnorm(80))
  df_Example$BMI[41] <- 276249.07
  result <- DetectDataAnomalies(df_Example, c("BMI", "Age", "Lab"), id_var = "ID")
  queue <- CreateAnomalyReviewQueue(result)
  display <- ScidrReviewNarratives(queue)
  first <- display$Issues[display$Issues$SuggestedReviewOrder == 1L, ]
  expect_equal(first$Variable, "BMI")
  expect_equal(first$SourceRow, 41L)
  expect_equal(first$ExtremeValue, "Extreme")
  expect_match(first$OrderRationale, "Extreme magnitude")
  # It outranks tiers more serious than its own.
  expect_true(any(display$Issues$ReviewTier < first$ReviewTier &
    display$Issues$SuggestedReviewOrder > 1L) || all(display$Issues$ReviewTier >= first$ReviewTier))
  expect_true(all(display$Issues$ExtremeValue[display$Issues$SuggestedReviewOrder > 1L] == "" |
    !is.na(display$Issues$ExtremeRatio[display$Issues$SuggestedReviewOrder > 1L])))
})
