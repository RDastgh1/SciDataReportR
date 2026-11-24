ReportExampleData <- function(n = 200) {
  set.seed(24)
  df_Data <- data.frame(ID = seq_len(n), x = stats::rnorm(n),
    Visit = as.Date("2024-01-01") + c(seq_len(n - 1) - 1, 2000))
  df_Data$y <- df_Data$x + stats::rnorm(n, sd = 0.1)
  df_Data$y[1] <- 10
  df_Data
}

test_that("default percentile cutoffs flag the top scores without overriding other cutoffs", {
  skip_if_not_installed("isotree")
  skip_if_not_installed("robustbase")
  df_Data <- ReportExampleData()
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "ID",
    methods = c("isolation", "distance"), explain_isolation = FALSE)
  for (method in c("isolation", "distance")) {
    records <- result$Results[result$Results$Method == method, ]
    expect_equal(sum(records$Flag %in% TRUE), 2L)
    expect_true(all(records$Flag[is.finite(records$Score)] ==
      (records$Score[is.finite(records$Score)] >= records$Threshold[is.finite(records$Score)])))
  }
  expect_true(any(result$Results$.Row == 1 & result$Results$Flag %in% TRUE))
  expect_setequal(result$Metadata$FlagPercentileMethods, c("isolation", "distance"))
  expect_match(result$Metadata$ThresholdProvenance, "default top 1% percentile")

  none <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "isolation",
    explain_isolation = FALSE, flag_percentile = NULL)
  expect_true(all(is.na(none$Results$Flag)))
  expect_null(none$Metadata$FlagPercentile)

  explicit <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "isolation",
    explain_isolation = FALSE, thresholds = c(isolation = 0.5))
  expect_equal(unique(explicit$Results$Threshold), 0.5)
  expect_length(explicit$Metadata$FlagPercentileMethods, 0L)

  expect_error(DetectDataAnomalies(df_Data, "x", methods = "isolation", flag_percentile = 1),
    "strictly between 0 and 1")
})

test_that("native date and conditional flags are untouched by the percentile default", {
  df_Data <- ReportExampleData(60)
  result <- DetectDataAnomalies(df_Data, "Visit", id_var = "ID", methods = "dates")
  expect_false("dates" %in% result$Metadata$FlagPercentileMethods)
  expect_true(any(result$Results$Flag %in% TRUE))
})

test_that("percentile cutoffs apply within each stratum", {
  skip_if_not_installed("isotree")
  df_Data <- ReportExampleData()
  df_Data$Site <- rep(c("A", "B"), each = 100)
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "ID", stratify_by = "Site",
    methods = "isolation", explain_isolation = FALSE)
  thresholds <- unique(result$Results[c("Stratum", "Threshold")])
  expect_equal(nrow(thresholds), 2L)
  expect_false(thresholds$Threshold[1] == thresholds$Threshold[2])
  expect_equal(sum(result$Results$Flag %in% TRUE), 2L)
})

test_that("print and summary describe the screen without dumping the object", {
  df_Data <- ReportExampleData(60)
  result <- DetectDataAnomalies(df_Data, "Visit", id_var = "ID", methods = "dates")
  output <- capture.output(print(result))
  expect_true(any(grepl("Cross-sectional data anomaly screen", output)))
  expect_true(any(grepl("Records screened: 60", output)))
  expect_true(any(grepl("review candidates, not confirmed errors", output)))
  expect_true(any(grepl("CreateDataAnomalyReport", output)))
  expect_false(any(grepl("ReviewData", output)))

  df_Summary <- summary(result)
  expect_true(all(c("CutoffSource", "ScoreDefinition", "Threshold") %in% names(df_Summary)))
  expect_equal(nrow(df_Summary), nrow(result$Summary))

  queue <- CreateAnomalyReviewQueue(result)
  output <- capture.output(print(queue))
  expect_true(any(grepl("Anomaly review queue", output)))
  expect_true(any(grepl("Decisions imported: no", output)))
})

test_that("SummarizeAnomalyDecisions counts decisions and never modifies the queue", {
  df_Data <- ReportExampleData(60)
  df_Data$Visit[2] <- as.Date("1990-01-01")
  result <- DetectDataAnomalies(df_Data, "Visit", id_var = "ID", methods = "dates")
  queue <- CreateAnomalyReviewQueue(result)
  expect_gte(nrow(queue$Issues), 2L)

  empty <- SummarizeAnomalyDecisions(queue)
  expect_equal(empty$Progress$Reviewed, 0L)
  expect_equal(empty$Progress$Outstanding, nrow(queue$Issues))
  expect_equal(nrow(empty$Corrections), 0L)

  decided <- queue
  decided$Issues$ReviewStatus[1:2] <- c("Correction needed", "Confirmed correct")
  decided$Issues$ProposedCorrection[1] <- "2024-03-01"
  decided$Issues$Reviewer[1:2] <- "AB"
  df_Decided <- SummarizeAnomalyDecisions(decided)
  expect_equal(df_Decided$Progress$Reviewed, 2L)
  expect_equal(df_Decided$Progress$WithProposedCorrection, 1L)
  expect_equal(nrow(df_Decided$Corrections), 1L)
  expect_equal(df_Decided$Corrections$ProposedCorrection, "2024-03-01")
  expect_equal(sum(df_Decided$Counts$Issues), nrow(queue$Issues))
  expect_setequal(df_Decided$Counts$ReviewStatus,
    unique(ifelse(decided$Issues$ReviewStatus == "", "Unreviewed", decided$Issues$ReviewStatus)))
  expect_equal(names(SummarizeAnomalyDecisions(decided, by = "variable"))[1], "Counts")
  expect_true("Reviewer" %in% names(SummarizeAnomalyDecisions(decided, by = "reviewer")$Counts))
  expect_true(all(queue$Issues$ReviewStatus == "Unreviewed"))
  expect_true(all(is.na(queue$Issues$ProposedCorrection)))
  expect_error(SummarizeAnomalyDecisions(result), "must be an AnomalyReviewQueue")
})

test_that("CreateDataAnomalyReport writes shareable files and returns their paths", {
  skip_if_not_installed("DT")
  skip_if_not_installed("htmltools")
  directory <- withr::local_tempdir()
  df_Data <- ReportExampleData(60)
  review <- CreateDataAnomalyReport(df_Data, c("x", "y", "Visit"), id_var = "ID",
    methods = "dates", workbook = FALSE, output_dir = directory, Title = "Review test")

  expect_s3_class(review$Queue, "AnomalyReviewQueue")
  expect_named(review$Results, "Cross-sectional")
  expect_true(file.exists(review$Files[["Explorer"]]))
  expect_match(review$Files[["Explorer"]], "DataAnomalyReview_")
  page <- paste(readLines(review$Files[["Explorer"]], warn = FALSE), collapse = "\n")
  expect_true(grepl("Review test", page, fixed = TRUE))
  expect_true(grepl("How to read this page", page, fixed = TRUE))

  skip_if_not(rmarkdown::pandoc_available("1.17") || isTRUE(ScidrAnomalyFindPandoc()),
    "pandoc unavailable for self-contained HTML")
  expect_false(grepl("_files/", page, fixed = TRUE))
  expect_false(dir.exists(file.path(directory, "DataAnomalyReview_files")))

  expect_error(CreateDataAnomalyReport(df_Data, "Visit", methods = "dates",
    workbook = FALSE, output_dir = directory), "overwrite = TRUE")
  expect_error(CreateDataAnomalyReport(df_Data, "Visit", explorer = FALSE, workbook = FALSE),
    "at least one")
  expect_error(CreateDataAnomalyReport(df_Data, "Visit", output_dir = file.path(directory, "absent")),
    "Directory does not exist")
  expect_error(CreateDataAnomalyReport(df_Data, "Visit", longitudinal = list(1),
    output_dir = directory), "uniquely named list")
})

test_that("CreateDataAnomalyReport adds a longitudinal analysis and the Excel workbook", {
  skip_if_not_installed("DT")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("jsonlite")
  directory <- withr::local_tempdir()
  set.seed(24)
  df_Long <- data.frame(ID = rep(seq_len(30), each = 3), Visit = rep(1:3, times = 30),
    score = stats::rnorm(90))
  df_Long$score[5] <- 40
  review <- CreateDataAnomalyReport(df_Long, "score", id_var = "ID",
    methods = "distance", output_dir = directory,
    longitudinal = list(variables = "score", time_var = "Visit", methods = "change"))

  expect_named(review$Results, c("Cross-sectional", "Longitudinal"))
  expect_setequal(names(review$Queue$Sources), c("Cross-sectional", "Longitudinal"))
  expect_true(all(file.exists(review$Files)))
  expect_setequal(names(review$Files), c("Explorer", "Workbook"))
  expect_match(review$Files[["Workbook"]], "\\.xlsx$")
  sheets <- openxlsx::getSheetNames(review$Files[["Workbook"]])
  expect_true("Review" %in% sheets)
})

test_that("summary names the cutoff source for native, user, default and unflagged methods", {
  skip_if_not_installed("isotree")
  df_Data <- ReportExampleData(150)
  native <- summary(DetectDataAnomalies(df_Data, c("Visit", "x"), id_var = "ID",
    methods = c("dates", "isolation"), explain_isolation = FALSE))
  expect_equal(native$CutoffSource[native$Method == "dates"], "Native cutoff")
  expect_equal(native$CutoffSource[native$Method == "isolation"], "Default percentile cutoff")

  user <- summary(DetectDataAnomalies(df_Data, c("x", "y"), methods = "isolation",
    explain_isolation = FALSE, thresholds = c(isolation = 0.5)))
  expect_equal(user$CutoffSource, "User cutoff")

  none <- summary(DetectDataAnomalies(df_Data, c("x", "y"), methods = "isolation",
    explain_isolation = FALSE, flag_percentile = NULL))
  expect_equal(none$CutoffSource, "No cutoff; scores only")
})

test_that("the explorer says when a method stopped before finishing", {
  skip_if_not_installed("DT")
  set.seed(5)
  df_Example <- data.frame(ID = 1:80, BMI = round(stats::rnorm(80, 27, 4), 2),
    Age = round(stats::rnorm(80, 50, 10)), Lab = stats::rnorm(80))
  df_Example$BMI[41] <- 276249.07
  result <- DetectDataAnomalies(df_Example, c("BMI", "Age", "Lab"), id_var = "ID",
    methods = "conditional")
  html <- paste(as.character(ExploreDataAnomalies(result, IncludePlots = FALSE)), collapse = "")
  expect_false(grepl("Incomplete screening", html))
  result$Diagnostics <- rbind(as.data.frame(result$Diagnostics),
    data.frame(Method = "conditional", Stratum = 1L, Issue = "backend_failure",
      Detail = "required resulting string length 874597747 is greater than maximal 8192",
      Rows = I(list(integer()))))
  html <- paste(as.character(ExploreDataAnomalies(result, IncludePlots = FALSE)), collapse = "")
  expect_true(grepl("Incomplete screening", html))
  expect_true(grepl("874597747", html))
})

test_that("both review tables carry dropdown filters for their categorical columns", {
  skip_if_not_installed("DT")
  set.seed(9)
  n <- 90
  df_Example <- data.frame(ID = sprintf("P%02d", 1:n), BMI = round(stats::rnorm(n, 27, 4), 2),
    Age = round(stats::rnorm(n, 50, 10)), Lab = stats::rnorm(n))
  df_Example$BMI[41] <- 276249.07
  df_Example$Lab[c(3, 9, 15)] <- NA
  df_Codebook <- data.frame(Variable = c("BMI", "Age", "Lab"),
    Label = c("Body mass index", "Age", "Lab value"),
    Category = c("Anthropometry", "Demographics", "Chemistry"))
  result <- DetectDataAnomalies(df_Example, c("BMI", "Age", "Lab"), id_var = "ID", codebook = df_Codebook)
  html <- paste(as.character(ExploreDataAnomalies(result, IncludePlots = FALSE, codebook = df_Codebook)), collapse = "")
  expect_true(grepl("addFilters", html, fixed = TRUE))
  for (title in c("Extreme value", "Anomaly type", "Severity tier", "Codebook category", "Method")) {
    expect_true(grepl(paste0("title: '", title, "'"), html, fixed = TRUE), info = title)
  }
  # The variables table drives filtering from its own columns, not the queue's.
  expect_true(grepl("name: 'CodebookCategory'", html, fixed = TRUE))
  expect_true(grepl("name: 'CodebookCategoryKeys'", html, fixed = TRUE))
})
