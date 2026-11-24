test_that("default screening matches direct native replacement predictions", {
  skip_if_not_installed("isotree")
  set.seed(7)
  df_Data <- tibble::tibble(`value x` = rnorm(35), y = rnorm(35),
    Category = factor(rep(c("A", "B", NA), length.out = 35)),
    Date = as.POSIXct("2020-01-01", tz = "America/New_York") + seq_len(35) * 86400,
    Constant = 1)
  df_Data$y[2] <- NA
  before <- df_Data
  rng <- .Random.seed
  Detect <- function(...) DetectDataAnomalies(df_Data, methods = "isolation", seed = 27,
    method_args = list(isolation = list(ntrees = 20, nthreads = 1)), ...)
  result <- Detect(explanation_args = list(rows = c(1, 2), batch_size = 1))
  expect_identical(.Random.seed, rng)
  expect_identical(df_Data, before)
  expect_identical(result$ReviewData$Date, df_Data$Date)
  expect_equal(nrow(result$FeatureContributions), 0)
  expect_equal(result$ExplanationMetadata$Method, "sensitivity")
  expect_true(all(is.na(result$ExplanationSummary$Baseline)))
  expect_true(all(is.na(result$ExplanationSummary$ReconstructionResidual)))
  expect_setequal(result$ExplanationSummary$.Row[result$ExplanationSummary$Status == "Explained"], c(1, 2))
  expect_false("Constant" %in% result$FeatureSensitivity$Variable)
  metadata <- result$ExplanationMetadata$Strata[[1]]
  df_Model <- as.data.frame(df_Data[metadata$ModelingVariables])
  df_Model$Date <- ScidrAnomalyDateDays(df_Model$Date) - result$Metadata$DateConversions$Date$ReferenceDays
  fit <- isotree::isolation.forest(df_Model, seed = 27, ntrees = 20, nthreads = 1)
  for (i in c(1, 2)) for (variable in metadata$ModelingVariables) {
    df_Replacement <- df_Model[rep(i, length(metadata$ReferenceRows)), , drop = FALSE]
    df_Replacement[[variable]] <- df_Model[[variable]][metadata$ReferenceRows]
    score <- as.numeric(predict(fit, df_Replacement, type = "score", nthreads = 1))
    evidence <- result$FeatureSensitivity[result$FeatureSensitivity$.Row == i & result$FeatureSensitivity$Variable == variable, ]
    native <- result$ExplanationSummary$Score[result$ExplanationSummary$.Row == i]
    expect_equal(evidence$ScoreDecrease, native - mean(score))
    expect_equal(evidence$ReplacementScoreMin, min(score))
    expect_equal(evidence$ReplacementScoreMax, max(score))
    expect_equal(evidence$ReplacementN, 3)
  }
  batched <- Detect(explanation_args = list(rows = c(1, 2), batch_size = 64))
  expect_equal(result$FeatureSensitivity, batched$FeatureSensitivity)
  repeated <- Detect(explanation_args = list(rows = c(1, 2), batch_size = 1))
  expect_identical(result$ExplanationMetadata$Fingerprint, repeated$ExplanationMetadata$Fingerprint)
  disabled <- Detect(explain_isolation = FALSE)
  expect_equal(result$Results, disabled$Results)
  expect_identical(ScidrReviewSnapshot(result), ScidrReviewSnapshot(disabled))
  expect_identical(CreateAnomalyReviewQueue(result)$Issues$IssueID, CreateAnomalyReviewQueue(disabled)$Issues$IssueID)
})

test_that("screening supports single features, strata and unavailable records", {
  skip_if_not_installed("isotree")
  df_Data <- data.frame(ID = rep(c("001", NA, "003"), 10), Group = rep(c("A", "B", NA), each = 10), x = seq_len(30))
  df_Data$x[1] <- NA
  result <- DetectDataAnomalies(df_Data, id_var = "ID", stratify_by = "Group", methods = "isolation",
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(rows = c(1, 2, 12, 22), replacements = 100))
  expect_setequal(result$FeatureSensitivity$.Row, c(2, 12))
  expect_true(all(result$FeatureSensitivity$ReplacementN <= 10))
  expect_true(all(result$FeatureSensitivity$Variable == "x"))
  expect_true(all(result$ExplanationSummary$Status[c(1, 22)] == "Unavailable"))
  for (metadata in result$ExplanationMetadata$Strata) {
    expect_length(unique(df_Data$Group[metadata$ReferenceRows]), 1)
    expect_true(all(metadata$ReferenceRows %in% metadata$BackgroundRows))
  }
  expect_error(DetectDataAnomalies(df_Data, explanation_args = list(method = "other")), "method")
  expect_error(DetectDataAnomalies(df_Data, explanation_args = list(batch_size = 0)), "batch_size")
  expect_error(DetectDataAnomalies(df_Data, explanation_args = list(replacements = 0)), "replacements")
})

test_that("screening failures preserve scores and do not require fastshap", {
  skip_if_not_installed("isotree")
  original <- stats::predict
  Predict <- function(object, newdata, ...) {
    if (nrow(newdata) == 3) stop("Synthetic replacement failure")
    original(object, newdata, ...)
  }
  testthat::local_mocked_bindings(predict = Predict, .package = "stats")
  df_Data <- data.frame(x = seq_len(20))
  result <- DetectDataAnomalies(df_Data, methods = "isolation", method_args = list(isolation = list(ntrees = 20)), explanation_args = list(top_n = 1))
  expect_true(all(result$Results$Status == "evaluated"))
  expect_equal(sum(result$ExplanationSummary$Status == "Unavailable"), 1)
  expect_equal(nrow(result$FeatureSensitivity), 0)
  expect_true("explanation_failure" %in% result$Diagnostics$Issue)
})

test_that("screening Excel evidence, method labels and legacy workbooks remain safe", {
  skip_if_not_installed("isotree")
  skip_if_not_installed("openxlsx")
  df_Data <- data.frame(ID = sprintf("%03d", 1:30), x = seq_len(30), y = sin(seq_len(30)))
  result <- DetectDataAnomalies(df_Data, id_var = "ID", methods = "isolation", method_args = list(isolation = list(ntrees = 20)), explanation_args = list(top_n = 1))
  queue <- CreateAnomalyReviewQueue(result, top_n = 1)
  meta <- ScidrReviewExplanationMeta(queue)[[1]]
  expect_equal(meta$method, "sensitivity")
  expect_true(length(meta$values) > 0)
  expect_true(is.finite(meta$values[[1]]$ScoreDecrease))
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  exported <- openxlsx::read.xlsx(file, "Contributions", check.names = FALSE, sep.names = " ")
  expect_true(all(exported[["Explanation method"]] == "sensitivity"))
  EditReviewWorkbook(file, function(d) { d[["Score decrease after replacement"]][1] <- "999"; d }, "Contributions")
  expect_error(ImportAnomalyReviewDecisions(queue, file), "Contributions.*changed")
  changed <- queue
  changed$Sources$Analysis$FeatureSensitivity$ScoreDecrease[1] <- 999
  expect_error(ScidrReviewContributions(changed), "fingerprint changed")
  ExportAnomalyReviewQueue(queue, file, overwrite = TRUE, layout = "full")
  wb <- openxlsx::loadWorkbook(file)
  openxlsx::removeWorksheet(wb, "Contributions")
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  html <- as.character(htmltools::renderTags(ExploreDataAnomalies(queue))$html)
  expect_match(html, "Columns worth checking", fixed = TRUE)
  expect_match(html, "Effects are not additive", fixed = TRUE)
})

test_that("signed screening effects are not additive attributions", {
  df_Data <- data.frame(x = c(10, 1, 2), y = c(10, 1, 2))
  setup <- list(seed = 2024, data = df_Data, labels = c(x = "X label", y = "Y label"),
    date_conversions = list(), Add = function(...) stop("Unexpected diagnostic"))
  testthat::local_mocked_bindings(predict = function(object, newdata, ...) newdata$x - newdata$y, .package = "stats")
  args <- ScidrAnomalyExplanationArgs(list(rows = 1), TRUE, "isolation", 3)
  result <- ScidrAnomalyIsolationSensitivity(NULL, df_Data, 1:3, rep(0, 3), 1, setup, args)
  expect_equal(result$Sensitivity$ScoreDecrease, c(17 / 3, -17 / 3))
  expect_equal(result$Sensitivity$PositiveRank[1], 1)
  expect_true(is.na(result$Sensitivity$PositiveRank[2]))
  expect_true(all(is.na(result$Summary$SumContributions)))
})
