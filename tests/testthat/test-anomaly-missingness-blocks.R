test_that("variables that go missing together are recovered as groups", {
  set.seed(4)
  n <- 100
  missing <- cbind(
    matrix(FALSE, n, 4, dimnames = list(NULL, c("wbc", "rbc", "hgb", "hct"))),
    matrix(FALSE, n, 3, dimnames = list(NULL, c("sodium", "potassium", "chloride"))),
    matrix(FALSE, n, 2, dimnames = list(NULL, c("age", "sex"))))
  missing[1:30, c("wbc", "rbc", "hgb", "hct")] <- TRUE
  missing[31:32, "wbc"] <- TRUE
  missing[60:79, c("sodium", "potassium", "chloride")] <- TRUE
  missing[5, "age"] <- TRUE
  args <- list(min_block_events = 5L, min_block_size = 3L, block_similarity = 0.8)
  blocks <- ScidrAnomalyMissingnessBlocks(missing, args)
  expect_equal(length(blocks), 2L)
  expect_setequal(blocks[[1]], c("hct", "hgb", "rbc", "wbc"))
  expect_setequal(blocks[[2]], c("chloride", "potassium", "sodium"))
  # Two unrelated forms must not chain into one through shared skippers.
  expect_false(any(vapply(blocks, function(block) "wbc" %in% block && "sodium" %in% block, logical(1))))
  # A lone blank belongs to nothing.
  expect_false("age" %in% unlist(blocks))
  expect_equal(ScidrAnomalyMissingnessBlocks(missing, utils::modifyList(args, list(min_block_size = 5L))), list())
})

test_that("a gap inside a form is a finding and a skipped form is not", {
  set.seed(7)
  n <- 120
  df_Data <- data.frame(ID = sprintf("P%03d", seq_len(n)), Age = round(stats::rnorm(n, 50, 10)),
    wbc = stats::rnorm(n, 7), rbc = stats::rnorm(n, 5), hgb = stats::rnorm(n, 14), hct = stats::rnorm(n, 42))
  df_Data[1:35, c("wbc", "rbc", "hgb", "hct")] <- NA
  df_Data$hgb[100] <- NA
  df_Data$Age[50] <- NA
  result <- suppressMessages(DetectDataAnomalies(df_Data,
    c("Age", "wbc", "rbc", "hgb", "hct"), id_var = "ID", methods = "missingness"))
  flagged <- subset(result$CellResults, Flag %in% TRUE)
  expect_equal(nrow(flagged), 1L)
  expect_equal(flagged$.Row, 100L)
  expect_equal(flagged$Variable, "hgb")
  expect_match(flagged$Explanation, "3 of the 4 variables in a group that go missing together")
  expect_match(flagged$Explanation, "hct, hgb, rbc, wbc", fixed = TRUE)
  statistics <- flagged$GroupStatistics[[1]]
  expect_equal(statistics$block_size, 4L)
  expect_equal(statistics$observed_in_block, 3L)
  # Whole-form skippers are not findings, and neither is the unplaced lone blank.
  expect_false(any(subset(result$CellResults, .Row %in% 1:35)$Flag %in% TRUE))
  expect_false("Age" %in% result$CellResults$Variable)
  blocks <- result$Metadata$Settings$missingness$Blocks
  expect_equal(nrow(blocks), 1L)
  expect_equal(blocks$Members, "hct, hgb, rbc, wbc")
  expect_equal(blocks$Started, 85L)
})

test_that("a value rule still wins over the group it would otherwise share", {
  set.seed(8)
  n <- 120
  df_Data <- data.frame(ID = sprintf("P%03d", seq_len(n)),
    Sex = rep(c("F", "M"), length.out = n), Menopause = stats::rnorm(n, 50),
    wbc = stats::rnorm(n, 7), rbc = stats::rnorm(n, 5), hgb = stats::rnorm(n, 14))
  df_Data$Menopause[df_Data$Sex == "M"] <- NA
  df_Data$Menopause[which(df_Data$Sex == "F")[1:3]] <- NA
  result <- suppressMessages(DetectDataAnomalies(df_Data,
    c("Sex", "Menopause", "wbc", "rbc", "hgb"), id_var = "ID", methods = "missingness"))
  flagged <- subset(result$CellResults, Flag %in% TRUE)
  expect_true(all(flagged$Variable == "Menopause"))
  expect_equal(nrow(flagged), 3L)
  expect_match(flagged$Explanation[1], "Sex")
  # No man is ever a finding, however complete the column looks overall.
  expect_false(any(df_Data$Sex[flagged$.Row] == "M"))
})

test_that("the page shows the groups and the completeness of every blank column", {
  skip_if_not_installed("DT")
  set.seed(9)
  n <- 100
  df_Data <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n), c = stats::rnorm(n), d = stats::rnorm(n))
  df_Data[1:30, c("a", "b", "c")] <- NA
  df_Data$a[31] <- NA
  df_Data$d[40] <- NA
  result <- suppressMessages(DetectDataAnomalies(df_Data, c("a", "b", "c", "d"), methods = "missingness"))
  html <- paste(as.character(ExploreDataAnomalies(result, IncludePlots = FALSE)), collapse = "")
  expect_true(grepl("Groups of variables that go missing together", html, fixed = TRUE))
  expect_true(grepl("Records carrying any", html, fixed = TRUE))
  expect_true(grepl("Completeness (", html, fixed = TRUE))
  # d has a blank, produced no finding, and must still be visible as a column.
  expect_true(grepl("Not grouped", html, fixed = TRUE))
})

test_that("block settings are validated", {
  df_Data <- data.frame(x = c(NA, 1:9), y = 1:10)
  for (bad in list(list(min_block_size = 1), list(min_block_events = 1)))
    expect_error(DetectDataAnomalies(df_Data, methods = "missingness", method_args = list(missingness = bad)), "integer")
  expect_error(DetectDataAnomalies(df_Data, methods = "missingness",
    method_args = list(missingness = list(block_similarity = 0))), "block_similarity")
  expect_error(DetectDataAnomalies(df_Data, methods = "missingness",
    method_args = list(missingness = list(max_block_gap = 2))), "max_block_gap")
})
