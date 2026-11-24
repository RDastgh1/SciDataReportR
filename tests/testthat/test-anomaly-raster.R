test_that("raster categories separate direction, blanks and conditional findings", {
  expect_equal(ScidrAnomalyRasterCategoryKey("Unexpected missingness", NA, NA, NA), "missing")
  expect_equal(ScidrAnomalyRasterCategoryKey("Distribution outlier", 10, 1, 5), "high")
  expect_equal(ScidrAnomalyRasterCategoryKey("Distribution outlier", 0, 1, 5), "low")
  expect_equal(ScidrAnomalyRasterCategoryKey("Conditional outlier", NA, NA, NA), "conditional")
  expect_equal(ScidrAnomalyRasterCategoryKey("Uncommon transition", NA, NA, NA), "change")
  expect_equal(ScidrAnomalyRasterCategoryKey("Uncommon combination", NA, NA, NA), "other")
  expect_true(all(ScidrAnomalyRasterCategories()$Key %in%
    c("missing", "high", "low", "conditional", "change", "other")))
})

test_that("severity tiers cover every anomaly type exactly once", {
  df_Severity <- ScidrAnomalySeverity()
  expect_setequal(df_Severity$Key, ScidrAnomalyTypes()$Key)
  df_Table <- ScidrAnomalySeverityTable()
  expect_equal(df_Table$Tier, sort(unique(df_Severity$Severity)))
  expect_setequal(unlist(strsplit(df_Table$Types, "; ")), ScidrAnomalyTypes()$Type)
})

test_that("the anomaly map places findings on variable/record cells", {
  skip_if_not_installed("DT")
  set.seed(31)
  df_Example <- data.frame(ID = sprintf("P%03d", 1:80),
    Trig = round(rlnorm(80, 4.5, 0.3), 2), Lab = rnorm(80))
  df_Example$Trig[5] <- 900
  df_Example$Lab[c(3, 9)] <- NA
  result <- DetectDataAnomalies(df_Example, c("Trig", "Lab"), id_var = "ID",
    methods = c("conditional", "missingness"))
  queue <- CreateAnomalyReviewQueue(result)
  display <- ScidrReviewNarratives(queue)
  map <- ScidrAnomalyRasterData(queue, display, names(queue$Sources)[1])
  expect_equal(unlist(map$variables), c("Trig", "Lab"))
  expect_length(map$rows, 80)
  expect_true(length(map$cells) > 0)
  for (cell in map$cells) {
    expect_true(cell$row >= 0 && cell$row < 80)
    expect_true(cell$variable %in% c(0L, 1L))
    expect_true(cell$category >= 1 && cell$category <= nrow(ScidrAnomalyRasterCategories()))
  }
  high <- Filter(function(cell) cell$variable == 0L && cell$row == 4L, map$cells)
  expect_length(high, 1)
  expect_equal(ScidrAnomalyRasterCategories()$Key[high[[1]]$category], "high")
  # Two cells cannot describe the same value twice.
  keys <- vapply(map$cells, function(cell) paste(cell$row, cell$variable), character(1))
  expect_false(anyDuplicated(keys) > 0)
})

test_that("record-level longitudinal findings are typed as changes, not combinations", {
  expect_equal(ScidrAnomalyTypeKey("change", NA, FALSE, FALSE), "change")
  expect_equal(ScidrAnomalyTypeKey("trajectory", NA, FALSE, FALSE), "change")
  expect_equal(ScidrAnomalyTypeKey("transition", NA, FALSE, FALSE), "transition")
  # Isolation screening keeps the exploratory whole-record type.
  expect_equal(ScidrAnomalyTypeKey("isolation", NA, FALSE, FALSE), "combination")
  set.seed(3)
  n <- 60
  df_Example <- data.frame(ID = rep(sprintf("P%02d", 1:30), each = 2), Visit = rep(1:2, 30),
    Trig = round(rlnorm(n, 4.5, 0.3), 2), Lab = rnorm(n))
  df_Example$Trig[7] <- 999
  result <- DetectLongitudinalAnomalies(df_Example, c("Trig", "Lab"), id_var = "ID", time_var = "Visit")
  queue <- CreateAnomalyReviewQueue(result)
  display <- ScidrReviewNarratives(queue)
  expect_true(nrow(display$Issues) > 0)
  expect_true(all(display$Issues$FindingType %in%
    c("Unexpected change over time", "Uncommon transition")))
  expect_true(all(display$Issues$ReviewTier == 5L))
})
