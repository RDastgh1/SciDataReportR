DateAnomalyFixture <- function() {
  df_Data <- data.frame(ID = seq_len(40), Birth = as.Date("1950-01-01") + seq_len(40) * 365,
    Visit = as.Date("2020-01-01") + seq_len(40))
  attr(df_Data$Birth, "label") <- "Date of birth"
  attr(df_Data$Visit, "label") <- "Assessment date"
  df_Data
}

test_that("date screening preserves inputs and excludes selected ID from features", {
  df_Data <- DateAnomalyFixture()
  df_Data$Visit[1] <- as.Date("2099-01-01")
  original <- df_Data
  r <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID", methods = "dates")
  expect_identical(df_Data, original)
  expect_identical(r$ReviewData$Visit, original$Visit)
  expect_equal(r$Metadata$Variables, c("Birth", "Visit"))
  expect_equal(r$Metadata$RequestedVariables, names(df_Data))
  expect_true("id_excluded_from_model" %in% r$Diagnostics$Issue)
  expect_equal(nrow(r$DateResults), 120)
  expect_true(all(r$DateResults$ID == r$DateResults$.Row))
  expect_true(any(r$DateResults$Flag[r$DateResults$.Row == 1] %in% TRUE))
  expect_equal(r$Results$Flag[1], TRUE)
  expect_true(any(grepl("2099-01-01", r$DateResults$Observed)))
  expect_equal(r$Metadata$DateConversions$Visit$ReferenceDate, "2020-01-03")
})

test_that("signed intervals follow selected order and explicit pairs restrict checks", {
  df_Data <- DateAnomalyFixture()
  df_Data$Followup <- df_Data$Visit + 100
  df_Data$Followup[1] <- df_Data$Visit[1] + 1000
  r <- DetectDataAnomalies(df_Data, c("Visit", "Followup", "Birth"), methods = "dates",
    date_pairs = matrix(c("Followup", "Visit"), ncol = 2))
  p <- r$DateResults[r$DateResults$CheckType == "interval", ]
  expect_equal(nrow(p), 40)
  expect_true(all(p$Variable == "Visit" & p$OtherVariable == "Followup"))
  expect_equal(p$ValueDays, as.numeric(df_Data$Followup - df_Data$Visit))
  expect_equal(p$Rule[1], "learned_majority")
  expect_true(p$Flag[1])
  expect_false(any(p$Flag[-1]))
  reverse <- DetectDataAnomalies(df_Data, c("Followup", "Visit"), methods = "dates")
  p <- reverse$DateResults[reverse$DateResults$CheckType == "interval", ]
  expect_true(all(p$ValueDays < 0))
  expect_false(any(p$Flag[-1]))
})

test_that("Tukey bounds agree with independent date calculations", {
  df_Data <- DateAnomalyFixture()
  df_Data$Visit[40] <- df_Data$Visit[40] + 10000
  r <- DetectDataAnomalies(df_Data, "Visit", methods = "dates")
  x <- as.numeric(df_Data$Visit - min(df_Data$Visit))
  q <- quantile(x, c(.25, .75), names = FALSE)
  lower <- q[1] - 3 * diff(q)
  upper <- q[2] + 3 * diff(q)
  expect_equal(r$DateResults$LowerFenceDays, rep(lower, 40))
  expect_equal(r$DateResults$UpperFenceDays, rep(upper, 40))
  expect_equal(r$DateResults$Flag, x < lower | x > upper)
  expect_equal(r$DateResults$Score, pmax(lower - x, x - upper, 0) / max(diff(q), 1))
  expect_equal(r$DateResults$UpperDate[1], as.character(min(df_Data$Visit) + upper))
})

test_that("POSIXct uses elapsed days and preserves time zones", {
  df_Data <- data.frame(ID = seq_len(40),
    Start = as.POSIXct("2024-03-01 12:00:00", tz = "America/New_York") + seq_len(40) * 86400)
  df_Data$End <- df_Data$Start + 3600
  df_Data$End[1] <- df_Data$Start[1] + 7200
  r <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID", methods = "dates")
  expect_identical(r$ReviewData$Start, df_Data$Start)
  expect_identical(attr(r$ReviewData$Start, "tzone"), "America/New_York")
  p <- r$DateResults[r$DateResults$CheckType == "interval", ]
  expect_equal(p$ValueDays, c(2, rep(1, 39)) / 24)
  expect_true(p$Flag[1])
  expect_true(grepl("EST|EDT", p$Observed[1]))
})

test_that("missing dates, insufficient strata and ambiguous zero IQR remain unavailable", {
  df_Data <- DateAnomalyFixture()
  df_Data$Birth[1] <- NA
  df_Data$Visit[] <- as.Date(NA)
  r <- DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates")
  expect_equal(r$Results$Status[1], "unavailable")
  expect_true(all(r$DateResults$Status[r$DateResults$Variable == "Visit"] == "unavailable"))
  df_Data$Group <- rep(1:4, each = 10)
  r <- DetectDataAnomalies(df_Data, "Birth", stratify_by = "Group", methods = "dates")
  expect_true(all(r$Results$Status == "unavailable"))
  df_Data$Visit <- as.Date("2020-01-01") + c(rep(0, 32), rep(1, 8))
  r <- DetectDataAnomalies(df_Data, "Visit", methods = "dates")
  expect_true(all(r$Results$Status == "unavailable"))
  r <- DetectDataAnomalies(df_Data, "Visit", methods = "dates", method_args = list(dates = list(majority_fraction = .8)))
  expect_equal(sum(r$Results$Flag), 8)
})

test_that("date controls validate pairs and cannot silently override native rules", {
  df_Data <- DateAnomalyFixture()
  expect_error(DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates", date_pairs = c("Birth", "Visit")), "matrix")
  expect_error(DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates", date_pairs = matrix(c("Birth", "Birth"), ncol = 2)), "itself")
  pairs <- rbind(c("Birth", "Visit"), c("Visit", "Birth"))
  expect_error(DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates", date_pairs = pairs), "duplicate")
  expect_error(DetectDataAnomalies(df_Data, "Birth", methods = "dates", thresholds = c(dates = 1)), "thresholds")
  expect_error(DetectDataAnomalies(df_Data, "Birth", methods = "dates", method_args = list(dates = list(min_observations = 1))), "Invalid date settings")
  expect_error(DetectDataAnomalies(df_Data, "Birth", methods = "dates", method_args = list(dates = list(iqr_multiplier = NA_real_))), "finite numeric")
  r <- DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates", date_pairs = matrix(character(), ncol = 2))
  expect_true(all(r$DateResults$CheckType == "column"))
})

test_that("many pairs warn and retain every comparison", {
  df_Data <- as.data.frame(setNames(rep(list(as.Date("2020-01-01") + seq_len(20)), 16), paste0("d", seq_len(16))))
  expect_warning(r <- DetectDataAnomalies(df_Data, names(df_Data), methods = "dates"), "120 pairs")
  expect_equal(nrow(r$DateResults), 20 * (16 + 120))
})

test_that("date absence and zero-row input are safe without optional backends", {
  r <- DetectDataAnomalies(data.frame(x = 1:30), "x", methods = "dates")
  expect_true(all(r$Results$Status == "not applicable"))
  expect_equal(nrow(r$DateResults), 0)
  df_Data <- DateAnomalyFixture()[FALSE, ]
  r <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID", methods = "dates")
  expect_equal(nrow(r$Results), 0)
  expect_equal(nrow(r$DateResults), 0)
})

test_that("existing methods use local numeric date copies", {
  skip_if_not_installed("dbscan")
  df_Data <- DateAnomalyFixture()
  df_Data$Visit[1] <- df_Data$Visit[1] + 1000
  r <- DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "local")
  x <- data.frame(Birth = as.numeric(df_Data$Birth - min(df_Data$Birth)),
    Visit = as.numeric(df_Data$Visit - min(df_Data$Visit)))
  expect_equal(r$Results$Score, dbscan::lof(scale(x), minPts = 5))
  expect_identical(r$ReviewData$Birth, df_Data$Birth)
  expect_true(all(r$Results$Status == "evaluated"))
  skip_if_not_installed("isotree")
  r <- DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "isolation", method_args = list(isolation = list(ntrees = 20)))
  fit <- isotree::isolation.forest(x, seed = 2024, nthreads = 1, ntrees = 20)
  expect_equal(r$Results$Score, as.numeric(predict(fit, x, type = "score")))
})

test_that("conditional date findings use readable dates and original values", {
  skip_if_not_installed("outliertree")
  df_Data <- data.frame(Visit = as.Date("2020-01-01") + seq_len(150))
  df_Data$Visit[1] <- as.Date("2099-01-01")
  r <- DetectDataAnomalies(df_Data, "Visit", methods = "conditional")
  expect_true(r$Results$Flag[1])
  expect_match(r$CellResults$Observed[1], "2099-01-01")
  expect_match(r$CellResults$Explanation[1], "2099-01-01")
})

test_that("date explorer shows date evidence and supports old results", {
  skip_if_not_installed("reactable")
  df_Data <- DateAnomalyFixture()
  df_Data$Birth[1] <- as.Date("1800-01-01")
  r <- DetectDataAnomalies(df_Data, c("Birth", "Visit"), methods = "dates")
  html <- htmltools::renderTags(ExploreDataAnomalies(r))$html
  expect_match(html, "Implausible date")
  expect_match(html, "Date of birth")
  r$DateResults <- NULL
  expect_s3_class(ExploreDataAnomalies(r), "shiny.tag.list")
})

test_that("missing and duplicate IDs do not pool within-row date comparisons", {
  df_Data <- DateAnomalyFixture()
  df_Data$ID <- rep("same", 40)
  df_Data$ID[1] <- NA
  r <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID", methods = "dates")
  expect_equal(nrow(r$Results), 40)
  expect_equal(nrow(r$DateResults), 120)
  expect_true(all(c("repeated_id", "missing_id") %in% r$Diagnostics$Issue))
  p <- r$DateResults[r$DateResults$CheckType == "interval", ]
  expect_equal(p$ValueDays, as.numeric(df_Data$Visit - df_Data$Birth))
})

test_that("date models keep backend cell evidence readable", {
  skip_if_not_installed("cellWise")
  set.seed(18)
  df_Data <- data.frame(First = as.Date("2020-01-01") + sample.int(100, 80, replace = TRUE),
    Second = as.Date("2021-01-01") + sample.int(100, 80, replace = TRUE))
  df_Data$Second[1] <- as.Date("2099-01-01")
  r <- DetectDataAnomalies(df_Data, names(df_Data), methods = "cellwise")
  expect_true(all(r$Results$Status == "evaluated"))
  expect_true(all(grepl("^[0-9]{4}-", r$CellResults$Observed)))
  expect_true(all(grepl("^[0-9]{4}-", r$CellResults$ExpectedDisplay)))
  expect_true(any(grepl("2099-01-01", r$CellResults$Observed)))
})

test_that("default methods include dates and conditional date conditions are readable", {
  skip_if_not_installed("outliertree")
  skip_if_not_installed("isotree")
  df_Data <- DateAnomalyFixture()
  r <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID")
  expect_equal(r$Metadata$Methods, c("conditional", "isolation", "dates", "missingness"))
  expect_true(all(c("conditional", "isolation", "dates") %in% r$Summary$Method))
  df_Data <- data.frame(Visit = as.Date("2020-01-01") + rep(c(0, 1000), each = 150))
  set.seed(1)
  df_Data$Value <- c(rnorm(150, 10, .1), rnorm(150, 30, .1))
  df_Data$Value[1] <- 20
  r <- DetectDataAnomalies(df_Data, c("Visit", "Value"), methods = "conditional")
  expect_true(r$Results$Flag[1])
  expect_match(r$CellResults$Explanation[r$CellResults$.Row == 1], "2020-01-01")
})

test_that("constant dates evaluate normally and nonfinite dates stay unavailable", {
  df_Data <- data.frame(Visit = rep(as.Date("2020-01-01"), 25))
  df_Data$Visit[1] <- structure(Inf, class = "Date")
  r <- DetectDataAnomalies(df_Data, "Visit", methods = "dates")
  expect_equal(r$Results$Status[1], "unavailable")
  expect_true(all(r$Results$Status[-1] == "evaluated"))
  expect_false(any(r$Results$Flag %in% TRUE))
  expect_equal(unique(r$DateResults$MinimumDate), "2020-01-01")
})
