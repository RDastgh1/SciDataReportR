AnomalyFixture <- function() {
  set.seed(12)
  df_Data <- data.frame(ID = 1:150, x = rnorm(150), y = rnorm(150),
    Group = rep(c("A", "B"), each = 75))
  df_Data$y <- df_Data$x + df_Data$y / 5
  df_Data$y[1] <- 30
  attr(df_Data$x, "label") <- "Measurement"
  attr(df_Data$y, "label") <- "Measurement"
  df_Data
}

LongAnomalyFixture <- function() {
  set.seed(13)
  df_Data <- expand.grid(ID = 1:15, Time = c(0, 1, 3, 6))
  df_Data$x <- rnorm(60) + df_Data$Time + rep(rnorm(15), 4)
  df_Data$Category <- factor(ifelse(df_Data$ID == 1 & df_Data$Time == 6, "B", "A"))
  df_Data[sample(nrow(df_Data)), ]
}

test_that("validation rejects unsafe model specifications and backend overrides", {
  df_Data <- AnomalyFixture()
  expect_error(DetectDataAnomalies(df_Data, c("x", "x")), "unique")
  expect_error(suppressWarnings(DetectDataAnomalies(df_Data, "missing")), "No requested variables are present")
  expect_error(DetectDataAnomalies(df_Data, "x", id_var = "absent"), "not found")
  expect_error(DetectDataAnomalies(df_Data, "ID", id_var = "ID"), "No model variables")
  expect_error(DetectDataAnomalies(df_Data, "x", thresholds = c(conditional = 1)), "thresholds")
  skip_if_not_installed("isotree")
  expect_error(DetectDataAnomalies(df_Data, "x", methods = "isolation", method_args = list(isolation = list(seed = 4))), "wrapper-controlled")
  expect_error(DetectDataAnomalies(df_Data, "x", methods = "isolation", method_args = list(isolation = list(typo = 4))), "wrapper-controlled")
})

test_that("absent screened columns are skipped with a warning rather than aborting the run", {
  df_Data <- AnomalyFixture()
  expect_warning(result <- DetectDataAnomalies(df_Data, c("x", "y", "CD8 % naive", "IFN-b"),
    methods = "conditional"), "not found and skipped")
  # The run completes on the columns that exist, and the two stale names are
  # dropped from the model rather than screened as all-missing columns.
  expect_true(nrow(result$Results) > 0)
  expect_equal(result$Metadata$Variables, c("x", "y"))
  diagnostics <- result$Diagnostics
  expect_true("absent_columns" %in% diagnostics$Issue)
  detail <- diagnostics$Detail[diagnostics$Issue == "absent_columns"]
  expect_match(detail, "CD8 % naive")
  expect_match(detail, "IFN-b")
})

test_that("an absent covariate is skipped while the remaining adjustment stands", {
  df_Data <- LongAnomalyFixture()
  expect_warning(result <- DetectLongitudinalAnomalies(df_Data, "x", id_var = "ID", time_var = "Time",
    covariates = c("Category", "absent"), methods = "change"), "not found and skipped")
  expect_true(nrow(result$Results) > 0)
  expect_equal(result$Metadata$Covariates, "Category")
})

test_that("distance adapter preserves labels, rows, thresholds and caller state", {
  skip_if_not_installed("robustbase")
  df_Data <- AnomalyFixture()
  original <- df_Data
  state <- .Random.seed
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "ID", methods = "distance", thresholds = c(distance = 4))
  expect_identical(df_Data, original)
  expect_identical(.Random.seed, state)
  expect_s3_class(result, "DataAnomalyResult")
  expect_equal(result$Results$.Row, seq_len(150))
  expect_identical(result$ReviewData$x, df_Data$x)
  expect_equal(unname(result$Metadata$Labels[c("x", "y")]), rep("Measurement", 2))
  set.seed(2024)
  fit <- robustbase::covMcd(df_Data[c("x", "y")])
  expected <- sqrt(mahalanobis(df_Data[c("x", "y")], fit$center, fit$cov))
  expect_equal(result$Results$Score, expected)
  expect_equal(result$Results$Flag, expected >= 4)
  expect_equal(result$ParticipantSummary$FlaggedRecords, as.integer(expected >= 4))
})

test_that("LOF and isolation scores agree with their backends", {
  df_Data <- AnomalyFixture()
  skip_if_not_installed("dbscan")
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "local")
  expect_equal(r$Results$Score, dbscan::lof(scale(df_Data[c("x", "y")]), minPts = 5))
  # Score-only methods take the default percentile cutoff; NULL restores unflagged scores.
  expect_equal(r$Results$Flag, r$Results$Score >= stats::quantile(r$Results$Score, 0.99, names = FALSE))
  expect_true(all(is.na(DetectDataAnomalies(df_Data, c("x", "y"), methods = "local",
    flag_percentile = NULL)$Results$Flag)))
  skip_if_not_installed("isotree")
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "isolation", method_args = list(isolation = list(ntrees = 20)))
  fit <- isotree::isolation.forest(df_Data[c("x", "y")], seed = 2024, nthreads = 1, ntrees = 20)
  expect_equal(r$Results$Score, as.numeric(predict(fit, df_Data[c("x", "y")], type = "score")))
})

test_that("mixed LOF uses Gower and stratification never pools missing groups", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("cluster")
  df_Data <- AnomalyFixture()
  r <- DetectDataAnomalies(df_Data, c("x", "Group"), methods = "local")
  expect_equal(r$Results$Score, dbscan::lof(cluster::daisy(df_Data[c("x", "Group")] %>%
    dplyr::mutate(Group = factor(Group)), metric = "gower"), minPts = 5))
  df_Data$Group[1] <- NA
  df_Data$x[2] <- NA
  r <- DetectDataAnomalies(df_Data, c("x", "y"), stratify_by = "Group", methods = "local")
  expect_true(all(r$Results$Status[1:2] == "unavailable"))
  expect_true("missing_stratum" %in% r$Diagnostics$Issue)
  for (s in unique(na.omit(r$Results$Stratum))) {
    i <- which(r$Results$Stratum == s & r$Results$Status == "evaluated")
    expect_equal(r$Results$Rank[i], rank(-r$Results$Score[i], ties.method = "min"))
  }
})

test_that("conditional targets and DDC cells retain backend evidence", {
  df_Data <- AnomalyFixture()
  skip_if_not_installed("outliertree")
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "conditional")
  expect_true(all(r$Results$Status == "evaluated"))
  expect_true(all(is.na(r$Results$Score)))
  expect_true(nrow(r$CellResults) > 0)
  expect_true(all(r$CellResults$Variable %in% c("x", "y")))
  skip_if_not_installed("cellWise")
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "cellwise")
  fit <- cellWise::DDC(as.matrix(df_Data[c("x", "y")]), DDCpars = list(silent = TRUE))
  expect_equal(r$CellResults$Residual, as.numeric(fit$stdResid))
  native <- matrix(FALSE, 150, 2)
  native[fit$indcells] <- TRUE
  expect_equal(r$CellResults$Flag, as.vector(native))
})

test_that("DDC filtering maps back to original rows and variables", {
  skip_if_not_installed("cellWise")
  df_Data <- AnomalyFixture()
  df_Data$x[1:90] <- NA
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "cellwise")
  if (nrow(r$CellResults)) {
    expect_true(all(r$CellResults$Observed == as.character(df_Data$y[r$CellResults$.Row])))
    expect_true(all(r$CellResults$Variable == "y"))
  }
  expect_true(any(grepl("ddc_excluded|backend_failure", r$Diagnostics$Issue)))
})

test_that("longitudinal change scores use real intervals and both endpoints", {
  df_Data <- LongAnomalyFixture()
  r <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time", methods = "change", thresholds = c(change = 2))
  pairs <- r$IntervalResults
  rates <- (df_Data$x[pairs$.Row] - df_Data$x[pairs$PreviousRow]) /
    (df_Data$Time[pairs$.Row] - df_Data$Time[pairs$PreviousRow])
  expect_equal(pairs$Rate, rates)
  expect_equal(pairs$Score, abs(rates - median(rates)) / mad(rates))
  expect_equal(df_Data$ID[pairs$.Row], df_Data$ID[pairs$PreviousRow])
  expect_true(all(df_Data$Time[pairs$.Row] > df_Data$Time[pairs$PreviousRow]))
  for (i in seq_len(nrow(df_Data))) {
    score <- max(pairs$Score[pairs$.Row == i | pairs$PreviousRow == i])
    expect_equal(r$Results$Score[i], score)
  }
  expect_equal(r$Results$.Row, seq_len(60))
  dates <- df_Data
  dates$Time <- as.POSIXct("2020-01-01", tz = "UTC") + dates$Time * 86400
  d <- DetectLongitudinalAnomalies(dates, "x", "ID", "Time", methods = "change")
  expect_equal(d$IntervalResults$Score, pairs$Score)
})

test_that("transitions exclude all data from the scored participant", {
  df_Data <- LongAnomalyFixture()
  r <- DetectLongitudinalAnomalies(df_Data, "Category", "ID", "Time", methods = "transition")
  p <- r$IntervalResults
  for (i in seq_len(nrow(p))) {
    reference <- p[p$Participant != p$Participant[i] & p$PreviousValue == p$PreviousValue[i], ]
    expect_equal(p$ConditioningCount[i], nrow(reference))
    expect_equal(p$Score[i], 1 - mean(reference$Observed == p$Observed[i]))
  }
  expect_true(all(is.na(DetectLongitudinalAnomalies(df_Data, "Category", "ID", "Time",
    methods = "transition", flag_percentile = NULL)$Results$Flag)))
})

test_that("sparse, missing, duplicate and degenerate inputs are explicit", {
  df_Data <- LongAnomalyFixture()
  expect_error(DetectLongitudinalAnomalies(rbind(df_Data, df_Data[1, ]), "x", "ID", "Time"), "Duplicate")
  df_Data$ID[1] <- NA
  df_Data$Time[2] <- NA
  r <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time")
  expect_true(all(r$Results$Status[r$Results$.Row %in% 1:2] == "unavailable"))
  r <- DetectLongitudinalAnomalies(df_Data[1:3, ], "x", "ID", "Time", methods = "change")
  expect_true(all(is.na(r$Results$Score)))
  df_Data <- LongAnomalyFixture()
  df_Data$x <- df_Data$Time
  r <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time", methods = "change")
  expect_true("zero_mad" %in% r$Diagnostics$Issue)
  expect_true(all(is.na(r$Results$Score)))
  skip_if_not_installed("dbscan")
  df_Data$x <- 1
  r <- DetectDataAnomalies(df_Data, "x", methods = "local")
  expect_true(all(r$Results$Status == "unavailable"))
})

test_that("trajectory adapter matches nlme residuals and handles failed fits", {
  skip_if_not_installed("nlme")
  df_Data <- LongAnomalyFixture()
  r <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time", methods = "trajectory")
  cells <- r$CellResults
  expect_equal(nrow(cells), 60)
  sorted <- df_Data[order(df_Data$ID, df_Data$Time), ]
  fit <- nlme::lme(x ~ Time, random = ~ 1 | ID, data = sorted)
  expect_equal(cells$Score, abs(as.numeric(residuals(fit, type = "pearson"))))
  expect_equal(cells$Expected, as.numeric(fitted(fit)))
  failed <- DetectLongitudinalAnomalies(df_Data[1:4, ], "x", "ID", "Time", methods = "trajectory")
  expect_true(all(failed$Results$Status == "unavailable"))
  expect_true("backend_failure" %in% failed$Diagnostics$Issue)
})

test_that("explorers render both classes, empty evidence and escaped text", {
  skip_if_not_installed("reactable")
  skip_if_not_installed("dbscan")
  r <- DetectDataAnomalies(AnomalyFixture(), c("x", "y"), methods = "local")
  before <- r
  html <- htmltools::renderTags(ExploreDataAnomalies(r, Title = "<script>unsafe</script>"))$html
  expect_match(html, "&lt;script&gt;")
  expect_identical(r, before)
  expect_error(ExploreDataAnomalies(r, TopN = 0), "positive")
  skip_if_not_installed("svglite")
  r <- DetectLongitudinalAnomalies(LongAnomalyFixture(), "x", "ID", "Time", methods = "change")
  expect_s3_class(ExploreDataAnomalies(r, TopN = 1), "shiny.tag.list")
  expect_true(nchar(htmltools::renderTags(ExploreDataAnomalies(r, TopN = 1))$html) > 0)
})

test_that("empty inputs, tied scores and missing RNG state remain valid", {
  skip_if_not_installed("dbscan")
  df_Data <- AnomalyFixture()
  r <- DetectDataAnomalies(df_Data[FALSE, ], c("x", "y"), methods = "local")
  expect_equal(nrow(r$Results), 0)
  expect_equal(nrow(r$IntervalResults), 0)
  expect_s3_class(ExploreDataAnomalies(r), "shiny.tag.list")
  df_Data$x <- rep(1:3, 50)
  df_Data$y <- rep(2:4, 50)
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "local")
  expect_equal(r$Results$Rank, rank(-r$Results$Score, ties.method = "min", na.last = "keep"))
  old <- .Random.seed
  on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(DetectDataAnomalies(df_Data, "x", methods = "local"))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("cellwise thresholds override native flags without losing missingness", {
  skip_if_not_installed("cellWise")
  df_Data <- AnomalyFixture()
  r <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "cellwise", thresholds = c(cellwise = 4))
  expect_equal(r$CellResults$Flag, r$CellResults$Score >= 4)
  expect_equal(r$Results$Flag, r$Results$Score >= 4)
})

test_that("categorical conditional screening returns learned conditions", {
  skip_if_not_installed("outliertree")
  df_Data <- data.frame(History = factor(rep(c("A", "B"), each = 1000)))
  df_Data$Response <- df_Data$History
  df_Data$Response[1] <- "B"
  r <- DetectDataAnomalies(df_Data, c("History", "Response"), methods = "conditional",
    method_args = list(conditional = list(min_gain = 0, z_outlier = 4)))
  expect_true(all(r$Results$Status == "evaluated"))
  expect_true(any(r$Results$Flag))
  expect_true(any(lengths(r$CellResults$Conditions) > 0))
})

test_that("transition support and covariate model coverage are explicit", {
  df_Data <- LongAnomalyFixture()
  small <- df_Data[df_Data$ID < 4, ]
  r <- DetectLongitudinalAnomalies(small, "Category", "ID", "Time", methods = "transition")
  expect_true(all(is.na(r$IntervalResults$Score)))
  expect_true("insufficient_support" %in% r$Diagnostics$Issue)
  skip_if_not_installed("nlme")
  df_Data$Age <- rep(seq_len(15), 4)
  df_Data$Age[1] <- NA
  r <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time", methods = "trajectory", covariates = "Age")
  expect_equal(r$Results$Status[1], "unavailable")
  expect_true("excluded_rows" %in% r$Diagnostics$Issue)
  expect_equal(nrow(r$CellResults), 59)
})

test_that("context columns cannot overwrite scores and repeated roles stay unique", {
  df_Data <- AnomalyFixture()
  df_Data$Score <- df_Data$ID
  expect_error(DetectDataAnomalies(df_Data, "x", id_var = "Score"), "conflict")
  skip_if_not_installed("dbscan")
  r <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "Group", stratify_by = "Group", methods = "local")
  expect_false(anyDuplicated(names(r$Results)) > 0)
  expect_equal(nrow(r$Results), nrow(df_Data))
})

test_that("DDC nonfinite residuals are unavailable, not anomaly flags", {
  skip_if_not_installed("cellWise")
  testthat::local_mocked_bindings(DDC = function(X, DDCpars) {
    residuals <- matrix(0, nrow(X), ncol(X))
    residuals[, 1] <- Inf
    list(stdResid = residuals, Xest = X, rowInAnalysis = seq_len(nrow(X)),
      colInAnalysis = seq_len(ncol(X)), indcells = seq_len(nrow(X)))
  }, .package = "cellWise")
  r <- DetectDataAnomalies(AnomalyFixture(), c("x", "y"), methods = "cellwise")
  expect_true(all(is.na(r$CellResults$Score[r$CellResults$Variable == "x"])))
  expect_true(all(is.na(r$CellResults$Flag[r$CellResults$Variable == "x"])))
  expect_true(all(r$Results$Score == 0))
  expect_false(any(r$Results$Flag))
  expect_true("nonfinite_residuals" %in% r$Diagnostics$Issue)
})

test_that("trajectory plots retain full histories without Plotly", {
  result <- DetectLongitudinalAnomalies(LongAnomalyFixture(), c("x", "Category"), "ID", "Time")
  plots <- PlotDataAnomalyTrajectories(result, participants = 1, TopN = 1)
  expect_length(plots, 2)
  expect_s3_class(plots[[1]], "ggplot")
  expect_equal(nrow(plots[[1]]$data), 4)
  expect_error(PlotDataAnomalyTrajectories(result, participants = -1), "participants")
  expect_error(PlotDataAnomalyTrajectories(result, variables = "absent"), "variables")
  if (!nzchar(system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly"))) skip("Plotly JavaScript unavailable")
  html <- htmltools::renderTags(ExploreDataAnomalies(result, TopN = 1, IncludePlots = TRUE))$html
  expect_match(html, "PlotDataAnomalyTrajectories")
})

test_that("zero-inflated count columns do not make conditional screening flag most records", {
  skip_if_not_installed("outliertree")
  set.seed(5)
  n <- 115
  df_Data <- data.frame(ID = paste0("S", seq_len(n)), stringsAsFactors = FALSE)
  for (j in 1:30) df_Data[[paste0("Cont", j)]] <- round(stats::rnorm(n, 100, 15), 1)
  # Counts that are usually zero: the shape that used to flag most of the sample,
  # because a near-zero spread makes any nonzero value look extreme.
  for (j in 1:30) df_Data[[paste0("Count", j)]] <- stats::rpois(n, 0.3)
  df_Data$Cont1[3] <- 900
  Result <- suppressMessages(DetectDataAnomalies(df_Data, methods = "conditional"))
  expect_lt(mean(Result$Results$Flag, na.rm = TRUE), 0.1)
  expect_true(isTRUE(Result$Results$Flag[3]))
  expect_true(any(Result$Diagnostics$Issue == "conditional_discrete_as_categorical"))
  Unfixed <- suppressMessages(DetectDataAnomalies(df_Data, methods = "conditional",
    method_args = list(conditional = list(min_distinct = 0L))))
  expect_gt(mean(Unfixed$Results$Flag, na.rm = TRUE), 0.3)
})

test_that("a backend formatting error keeps the finding instead of aborting the method", {
  skip_if_not_installed("outliertree")
  set.seed(11)
  df_Example <- data.frame(ID = 1:60, A = c(stats::rnorm(59, 10, 1), 500000),
    B = stats::rnorm(60), C = stats::rnorm(60))
  # outliertree can throw while printing its own report for a finding (an
  # invalid decimal count, or a field width beyond sprintf's limit). Detection
  # already succeeded at that point, so the finding must survive.
  result <- testthat::with_mocked_bindings(
    DetectDataAnomalies(df_Example, c("A", "B", "C"), id_var = "ID", methods = "conditional"),
    capture.output = function(...) stop("invalid format '%.NA'; use format %f"),
    .package = "utils")
  expect_true(nrow(result$CellResults) >= 1L)
  expect_true(any(result$CellResults$Variable == "A" & result$CellResults$.Row == 60L))
  expect_false(any(result$Diagnostics$Issue == "backend_failure"))
  expect_true(any(result$Diagnostics$Issue == "conditional_explanation_rebuilt"))
  explanation <- result$CellResults$Explanation[result$CellResults$Variable == "A"][1]
  expect_match(explanation, "unusual within the learned subgroup")
  expect_match(explanation, "could not format its own report")
})

test_that("an extreme magnitude outlier is reported on every run", {
  skip_if_not_installed("outliertree")
  set.seed(12)
  df_Example <- data.frame(ID = 1:80, BMI = round(stats::rnorm(80, 27, 4), 2),
    Age = round(stats::rnorm(80, 50, 10)), Lab = stats::rnorm(80))
  df_Example$BMI[41] <- 276249.07
  for (attempt in 1:3) {
    result <- DetectDataAnomalies(df_Example, c("BMI", "Age", "Lab"), id_var = "ID",
      methods = "conditional")
    expect_true(any(result$CellResults$Variable == "BMI" & result$CellResults$.Row == 41L),
      info = paste("attempt", attempt))
  }
})

test_that("a column name mangled by an encoding round trip is identified, not just skipped", {
  beta <- "IFN-β"
  gamma <- "IFN-γ"
  naive <- "CD8 % naïve"
  # What UTF-8 becomes when it is read back as MacRoman, which is what a CSV
  # opened and resaved in Excel on a Mac produces.
  mangled <- c("IFN-Œ≤", "IFN-Œ≥", "CD8 % na√Øve")
  df_Data <- AnomalyFixture()
  df_Data[[beta]] <- stats::rnorm(nrow(df_Data))
  df_Data[[gamma]] <- stats::rnorm(nrow(df_Data))
  df_Data[[naive]] <- stats::rnorm(nrow(df_Data))
  expect_warning(result <- DetectDataAnomalies(df_Data, c("x", mangled), methods = "conditional"),
    "did you mean")
  detail <- result$Diagnostics$Detail[result$Diagnostics$Issue == "absent_columns"]
  # Beta and gamma differ only in the mangled tail, so a suggestion that resolved
  # them by stripping the non-ASCII characters would have been ambiguous.
  expect_match(detail, beta, fixed = TRUE)
  expect_match(detail, gamma, fixed = TRUE)
  expect_match(detail, naive, fixed = TRUE)
})

test_that("the suggestion works when the data frame carries the mangled name instead", {
  df_Data <- AnomalyFixture()
  df_Data[["IFN-Œ≤"]] <- stats::rnorm(nrow(df_Data))
  expect_warning(DetectDataAnomalies(df_Data, c("x", "IFN-β"), methods = "conditional"),
    "did you mean")
})

test_that("an absent name that is not an encoding problem is reported without a guess", {
  df_Data <- AnomalyFixture()
  expect_warning(result <- DetectDataAnomalies(df_Data, c("x", "nothing_like_this"),
    methods = "conditional"), "nothing_like_this")
  detail <- result$Diagnostics$Detail[result$Diagnostics$Issue == "absent_columns"]
  expect_false(grepl("did you mean", detail, fixed = TRUE))
})

test_that("plain ASCII names are never repaired into a suggestion", {
  suggestions <- SciDataReportR:::ScidrNameSuggestions(c("Weight", "height"), c("weight", "Height", "x"))
  expect_true(all(is.na(suggestions)))
})
