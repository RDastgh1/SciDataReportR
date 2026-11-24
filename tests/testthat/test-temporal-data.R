TemporalDataFixture <- function() {
  set.seed(11)
  n <- 150L
  hours <- c(runif(n - 2, 22, 26) %% 24, 12, 13)
  df_Data <- data.frame(clock = hms::as_hms(hours * 3600), duration = as.difftime(rexp(n, 1/60), units = "secs"),
    x = rnorm(n), y = rnorm(n), group = factor(rep(c("A", "B"), length.out = n)))
  codebook <- data.frame(Variable = c("clock", "duration"), Type = c("ClOcK", "Duration"),
    Label = c("Assessment clock", "Task duration"), Category = c("Timing", "Performance"))
  list(Data = df_Data, Codebook = codebook)
}

test_that("temporal intake preserves typed values/classes/units and attaches meaning", {
  skip_if_not_installed("hms")
  df_Data <- data.frame(clock = hms::as_hms(c(86340, 60, NA)), duration = as.difftime(c(1, 2, NA), units = "mins"))
  before <- df_Data
  codebook <- data.frame(Variable = names(df_Data), Type = c("CLOCK", "Duration"))
  result <- RevalueData(df_Data, codebook)
  expect_identical(df_Data, before)
  for (variable in names(df_Data)) {
    expect_identical(class(result$RevaluedData[[variable]]), class(df_Data[[variable]]))
    expect_identical(attr(result$RevaluedData[[variable]], "units"), attr(df_Data[[variable]], "units"))
    expect_equal(as.numeric(result$RevaluedData[[variable]], units = "secs"), as.numeric(df_Data[[variable]], units = "secs"))
  }
  expect_equal(attr(result$RevaluedData$clock, "scidr_type"), "clock")
  expect_equal(attr(result$RevaluedData$duration, "scidr_type"), "duration")
  expect_true(all(result$TimeMetadata$Status == "Converted"))
  override <- RevalueData(df_Data, codebook, time_types = c(clock = "duration"))
  expect_equal(attr(override$RevaluedData$clock, "scidr_type"), "duration")
  expect_equal(override$TimeMetadata$Source[1], "override")
})

test_that("explicit numeric/text conversions handle units, missing codes and atomic failure", {
  skip_if_not_installed("hms")
  df_Data <- data.frame(clock = c(0.5, 0, NA), duration = c("25:01:02.5", "-00:01:00", "missing"), short = c("01:30.5", "00:20", NA))
  codebook <- data.frame(Variable = names(df_Data), Type = c("clock", "duration", "duration"),
    TimeUnit = c("days", NA, NA), TimeFormat = c(NA, "HMS", "MS"), MissingCode = c(NA, "missing", NA))
  result <- RevalueData(df_Data, codebook)
  expect_equal(as.numeric(result$RevaluedData$clock), c(43200, 0, NA))
  expect_equal(as.numeric(result$RevaluedData$duration), c(90062.5, -60, NA))
  expect_equal(as.numeric(result$RevaluedData$short), c(90.5, 20, NA))
  expect_s3_class(result$RevaluedData$duration, "difftime")
  expect_equal(result$TimeMetadata$OutputUnit, rep("secs", 3))
  bad <- data.frame(clock = c("23:59", "25:00"), duration = c("01:90", "bad"), value = 1:2)
  cb <- data.frame(Variable = names(bad), Type = c("clock", "duration", "continuous"), TimeFormat = c("HM", "MS", NA))
  original <- bad
  expect_error(RevalueData(bad, cb), "clock.*rows: 2")
  warned <- RevalueData(bad, cb, on_error = "warn")
  expect_identical(warned$RevaluedData$clock, original$clock)
  expect_identical(warned$RevaluedData$duration, original$duration)
  expect_equal(nrow(warned$errors), 2)
  expect_true(all(warned$TimeMetadata$Status == "Failed"))
  expect_identical(bad, original)
  expect_error(RevalueData(data.frame(t = 1:3), data.frame(Variable = "t", Type = "duration")), "TimeUnit")
  expect_error(RevalueData(data.frame(t = c("12:30", "13:30")), data.frame(Variable = "t", Type = "duration")), "TimeFormat")
  cb_map <- data.frame(Variable = "t", Type = "duration", Recode = "Yes", Code = "1=One", TimeUnit = "secs")
  expect_error(RevalueData(data.frame(t = 1:3), cb_map), "Categorical value-code")
})

test_that("declarations include typed duration defaults and warn on ambiguous hms", {
  skip_if_not_installed("hms")
  df_Data <- data.frame(t = hms::as_hms(1:30), d = as.difftime(1:30, units = "mins"), x = seq_len(30))
  # Ambiguous hms warns and continues: the column is compared as plain seconds
  # rather than blocking the whole analysis.
  expect_warning(Ambiguous <- DetectDataAnomalies(df_Data, methods = "dates"),
    "Ambiguous hms columns compared as plain seconds")
  expect_setequal(Ambiguous$Metadata$Variables, names(df_Data))
  expect_false("t" %in% Ambiguous$Metadata$TimeConversions$Variable)
  result <- DetectDataAnomalies(df_Data, methods = "dates", time_types = c(t = "clock"))
  expect_setequal(result$Metadata$Variables, names(df_Data))
  expect_equal(result$Metadata$TimeConversions$TimeType, c("clock", "duration"))
  expect_equal(result$Metadata$TimeConversions$Source, c("override", "difftime class"))
  expect_true(all(result$Results$Status == "not applicable"))
  expect_error(DetectDataAnomalies(df_Data, time_types = c(t = "time")), "clock or duration")
  expect_error(DetectDataAnomalies(df_Data, time_types = c(unknown = "clock")), "existing columns")
  declared <- RevalueData(df_Data, data.frame(Variable = "t", Type = "duration"))$RevaluedData
  detected <- DetectDataAnomalies(declared, methods = "dates")
  expect_equal(detected$Metadata$TimeConversions$TimeType, c("duration", "duration"))
})

test_that("circular encoding and Gower distance respect midnight and physical units", {
  df_Data <- data.frame(clock = c(86340, 60, 43200), value = c(1, 1, 1), group = factor(c("A", "A", "A")))
  declarations <- data.frame(Variable = "clock", TimeType = "clock")
  map <- ScidrTimeFeatureMap(names(df_Data), declarations)
  encoded <- ScidrTimeEncode(df_Data, map)
  distance <- as.matrix(stats::dist(encoded[unlist(map["clock"], use.names = FALSE)]))
  expect_lt(distance[1, 2], distance[1, 3] / 100)
  circular <- as.matrix(ScidrTimeGower(df_Data[c("clock", "group")], declarations))
  expect_equal(circular[1, 2], (120 / 43200) / 2)
  expect_equal(circular[2, 3], (43140 / 43200) / 2)
  expect_equal(ScidrTimeSeconds(as.difftime(2, units = "hours"), data.frame()), 7200)
  for (unit in c("secs", "mins", "hours", "days", "weeks")) {
    amount <- c(secs = 604800, mins = 10080, hours = 168, days = 7, weeks = 1)[unit]
    declaration <- data.frame(TimeUnit = unit)
    expect_equal(ScidrTimeSeconds(unname(amount), declaration), 604800)
  }
})

test_that("native models evaluate temporal inputs and preserve originals and RNG", {
  skip_if_not_installed("hms")
  for (backend in c("outliertree", "isotree", "dbscan", "robustbase", "cellWise")) skip_if_not_installed(backend)
  fixture <- TemporalDataFixture()
  original <- fixture$Data
  before <- .Random.seed
  for (method in c("conditional", "isolation", "local", "distance", "cellwise")) {
    result <- DetectDataAnomalies(fixture$Data, codebook = fixture$Codebook, methods = method,
      method_args = if (method == "isolation") list(isolation = list(ntrees = 20)) else list(),
      explanation_args = if (method == "isolation") list(rows = c(1, 150)) else list())
    expect_true(all(result$Results$Status == "evaluated"))
    expect_identical(.Random.seed, before)
    expect_identical(fixture$Data, original)
    expect_identical(result$ReviewData$clock, original$clock)
    expect_identical(result$ReviewData$duration, original$duration)
    if (nrow(result$CellResults)) {
      expect_true(all(result$CellResults$Variable %in% names(original)))
      if (method == "cellwise") {
        expect_false(anyDuplicated(result$CellResults[c(".Row", "Variable")]) > 0)
        clock <- result$CellResults[result$CellResults$Variable == "clock", ]
        expect_true(all(is.na(clock$Expected)))
        expect_true(any(lengths(clock$NativeComponents) > 0))
      }
    }
    queue <- CreateAnomalyReviewQueue(result, top_n = 2)
    if (nrow(queue$Issues)) {
      narratives <- ScidrReviewNarratives(queue)
      expect_false(any(grepl("Recorded clock time.*above.*boundary", narratives$Issues$WhyFlagged)))
    }
  }
})

test_that("sensitivity replaces complete clocks and matches native predictions", {
  skip_if_not_installed("isotree")
  skip_if_not_installed("hms")
  fixture <- TemporalDataFixture()
  result <- DetectDataAnomalies(fixture$Data, codebook = fixture$Codebook, methods = "isolation", seed = 27,
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(rows = 150))
  metadata <- result$ExplanationMetadata$Strata[[1]]
  df_Source <- fixture$Data[metadata$SourceVariables]
  for (variable in c("clock", "duration")) df_Source[[variable]] <- as.numeric(df_Source[[variable]], units = "secs")
  fit <- isotree::isolation.forest(ScidrTimePredictData(df_Source, metadata$ModelAdapter), ntrees = 20, nthreads = 1, seed = 27)
  df_Changed <- df_Source[rep(150, length(metadata$ReferenceRows)), , drop = FALSE]
  df_Changed$clock <- df_Source$clock[metadata$ReferenceRows]
  predicted <- predict(fit, ScidrTimePredictData(df_Changed, metadata$ModelAdapter), type = "score", nthreads = 1)
  native <- result$Results$Score[150]
  effect <- result$FeatureSensitivity[result$FeatureSensitivity$Variable == "clock", ]
  expect_equal(effect$ScoreDecrease, native - mean(predicted))
  expect_equal(nrow(result$FeatureSensitivity), ncol(fixture$Data))
  expect_match(result$FeatureSensitivity$Observed[result$FeatureSensitivity$Variable == "duration"], "seconds")
  expect_true(all(result$FeatureSensitivity$Variable %in% names(fixture$Data)))
})

test_that("SHAP explains source clocks through the exact encoded prediction wrapper", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  skip_if_not_installed("fastshap")
  fixture <- TemporalDataFixture()
  result <- DetectDataAnomalies(fixture$Data, codebook = fixture$Codebook, methods = "isolation", seed = 27,
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(method = "shap", rows = 150, nsim = 4, repeats = 2))
  metadata <- result$ExplanationMetadata$Strata[[1]]
  df_Source <- fixture$Data[metadata$SourceVariables]
  for (variable in c("clock", "duration")) df_Source[[variable]] <- as.numeric(df_Source[[variable]], units = "secs")
  fit <- isotree::isolation.forest(ScidrTimePredictData(df_Source, metadata$ModelAdapter), ntrees = 20, nthreads = 1, seed = 27)
  Predict <- function(object, newdata) as.numeric(predict(object, ScidrTimePredictData(newdata, metadata$ModelAdapter), type = "score", nthreads = 1))
  background <- df_Source[metadata$BackgroundRows, , drop = FALSE]
  baseline <- mean(Predict(fit, background))
  values <- vapply(metadata$Seeds[["150"]], function(seed) {
    set.seed(seed)
    as.numeric(fastshap::explain(fit, X = background, newdata = df_Source[150, , drop = FALSE], nsim = 4, adjust = TRUE,
      baseline = baseline, pred_wrapper = Predict, parallel = FALSE)[1, names(df_Source)])
  }, numeric(ncol(df_Source)))
  expect_equal(result$FeatureContributions$Contribution, rowMeans(values), tolerance = 1e-10)
  expect_setequal(result$FeatureContributions$Variable, names(df_Source))
  expect_lt(abs(result$ExplanationSummary$ReconstructionResidual[150]), 1e-6)
})

test_that("invalid clocks are diagnosed rather than wrapped and parsing failures stay unavailable", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  df_Data <- data.frame(clock = hms::as_hms(c(49*3600, seq_len(39)*100)), duration = as.difftime(seq_len(40), units = "secs"))
  before <- df_Data
  result <- DetectDataAnomalies(df_Data, methods = "isolation", time_types = c(clock = "clock"), method_args = list(isolation = list(ntrees = 20)))
  expect_true("invalid_temporal_value" %in% result$Diagnostics$Issue)
  expect_identical(result$ReviewData$clock, before$clock)
  expect_identical(df_Data, before)
  expect_equal(result$Diagnostics$Rows[[which(result$Diagnostics$Issue == "invalid_temporal_value")]], 1L)
  parsed <- DetectDataAnomalies(data.frame(clock = c("oops", "00:10"), x = c(1, 2)), methods = "dates",
    time_types = c(clock = "clock"), time_formats = c(clock = "HM"))
  expect_equal(parsed$Metadata$TimeConversions$Status, "Partially available")
  expect_equal(parsed$Metadata$TimeConversions$UnavailableRows[[1]], 1L)
  expect_true("temporal_parse_failure" %in% parsed$Diagnostics$Issue)
})

test_that("clock regions and temporal workbook evidence are readable and protected", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  skip_if_not_installed("openxlsx")
  region <- ScidrTimeRegion("cos", ">", cos(2*pi/24))
  expect_equal(nrow(region$Intervals), 2L)
  expect_match(region$Text, "00:00:00", fixed = TRUE)
  fixture <- TemporalDataFixture()
  result <- DetectDataAnomalies(fixture$Data, methods = "isolation", codebook = fixture$Codebook,
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(top_n = 2))
  queue <- CreateAnomalyReviewQueue(result, top_n = 2)
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  plot <- ScidrReviewPlotSources(queue)[[1]]
  expect_true(plot$clocks$clock)
  expect_true(plot$numeric$duration)
  expect_false(plot$numeric$clock)
  expect_equal(unlist(plot$values$duration), as.numeric(fixture$Data$duration))
  expect_match(ScidrReviewPlotJS(), "barpolar", fixed = TRUE)
  raw <- queue
  ExploreDataAnomalies(queue, codebook = transform(fixture$Codebook, Type = "duration"), IncludePlots = FALSE)
  expect_identical(queue, raw)
})

test_that("strata, duplicate labels and constant clock coordinates retain source mapping", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  set.seed(6)
  df_Data <- data.frame(ID = rep(c("001", "001", NA), 20), site = rep(c("A", "B"), each = 30),
    clock = hms::as_hms(rep(c(0, 43200), 30)), fixed = hms::as_hms(rep(600, 60)),
    duration = as.difftime(rnorm(60), units = "days"))
  codebook <- data.frame(Variable = c("clock", "fixed", "duration"), Type = c("clock", "clock", "duration"),
    Label = rep("Same label", 3), Category = c("Timing", "Timing", "Task"))
  result <- DetectDataAnomalies(df_Data, methods = "isolation", id_var = "ID", stratify_by = "site", codebook = codebook,
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(rows = c(1, 31)))
  expect_true(all(result$Results$Status == "evaluated"))
  expect_false("fixed" %in% result$FeatureSensitivity$Variable)
  expect_true("constant_clock_components" %in% result$Diagnostics$Issue)
  for (metadata in result$ExplanationMetadata$Strata) {
    expect_length(unique(df_Data$site[metadata$ReferenceRows]), 1)
    expect_length(metadata$ModelAdapter$FeatureMap$clock, 1)
  }
  expect_setequal(result$FeatureSensitivity$.Row, c(1, 31))
  expect_identical(result$ReviewData$ID, df_Data$ID)
})

test_that("explicit numeric units reproduce typed temporal native scores", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  fixture <- TemporalDataFixture()
  raw <- fixture$Data
  raw$clock <- as.numeric(raw$clock) / 3600
  raw$duration <- as.numeric(raw$duration) / 60
  book <- fixture$Codebook
  book$TimeUnit <- c("hours", "mins")
  typed <- DetectDataAnomalies(fixture$Data, methods = "isolation", codebook = fixture$Codebook,
    method_args = list(isolation = list(ntrees = 20)), explain_isolation = FALSE)
  numeric <- DetectDataAnomalies(raw, methods = "isolation", codebook = book,
    method_args = list(isolation = list(ntrees = 20)), explain_isolation = FALSE)
  expect_equal(typed$Results$Score, numeric$Results$Score, tolerance = 1e-12)
  expect_identical(numeric$ReviewData$clock, raw$clock)
  expect_error(RevalueData(data.frame(t = 1:3), data.frame(Variable = c("t", "t"), Type = c("clock", "duration"))), "Conflicting temporal Type")
})

test_that("temporal model semantics bind snapshots while legacy identity stays unchanged", {
  skip_if_not_installed("hms")
  df_Data <- data.frame(t = hms::as_hms(seq_len(30)))
  clock <- DetectDataAnomalies(df_Data, methods = "dates", time_types = c(t = "clock"))
  duration <- DetectDataAnomalies(df_Data, methods = "dates", time_types = c(t = "duration"))
  expect_identical(clock$Results, duration$Results)
  expect_false(identical(ScidrReviewSnapshot(clock), ScidrReviewSnapshot(duration)))
  plain <- DetectDataAnomalies(data.frame(x = 1:30), methods = "dates")
  old <- ScidrReviewHash(list(ReviewData = plain$ReviewData, Results = plain$Results, CellResults = plain$CellResults,
    DateResults = plain$DateResults, IntervalResults = plain$IntervalResults, Variables = plain$Metadata$Variables,
    Methods = plain$Metadata$Methods, Thresholds = plain$Metadata$Thresholds))
  expect_identical(ScidrReviewSnapshot(plain), old)
})

test_that("temporal inputs coexist with retained dates and timestamps", {
  skip_if_not_installed("hms")
  skip_if_not_installed("isotree")
  fixture <- TemporalDataFixture()
  df_Data <- fixture$Data[1:40, ]
  df_Data$visit <- as.Date("2020-01-01") + seq_len(40)
  df_Data$stamp <- as.POSIXct("2020-01-01 12:00:00", tz = "America/New_York") + seq_len(40) * 86400
  result <- DetectDataAnomalies(df_Data, methods = "isolation", codebook = fixture$Codebook, seed = 27,
    method_args = list(isolation = list(ntrees = 20)), explanation_args = list(rows = 1))
  metadata <- result$ExplanationMetadata$Strata[[1]]
  df_Source <- df_Data[metadata$SourceVariables]
  for (variable in c("clock", "duration")) df_Source[[variable]] <- as.numeric(df_Source[[variable]], units = "secs")
  for (variable in c("visit", "stamp")) df_Source[[variable]] <- ScidrAnomalyDateDays(df_Source[[variable]]) - result$Metadata$DateConversions[[variable]]$ReferenceDays
  fit <- isotree::isolation.forest(ScidrTimePredictData(df_Source, metadata$ModelAdapter), seed = 27, ntrees = 20, nthreads = 1)
  expect_equal(result$Results$Score, as.numeric(predict(fit, ScidrTimePredictData(df_Source, metadata$ModelAdapter), type = "score", nthreads = 1)))
  expect_identical(result$ReviewData$visit, df_Data$visit)
  expect_identical(result$ReviewData$stamp, df_Data$stamp)
  expect_identical(attr(result$ReviewData$stamp, "tzone"), "America/New_York")
  expect_equal(ScidrTimeSeconds("1:2:3.5", data.frame(TimeFormat = "HMS", TimeType = "clock")), 3723.5)
  excluded <- ScidrTimeRegion("cos", "!=", 0)
  expect_match(excluded$Text, "except", fixed = TRUE)
  expect_equal(length(excluded$ExcludedPoints), 2)
})

test_that("an undeclared hms column is compared as seconds instead of stopping the backends", {
  skip_if_not_installed("outliertree")
  skip_if_not_installed("isotree")
  set.seed(5)
  n <- 80
  seconds <- c(stats::rnorm(n - 1, 3600, 400), 20000)
  df_Data <- data.frame(ID = sprintf("P%02d", 1:n),
    score = c(stats::rnorm(n - 1, 50, 8), 200), age = stats::rnorm(n, 55, 10))
  df_Data$elapsed <- hms::as_hms(seconds)

  Run <- function(...) suppressWarnings(suppressMessages(
    DetectDataAnomalies(df_Data, id_var = "ID", methods = c("conditional", "isolation"), ...)))
  result <- Run()

  # The regression: outliertree and isotree reject an hms column outright, which
  # aborted both methods and returned no findings at all from either.
  expect_false(any(result$Diagnostics$Issue == "backend_failure"))
  expect_false(any(result$Results$Status == "unavailable"))
  expect_gt(sum(result$Results$Flag %in% TRUE), 0L)

  # The warning promises plain seconds, so the model must see exactly that: the
  # same findings a numeric seconds column produces.
  df_Numeric <- df_Data
  df_Numeric$elapsed <- as.numeric(seconds)
  numeric_result <- suppressWarnings(suppressMessages(
    DetectDataAnomalies(df_Numeric, id_var = "ID", methods = c("conditional", "isolation"))))
  expect_identical(result$Results$Flag, numeric_result$Results$Flag)

  # The conversion is recorded where a reader can see it, and names the fix.
  diagnostic <- result$Diagnostics[result$Diagnostics$Issue == "ambiguous_temporal_seconds", ]
  expect_equal(nrow(diagnostic), 1L)
  expect_match(diagnostic$Detail, "elapsed", fixed = TRUE)
  expect_match(diagnostic$Detail, "time_types", fixed = TRUE)

  # Declaring the column removes the diagnostic; a clock is modelled with
  # wrap-around handling rather than as raw seconds.
  for (declared in c("duration", "clock")) {
    typed <- Run(time_types = stats::setNames(declared, "elapsed"))
    expect_false(any(typed$Diagnostics$Issue == "backend_failure"))
    expect_false(any(typed$Diagnostics$Issue == "ambiguous_temporal_seconds"))
  }

  # The displayed value is still the original time, not the seconds the model saw.
  queue <- suppressMessages(CreateAnomalyReviewQueue(result))
  observed <- queue$Evidence$Observed[queue$Evidence$Variable == "elapsed"]
  observed <- observed[!is.na(observed) & nzchar(observed)]
  if (length(observed)) expect_match(observed, "^[0-9]+:[0-9]{2}:[0-9]{2}$", all = TRUE)
})

test_that("temporal columns reach a model as seconds whatever units they were stored in", {
  expect_identical(ScidrTimeModelNumeric(hms::as_hms(90)), 90)
  expect_identical(ScidrTimeModelNumeric(as.difftime(c(1, 2), units = "mins")), c(60, 120))
  expect_identical(ScidrTimeModelNumeric(as.difftime(2, units = "hours")), 7200)
  # Anything else is handed back untouched, class included.
  expect_identical(ScidrTimeModelNumeric(c(1.5, 2.5)), c(1.5, 2.5))
  expect_identical(ScidrTimeModelNumeric(letters[1:3]), letters[1:3])
})
