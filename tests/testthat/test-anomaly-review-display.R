NarrativeFixture <- function(method = "conditional") {
  result <- ReviewQueueFixture()
  queue <- CreateAnomalyReviewQueue(result)
  queue$Sources$Analysis$ReviewData$Value <- seq_len(nrow(result$ReviewData))
  queue$Sources$Analysis$Metadata$Labels["Value"] <- "Example measure"
  queue$Sources$Analysis$Metadata$Labels["Group"] <- "Comparison group"
  queue$Evidence <- queue$Evidence[1, ]
  queue$Issues <- queue$Issues[1, ]
  queue$Evidence$Method <- method
  queue$Evidence$Variable <- "Value"
  queue$Evidence$Observed <- "147"
  queue$Evidence$ReferenceN <- 80
  queue$Evidence$ReferenceLower <- NA_character_
  queue$Evidence$ReferenceUpper <- "60"
  queue$Evidence$GroupStatistics <- list(list(upper_thr = 60, mean = 25, n_obs = 80))
  queue$Evidence$Conditions <- list(list(list(column = "Group", value_this = "A", comparison = "=", value_comp = "A")))
  queue
}

test_that("conditional narratives explain direction and comparison without inventing rules", {
  skip_if_not_installed("digest")
  queue <- NarrativeFixture()
  before <- queue
  display <- ScidrReviewNarratives(queue)
  expect_match(display$Evidence$WhyFlagged, "147.*above.*60")
  expect_match(display$Evidence$ComparisonDetails, "Comparison group \\(Group\\)")
  expect_match(display$Evidence$ComparisonDetails, "Reference group mean: 25")
  expect_match(display$Evidence$ComparisonDetails, "Reference observations: 80")
  expect_identical(queue, before)
  queue$Evidence$Observed <- "-100"
  queue$Evidence$ReferenceLower <- "0"
  display <- ScidrReviewNarratives(queue)
  expect_match(display$Evidence$WhyFlagged, "below")
  queue$Evidence$Conditions <- list(list())
  expect_match(ScidrReviewNarratives(queue)$Evidence$WhyFlagged, "column")
  queue$Evidence$ReferenceLower <- queue$Evidence$ReferenceUpper <- NA_character_
  queue$Evidence$GroupStatistics <- list(list(categs_common = "B", pct_common = .99, n_obs = 80))
  expect_match(ScidrReviewNarratives(queue)$Evidence$ComparisonDetails, "Common response.*B.*99%")
  queue$Evidence$GroupStatistics <- list(list())
  queue$Evidence$Explanation <- "Original backend explanation"
  display <- ScidrReviewNarratives(queue)
  expect_equal(display$Evidence$WhyFlagged, "Original backend explanation")
  expect_match(display$Evidence$ComparisonDetails, "not retained")
})

test_that("date, longitudinal and record narratives describe their native evidence", {
  skip_if_not_installed("digest")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  expect_match(ScidrReviewNarratives(queue)$Issues$WhyFlagged[1], "learned date boundaries")
  queue <- NarrativeFixture("dates")
  queue$Evidence$Kind <- "date_interval"
  queue$Evidence$OtherVariable <- "Visit"
  queue$Evidence$ActualIntervalDays <- -50
  # Signed day counts now read as a direction rather than a negative number.
  expect_match(ScidrReviewNarratives(queue)$Issues$WhyFlagged, "falls 50 days before")
  queue <- NarrativeFixture("change")
  queue$Evidence$OtherObserved <- "100"
  queue$Evidence$Elapsed <- 2
  queue$Evidence$ElapsedUnit <- "days"
  queue$Evidence$Score <- 8
  queue$Evidence$NativeEvidence <- list(list(Rate = 23.5))
  expect_match(ScidrReviewNarratives(queue)$Evidence$ComparisonDetails, "100 to 147 over 2 days.*23.5")
  queue$Evidence$Method <- "transition"
  queue$Evidence$NativeEvidence <- list(list(TransitionCount = 1, ConditioningCount = 50, ReferenceParticipants = 20))
  expect_match(ScidrReviewNarratives(queue)$Evidence$ComparisonDetails, "1 times in 50.*excluding this participant")
  queue$Evidence$Method <- "trajectory"
  queue$Evidence$Expected <- "100"
  expect_match(ScidrReviewNarratives(queue)$Evidence$ComparisonDetails, "in-sample")
  queue$Evidence$Method <- "isolation"
  expect_match(ScidrReviewNarratives(queue)$Evidence$ComparisonDetails, "does not identify a single")
})

test_that("readable Excel projections preserve legacy imports and reject changed display evidence", {
  skip_if_not_installed("openxlsx")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture(), id_var = "ID")
  before <- queue
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  displayed <- openxlsx::read.xlsx(file, "Evidence", check.names = FALSE, sep.names = " ")
  expect_equal(displayed$`Why flagged`, ScidrReviewNarratives(queue)$Evidence$WhyFlagged)
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  EditReviewWorkbook(file, function(d) {d$`Why flagged`[1] <- "altered"; d}, "Evidence")
  expect_error(ImportAnomalyReviewDecisions(queue, file), "evidence.*changed")
  ExportAnomalyReviewQueue(queue, file, overwrite = TRUE, layout = "full")
  for (sheet in c("Review", "Evidence")) EditReviewWorkbook(file,
    function(d) d[setdiff(names(d), c("Why flagged", "Comparison details"))], sheet)
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  expect_identical(queue, before)
})

test_that("interactive plots are lazy, use learned groups and retain complete histories", {
  skip_if_not_installed("DT")
  if (!nzchar(Sys.which("node"))) skip("Node unavailable")
  queue <- NarrativeFixture()
  source <- queue$Sources$Analysis
  source$ReviewData$Group <- rep(c("A", "B"), length.out = nrow(source$ReviewData))
  source$ReviewData$Value[1] <- 147
  queue$Sources$Analysis <- source
  meta <- list(columns = names(queue$Issues), plotSources = ScidrReviewPlotSources(queue),
    plotConditions = stats::setNames(list(queue$Evidence$Conditions[[1]]), queue$Issues$IssueID))
  values <- as.list(queue$Issues[1, ])
  values$Variable <- "Value"
  fn_Input <- tempfile(fileext = ".json"); fn_Script <- tempfile(fileext = ".js")
  on.exit(unlink(c(fn_Input, fn_Script)), add = TRUE)
  jsonlite::write_json(list(renderer = ScidrReviewPlotJS(), meta = meta, values = unname(values)), fn_Input, auto_unbox = TRUE)
  writeLines(c("const x = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));",
    "let result;global.Plotly = {newPlot:(el,data,layout,config)=>{result={data,layout,config};}};",
    "global.document={createElement:()=>({style:{}})};const panel={prepend:()=>{}};",
    "const render=eval('('+x.renderer+')');render(panel,x.values,x.meta);",
    "if(result.data[0].x.length !== 20 || !result.data[0].x.includes(147)) throw Error('Incorrect learned group');",
    "if(!result.config.scrollZoom || !result.config.modeBarButtonsToRemove.includes('toImage')) throw Error('Interactive options');",
    "console.log('PASS');"), fn_Script)
  out <- system2(Sys.which("node"), c(shQuote(fn_Script), shQuote(fn_Input)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n")); expect_true("PASS" %in% out)
  df_Long <- expand.grid(ID = 1:12, Time = 0:3)
  set.seed(51); df_Long$Value <- rnorm(nrow(df_Long)); df_Long$Value[1] <- 20
  long <- DetectLongitudinalAnomalies(df_Long, "Value", "ID", "Time", methods = "change", thresholds = c(change = 3))
  long_queue <- CreateAnomalyReviewQueue(long)
  if (!nzchar(system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly"))) skip("Plotly JavaScript unavailable")
  table <- Filter(function(x) length(x$x$options$scidrMeta$evidence) > 0,
    ExplorerTables(ExploreDataAnomalies(long_queue, IncludePlots = TRUE)))[[1]]
  plot_meta <- table$x$options$scidrMeta
  expect_equal(length(plot_meta$plotSources$Analysis$rows), nrow(df_Long))
  expect_equal(length(plot_meta$plotSources$Analysis$values$Time), nrow(df_Long))
  expect_match(as.character(table$x$callback), "Plotly.newPlot", fixed = TRUE)
})

test_that("conditional plot conditions and learned bounds reach the panel for numeric columns", {
  skip_if_not_installed("outliertree")
  skip_if_not_installed("DT")
  set.seed(12)
  n <- 160
  Statin <- rep(c(0, 1), each = n / 2)
  df_Data <- data.frame(ID = paste0("S", seq_len(n)), Statin = Statin, stringsAsFactors = FALSE)
  for (j in 1:8) df_Data[[paste0("Lab", j)]] <- round(stats::rnorm(n, 100, 20), 1)
  df_Data$Trig <- round(ifelse(Statin == 0, stats::rnorm(n, 120, 10), stats::rnorm(n, 260, 10)))
  df_Data$Trig[5] <- 621
  Result <- suppressMessages(DetectDataAnomalies(df_Data, methods = "conditional"))
  Tables <- ExplorerTables(suppressMessages(ExploreDataAnomalies(Result)))
  Meta <- Filter(function(x) !is.null(x$x$options$scidrMeta$plotConditions), Tables)
  expect_gt(length(Meta), 0L)
  Conditions <- Meta[[1]]$x$options$scidrMeta$plotConditions
  Conditions <- Filter(function(x) length(x) > 0L, Conditions)
  expect_gt(length(Conditions), 0L)
  # The condition names a numeric column while its value is carried as text; the
  # panel casts both sides, so a strict comparison no longer empties the chart.
  expect_true(any(vapply(Conditions, function(conditions)
    any(vapply(conditions, function(condition) identical(condition$column, "Statin"), logical(1))), logical(1))))
  Bounds <- Filter(Negate(is.null), Meta[[1]]$x$options$scidrMeta$plotBounds)
  expect_gt(length(Bounds), 0L)
  expect_true(all(vapply(Bounds, function(b) is.numeric(b$upper) || is.numeric(b$lower), logical(1))))
  expect_match(ScidrReviewPlotJS(), "const cast = function", fixed = TRUE)
})

test_that("common reasons summarise a variable's findings and recorded values are rounded", {
  Reasons <- ScidrReviewCommonReasons(
    c("Unexpected missingness", "Unexpected missingness", "Conditional outlier"),
    c("Blank, although every other record has a value.", "Blank again.", "621.5 is above the boundary."))
  expect_match(Reasons, "Unexpected missingness (2 records)", fixed = TRUE)
  expect_match(Reasons, "Conditional outlier", fixed = TRUE)
  # The most common kind leads, and each kind carries one example.
  expect_lt(regexpr("Unexpected missingness", Reasons, fixed = TRUE),
            regexpr("Conditional outlier", Reasons, fixed = TRUE))
  expect_match(Reasons, "Blank, although every other record has a value.", fixed = TRUE)
  Many <- ScidrReviewCommonReasons(c("A", "B", "C", "D"), c("a", "b", "c", "d"), limit = 2L)
  expect_match(Many, "Plus 2 other kinds of finding.", fixed = TRUE)
  expect_equal(ScidrReviewCommonReasons("Only", NA_character_), "Only")

  expect_equal(ScidrReviewRoundedValue("4980.49871276928"), "4980.5")
  expect_equal(ScidrReviewRoundedValue("621"), "621")
  expect_equal(ScidrReviewRoundedValue("0.63612345"), "0.636")
  # Text, dates and categories are left exactly as they are.
  expect_equal(ScidrReviewRoundedValue("2020-01-01"), "2020-01-01")
  expect_equal(ScidrReviewRoundedValue("Female"), "Female")
})

test_that("the identifier column is headed with the ID column's own label", {
  set.seed(44)
  n <- 80
  df_Data <- data.frame(PID = sprintf("S%03d", seq_len(n)),
    x = stats::rnorm(n), y = stats::rnorm(n))
  df_Data$y[1] <- 25
  attr(df_Data$PID, "label") <- "Participant number"
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "PID", methods = "conditional")
  html <- as.character(ExploreDataAnomalies(result))
  expect_true(any(grepl(">Participant number<", html, fixed = TRUE)))
  expect_true(any(grepl(">Other Participant number<", html, fixed = TRUE)))
  expect_false(any(grepl("Study ID", html, fixed = TRUE)))
})

test_that("an unlabelled identifier column is headed with its column name", {
  set.seed(45)
  n <- 80
  df_Data <- data.frame(PID = sprintf("S%03d", seq_len(n)),
    x = stats::rnorm(n), y = stats::rnorm(n))
  df_Data$y[1] <- 25
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "PID", methods = "conditional")
  html <- as.character(ExploreDataAnomalies(result))
  expect_true(any(grepl(">PID<", html, fixed = TRUE)))
  expect_false(any(grepl("Study ID", html, fixed = TRUE)))
})

test_that("a codebook label for the identifier column reaches the heading", {
  set.seed(46)
  n <- 80
  df_Data <- data.frame(PID = sprintf("S%03d", seq_len(n)),
    x = stats::rnorm(n), y = stats::rnorm(n))
  df_Data$y[1] <- 25
  result <- DetectDataAnomalies(df_Data, c("x", "y"), id_var = "PID", methods = "conditional")
  codebook <- data.frame(Variable = c("PID", "x"), Label = c("Subject code", "Measure x"),
    Category = c("Identifiers", "Labs"))
  html <- as.character(ExploreDataAnomalies(result, codebook = codebook))
  expect_true(any(grepl(">Subject code<", html, fixed = TRUE)))
})

test_that("a detection with no identifier column keeps the generic heading", {
  set.seed(47)
  n <- 80
  df_Data <- data.frame(x = stats::rnorm(n), y = stats::rnorm(n))
  df_Data$y[1] <- 25
  result <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "conditional")
  expect_silent(html <- as.character(ExploreDataAnomalies(result)))
  expect_true(length(html) > 0L)
})
