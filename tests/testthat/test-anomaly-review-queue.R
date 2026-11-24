test_that("queues retain optional IDs and deterministic snapshot identity", {
  skip_if_not_installed("digest")
  result <- ReviewQueueFixture()
  original <- result
  queue <- CreateAnomalyReviewQueue(result)
  expect_s3_class(queue, "AnomalyReviewQueue")
  expect_equal(nrow(queue$Issues), 2)
  expect_true(all(is.na(queue$Issues$StudyID)))
  expect_equal(queue$Issues$SourceRow, 1:2)
  expect_equal(queue$Issues$Priority, rep("Unassigned", 2))
  expect_equal(queue$Issues$ReviewStatus, rep("Unreviewed", 2))
  expect_identical(queue$Issues$IssueID, CreateAnomalyReviewQueue(result)$Issues$IssueID)
  expect_identical(result, original)
  identified <- CreateAnomalyReviewQueue(result, id_var = "ID", context_vars = "Site")
  expect_equal(identified$Issues$StudyID, c("001", "002"))
  expect_match(identified$Issues$Context[1], "Assessment site: A")
  expect_equal(identified$Issues$IssueID, queue$Issues$IssueID)
  expect_equal(CreateAnomalyReviewQueue(ReviewQueueFixture(TRUE))$Issues$StudyID, c("001", "002"))
  expect_error(CreateAnomalyReviewQueue(result, id_var = "Absent"), "retained")
  changed <- result
  changed$ReviewData <- changed$ReviewData[40:1, ]
  expect_false(identical(CreateAnomalyReviewQueue(changed)$Metadata$Snapshots, queue$Metadata$Snapshots))
})

test_that("duplicate or missing study IDs never collapse rows or block review", {
  skip_if_not_installed("digest")
  result <- ReviewQueueFixture()
  result$ReviewData$ID[] <- "same"
  result$ReviewData$ID[2] <- NA
  queue <- CreateAnomalyReviewQueue(result, id_var = "ID")
  expect_equal(nrow(queue$Issues), 2)
  expect_equal(length(unique(queue$Issues$IssueID)), 2)
  expect_true(is.na(queue$Issues$StudyID[2]))
  expect_true(all(c("missing_id", "duplicate_id") %in% queue$Diagnostics$Issue))
})

test_that("analyses retain independent row namespaces", {
  skip_if_not_installed("digest")
  result <- ReviewQueueFixture()
  queue <- CreateAnomalyReviewQueue(list(Clinical = result, Dates = result))
  expect_equal(nrow(queue$Issues), 4)
  expect_equal(length(unique(queue$Issues$IssueID)), 4)
  expect_equal(length(queue$Sources), 2)
  expect_error(CreateAnomalyReviewQueue(list(result, result)), "uniquely named")
  expect_error(CreateAnomalyReviewQueue(setNames(list(result, result), c("same", "same"))), "uniquely named")
})

test_that("unclassified rankings default on, allow opt-out and are limited per stratum", {
  skip_if_not_installed("digest")
  result <- ReviewQueueFixture()
  result$Results$Method <- "isolation"
  result$Results$Flag <- NA
  result$Results$Score <- seq_len(40)
  result$Results$Stratum <- rep(1:2, each = 20)
  result$DateResults <- result$DateResults[FALSE, ]
  queue <- CreateAnomalyReviewQueue(result, include_ranked = FALSE)
  expect_equal(nrow(queue$Issues), 0)
  expect_equal(nrow(CreateAnomalyReviewQueue(result)$Issues), 20)
  ranked <- CreateAnomalyReviewQueue(result, include_ranked = TRUE, top_n = 3)
  expect_equal(nrow(ranked$Issues), 6)
  expect_setequal(ranked$Issues$SourceRow, c(18:20, 38:40))
  expect_true(all(ranked$Issues$ReviewCategory == "Uncommon combination"))
  expect_true(all(is.na(ranked$Issues$Variable)))
  result$Results$Flag[40] <- TRUE
  expect_equal(CreateAnomalyReviewQueue(result, include_ranked = FALSE)$Issues$SourceRow, 40)
})

test_that("matching cell evidence consolidates but distinct date pairs do not", {
  skip_if_not_installed("digest")
  result <- ReviewQueueFixture()
  conditional <- result$Results
  conditional$Method <- "conditional"
  conditional$Flag <- seq_len(40) == 1
  result$Results <- dplyr::bind_rows(result$Results, conditional)
  result$CellResults <- tibble::tibble(.Row = 1L, Method = "conditional", Variable = "Visit",
    Flag = TRUE, Explanation = "Conditional evidence", Conditions = list(list()))
  queue <- CreateAnomalyReviewQueue(result)
  issue <- queue$Issues[queue$Issues$SourceRow == 1, ]
  expect_equal(nrow(issue), 1)
  expect_equal(issue$EvidenceCount, 2)
  expect_equal(nrow(queue$Evidence[queue$Evidence$IssueID == issue$IssueID, ]), 2)
  df_Data <- data.frame(First = as.Date("2020-01-01") + 1:40)
  df_Data$Second <- df_Data$Third <- df_Data$First + 100
  df_Data$Second[1] <- df_Data$Third[1] <- df_Data$First[1] + 1000
  queue <- CreateAnomalyReviewQueue(DetectDataAnomalies(df_Data, names(df_Data), methods = "dates"))
  pairs <- queue$Issues[queue$Issues$ReviewCategory == "Implausible date interval", ]
  expect_equal(nrow(pairs), 2)
  expect_equal(length(unique(pairs$IssueID)), 2)
  expect_setequal(pairs$OtherVariable, c("Second", "Third"))
})

test_that("new structured conditional statistics and older findings are supported", {
  skip_if_not_installed("digest")
  skip_if_not_installed("outliertree")
  set.seed(12)
  df_Data <- data.frame(x = rnorm(150), y = rnorm(150))
  df_Data$y[1] <- 35
  result <- DetectDataAnomalies(df_Data, c("x", "y"), methods = "conditional")
  expect_true(is.list(result$CellResults$GroupStatistics[[1]]))
  queue <- CreateAnomalyReviewQueue(result)
  expect_true(any(!is.na(queue$Evidence$ReferenceN)))
  result$CellResults$GroupStatistics <- NULL
  old <- CreateAnomalyReviewQueue(result)
  expect_true(all(is.na(old$Evidence$ReferenceN)))
  expect_true(any(nzchar(old$Evidence$Explanation)))
})

test_that("longitudinal intervals retain both endpoint rows without duplicate record issues", {
  skip_if_not_installed("digest")
  df_Data <- expand.grid(ID = 1:12, Time = 0:3)
  set.seed(16)
  df_Data$x <- rnorm(nrow(df_Data))
  df_Data$x[1] <- 25
  result <- DetectLongitudinalAnomalies(df_Data, "x", "ID", "Time", methods = "change", thresholds = c(change = 2))
  queue <- CreateAnomalyReviewQueue(result)
  expect_true(nrow(queue$Issues) > 0)
  expect_true(all(!is.na(queue$Issues$OtherRow)))
  expect_true(all(queue$Evidence$Kind == "interval"))
  expect_true(all(queue$Issues$ReviewCategory == "Unexpected change over time"))
  expect_true(all(!is.na(queue$Evidence$Elapsed)))
})

test_that("Excel round trips preserve leading-zero IDs, blanks and partial decisions", {
  skip_if_not_installed("digest")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("jsonlite")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture(), id_var = "ID")
  original <- queue
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  expect_error(ExportAnomalyReviewQueue(queue, file, layout = "full"), "already exists")
  expect_identical(queue, original)
  wb <- openxlsx::loadWorkbook(file)
  expect_equal(openxlsx::sheetVisibility(wb)[4], "hidden")
  d <- openxlsx::read.xlsx(wb, "Review", check.names = FALSE, sep.names = " ")
  expect_equal(d$`Study ID`, c("001", "002"))
  EditReviewWorkbook(file, function(d) {
    d$`Review status`[1] <- "Correction needed"
    d$Priority[1] <- "High"
    d$`Proposed correction`[1] <- "Visit = 2020-01-02"
    d$`Reviewed on`[1] <- "2026-11-19"
    d$Notes[1] <- "Checked source"
    d[2:1, ]
  })
  updated <- ImportAnomalyReviewDecisions(queue, file)
  expect_equal(updated$Issues$ReviewStatus[1], "Correction needed")
  expect_equal(updated$Issues$Priority[1], "High")
  expect_equal(updated$Issues$ProposedCorrection[1], "Visit = 2020-01-02")
  expect_identical(updated$Sources, original$Sources)
  ExportAnomalyReviewQueue(updated, file, overwrite = TRUE, layout = "full")
  EditReviewWorkbook(file, function(d) { d$`Review status`[1] <- ""; d$Notes[1] <- ""; d[1, ] })
  unchanged <- ImportAnomalyReviewDecisions(updated, file)
  expect_equal(unchanged$Issues$ReviewStatus, updated$Issues$ReviewStatus)
  expect_equal(unchanged$Issues$Notes, updated$Issues$Notes)
  ExportAnomalyReviewQueue(updated, file, overwrite = TRUE, layout = "full")
  EditReviewWorkbook(file, function(d) { d$`Review status`[1] <- "Unreviewed"; d[1, ] })
  expect_equal(ImportAnomalyReviewDecisions(updated, file)$Issues$ReviewStatus[1], "Unreviewed")
})

test_that("imports reject changed evidence, IDs, decisions and snapshots", {
  skip_if_not_installed("digest")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("jsonlite")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  Mutate <- function(transform, pattern, sheet = "Review") {
    ExportAnomalyReviewQueue(queue, file, overwrite = TRUE, layout = "full")
    EditReviewWorkbook(file, transform, sheet)
    expect_error(ImportAnomalyReviewDecisions(queue, file), pattern)
  }
  Mutate(function(d) { d$`Recorded value`[1] <- "Changed"; d }, "evidence changed")
  Mutate(function(d) rbind(d, d[1, ]), "duplicate")
  Mutate(function(d) { d$`Issue ID`[1] <- "unknown"; d }, "unknown")
  Mutate(function(d) { d$`Review status`[1] <- "Fixed"; d }, "Invalid ReviewStatus")
  Mutate(function(d) { d$Priority[1] <- "Urgent"; d }, "Invalid Priority")
  Mutate(function(d) { d$`Reviewed on`[1] <- "2026-02-30"; d }, "YYYY-MM-DD")
  Mutate(function(d) { d$Score[1] <- "999"; d }, "Evidence original", "Evidence")
  Mutate(function(d) { d$Snapshot[1] <- "different"; d }, "Manifest original", "Manifest")
  ExportAnomalyReviewQueue(queue, file, overwrite = TRUE, layout = "full")
  altered <- queue
  altered$Sources[[1]]$ReviewData$Visit[1] <- as.Date("2020-01-01")
  expect_error(ImportAnomalyReviewDecisions(altered, file), "Source snapshot changed")
  altered <- queue
  altered$Issues$Observed[1] <- "changed"
  expect_error(ExportAnomalyReviewQueue(altered, file, overwrite = TRUE, layout = "full"), "Queue evidence changed")
  other <- CreateAnomalyReviewQueue(ReviewQueueFixture(), include_ranked = TRUE)
  # Same candidates and snapshots are compatible; changed source data are not.
  result <- ReviewQueueFixture()
  result$ReviewData$Visit[1] <- as.Date("2097-01-01")
  other <- CreateAnomalyReviewQueue(result)
  expect_error(ImportAnomalyReviewDecisions(other, file), "Manifest original")
})

test_that("empty queues, midnight timestamps and workbook validation survive export", {
  skip_if_not_installed("digest")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("jsonlite")
  df_Data <- data.frame(When = as.POSIXct("2020-01-01", tz = "America/New_York") + seq_len(40) * 86400)
  result <- DetectDataAnomalies(df_Data, "When", methods = "dates")
  queue <- CreateAnomalyReviewQueue(result)
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  expect_equal(nrow(queue$Issues), 0)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  expect_equal(ScidrReviewFormat(as.POSIXct("2020-01-01", tz = "America/New_York")), "2020-01-01")
  expect_identical(result$ReviewData$When, df_Data$When)
  expect_error(CreateAnomalyReviewQueue(result, top_n = 0), "positive")
})

test_that("explorer accepts queues, multiple analyses and read-only decisions", {
  skip_if_not_installed("digest")
  skip_if_not_installed("DT")
  queue <- CreateAnomalyReviewQueue(list(Dates = ReviewQueueFixture(), Other = ReviewQueueFixture()), id_var = "ID")
  queue$Issues$Notes[1] <- "<script>alert('x')</script>"
  queue$Issues$ReviewStatus[1] <- "Confirmed correct"
  original <- queue
  explorer <- ExploreDataAnomalies(queue)
  html <- htmltools::renderTags(explorer)$html
  expect_match(html, "Potential data issues")
  expect_match(html, "at most one finding per record")
  expect_match(html, "Coverage")
  expect_identical(queue, original)
})

test_that("workbook validation and freeze settings are emitted and evidence can be sorted", {
  skip_if_not_installed("digest")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("jsonlite")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture(), id_var = "ID")
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  con <- unz(file, "xl/worksheets/sheet1.xml")
  xml <- paste(readLines(con, warn = FALSE), collapse = "")
  close(con)
  expect_match(xml, 'dataValidations count="2"', fixed = TRUE)
  expect_match(xml, 'state="frozen"', fixed = TRUE)
  EditReviewWorkbook(file, function(d) d[nrow(d):1, ], "Evidence")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
})

test_that("queue row expansion exposes safe structured evidence and hides tracking noise", {
  skip_if_not_installed("digest")
  skip_if_not_installed("DT")
  skip_if_not_installed("jsonlite")
  if (!nzchar(Sys.which("node"))) skip("Node is unavailable for the queue renderer test")
  result <- ReviewQueueFixture()
  queue <- CreateAnomalyReviewQueue(result)
  queue$Issues$Notes[1] <- "<script>unsafe</script>"
  Collect <- function(node) {
    if (inherits(node, "datatables")) return(list(node))
    if (inherits(node, "shiny.tag")) return(unlist(lapply(node$children, Collect), recursive = FALSE))
    if (is.list(node)) return(unlist(lapply(node, Collect), recursive = FALSE))
    list()
  }
  tables <- Collect(ExploreDataAnomalies(queue))
  worklist <- Filter(function(x) length(x$x$options$scidrMeta$evidence) > 0, tables)[[1]]
  meta <- worklist$x$options$scidrMeta
  values <- lapply(worklist$x$data, function(x) x[1])
  input <- list(renderer = ScidrAnomalyDetailJS(), meta = meta, values = unname(values))
  fn_Input <- tempfile(fileext = ".json")
  fn_Script <- tempfile(fileext = ".js")
  on.exit(unlink(c(fn_Input, fn_Script)), add = TRUE)
  jsonlite::write_json(input, fn_Input, auto_unbox = TRUE)
  writeLines(c("const x = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));",
    "const render = eval('(' + x.renderer + ')');",
    "const html = render(x.values, x.meta);",
    "if (!html.includes('Supporting evidence') || !html.includes('Reference observations')) throw Error('Missing structured evidence');",
    "if (html.includes('<script>') || !html.includes('&lt;script&gt;')) throw Error('Unescaped note');",
    "if (!html.includes('<summary>Tracking details</summary>')) throw Error('Tracking fields not separated');",
    "console.log('PASS');"), fn_Script)
  out <- system2(Sys.which("node"), c(shQuote(fn_Script), shQuote(fn_Input)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(out, "status"))
  expect_true("PASS" %in% out)
})

test_that("the compact workbook is one sheet of scannable columns", {
  skip_if_not_installed("openxlsx")
  set.seed(21)
  df_Data <- data.frame(ID = sprintf("P%02d", 1:60),
    bmi = c(stats::rnorm(59, 27, 4), 9000), chol = stats::rnorm(60, 180, 30))
  queue <- suppressMessages(CreateAnomalyReviewQueue(
    suppressMessages(DetectDataAnomalies(df_Data, c("bmi", "chol"), id_var = "ID"))))
  compact <- withr::local_tempfile(fileext = ".xlsx")
  full <- withr::local_tempfile(fileext = ".xlsx")
  ExportAnomalyReviewQueue(queue, compact, overwrite = TRUE)
  ExportAnomalyReviewQueue(queue, full, overwrite = TRUE, layout = "full")

  workbook <- openxlsx::loadWorkbook(compact)
  columns <- names(openxlsx::read.xlsx(compact, "Review", check.names = FALSE, sep.names = " "))
  expect_equal(names(workbook), "Review")
  # No codebook was supplied, so the category column cannot apply; every other
  # compact column is present, in order.
  expect_equal(columns, setdiff(ScidrReviewCompactColumns(), ScidrReviewHeaders("CodebookCategory")))
  # Every row is one line high, whatever its text.
  expect_equal(unique(as.numeric(unlist(workbook$rowHeights[[1]])[-1])), 22)
  expect_true(file.size(compact) < file.size(full))

  # The long explanation columns are gone, not merely narrowed.
  for (dropped in c("Comparison details", "Order rationale", "Why review", "Columns supplying categories"))
    expect_false(dropped %in% columns)
  # Compact is a deliverable, not a round trip, and says so.
  expect_error(ImportAnomalyReviewDecisions(queue, compact), "compact workbook")
  expect_equal(nrow(ImportAnomalyReviewDecisions(queue, full)$Issues), nrow(queue$Issues))
  expect_error(ExportAnomalyReviewQueue(queue, compact, overwrite = TRUE, layout = "tiny"), "should be one of")
})

test_that("compact drops a column only when it cannot apply", {
  skip_if_not_installed("openxlsx")
  set.seed(22)
  df_Data <- data.frame(bmi = c(stats::rnorm(59, 27, 4), 44), chol = stats::rnorm(60, 180, 30))
  queue <- suppressMessages(CreateAnomalyReviewQueue(
    suppressMessages(DetectDataAnomalies(df_Data, c("bmi", "chol")))))
  file <- withr::local_tempfile(fileext = ".xlsx")
  ExportAnomalyReviewQueue(queue, file, overwrite = TRUE)
  columns <- names(openxlsx::read.xlsx(file, "Review", check.names = FALSE, sep.names = " "))
  # No ID variable was supplied, so Study ID cannot apply and is dropped.
  expect_false("Study ID" %in% columns)
  # The decision fields are empty by design and must survive.
  expect_true(all(ScidrReviewHeaders(ScidrReviewFields()) %in% columns))
  expect_true(all(c("Measure", "Recorded value", "Why flagged") %in% columns))
})

test_that("compact carries the codebook category only when a codebook assigned one", {
  skip_if_not_installed("openxlsx")
  set.seed(23)
  df_Data <- data.frame(ID = sprintf("P%02d", 1:60),
    bmi = c(stats::rnorm(59, 27, 4), 9000), chol = stats::rnorm(60, 180, 30))
  codebook <- data.frame(Variable = c("bmi", "chol"), Label = c("Body mass index", "Cholesterol"),
    Category = c("Anthropometry", "Lipids"))
  category <- ScidrReviewHeaders("CodebookCategory")
  Columns <- function(queue) {
    file <- withr::local_tempfile(fileext = ".xlsx", .local_envir = parent.frame())
    ExportAnomalyReviewQueue(queue, file, overwrite = TRUE)
    list(Names = names(openxlsx::read.xlsx(file, "Review", check.names = FALSE, sep.names = " ")),
      Data = openxlsx::read.xlsx(file, "Review", check.names = FALSE, sep.names = " "))
  }
  # A codebook given once at detection reaches the workbook with no second pass.
  detected <- suppressMessages(DetectDataAnomalies(df_Data, c("bmi", "chol"), id_var = "ID", codebook = codebook))
  from_detection <- Columns(suppressMessages(CreateAnomalyReviewQueue(detected)))
  expect_true(category %in% from_detection$Names)
  expect_true(any(grepl("Anthropometry", from_detection$Data[[category]], fixed = TRUE)))
  # It sits before the extreme value, not at the end among the decision fields.
  expect_lt(match(category, from_detection$Names), match(ScidrReviewHeaders("WhyFlagged"), from_detection$Names))
  expect_false(category %in% names(from_detection$Data)[match(ScidrReviewHeaders(ScidrReviewFields()), from_detection$Names)])

  # A late override on the queue reaches it too.
  plain <- suppressMessages(DetectDataAnomalies(df_Data, c("bmi", "chol"), id_var = "ID"))
  expect_false(category %in% Columns(suppressMessages(CreateAnomalyReviewQueue(plain)))$Names)
  late <- Columns(suppressMessages(CreateAnomalyReviewQueue(plain, codebook = codebook)))
  expect_true(category %in% late$Names)
  expect_true(any(grepl("Lipids", late$Data[[category]], fixed = TRUE)))

  # A codebook naming only unscreened columns assigns nothing, so the column
  # stays out rather than arriving full of "Uncategorized".
  unrelated <- data.frame(Variable = "weight", Label = "Weight", Category = "Anthropometry")
  expect_false(category %in% Columns(suppressMessages(CreateAnomalyReviewQueue(plain, codebook = unrelated)))$Names)
})
