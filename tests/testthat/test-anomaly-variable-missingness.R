test_that("missingness compares blank fields without counting the reviewed participant", {
  df_Data <- data.frame(ID = c("same", "same", paste0("p", 3:62)),
    Target = c(NA, 2, rep(1, 30), rep(NA, 30)), Group = rep(c("A", "B"), each = 32)[1:62])
  # A has one blank and 31 observed; B is a legitimate skipped group.
  df_Data$Group <- c(rep("A", 32), rep("B", 30))
  original <- df_Data
  set.seed(102); rng <- .Random.seed
  result <- DetectDataAnomalies(df_Data, id_var = "ID", methods = "missingness")
  target <- subset(result$CellResults, Variable == "Target")
  first <- target[target$.Row == 1, ]
  expect_true(first$Flag)
  expect_equal(first$GroupStatistics[[1]]$n_obs, 30)
  expect_equal(first$GroupStatistics[[1]]$n_observed, 30)
  expect_equal(first$Conditions[[1]][[1]]$column, "Group")
  expect_false(any(target$Flag[target$.Row > 32], na.rm = TRUE))
  expect_identical(df_Data, original)
  expect_identical(.Random.seed, rng)
  expect_false("missingness" %in% result$Metadata$FlagPercentileMethods)
  queue <- CreateAnomalyReviewQueue(result)
  expect_equal(queue$Issues$SourceRow, 1L)
  expect_match(queue$Issues$SuggestedCheck, "skip logic")
  expect_equal(queue$Issues$ReviewCategory, "Unexpected missingness")
  # Tier is now the severity of the finding type: missingness is 6 of 8.
  expect_equal(ScidrReviewNarratives(queue)$Issues$ReviewTier, 6L)
  expect_match(ScidrReviewNarratives(queue)$Issues$ComparisonDetails, "Group.*A")
  expect_match(ScidrReviewNarratives(queue)$Issues$ComparisonDetails, "exclude this participant")
})

test_that("a blank nothing explains is not a finding", {
  # One blank in an almost-complete column, belonging to no group of variables
  # that go missing together and predicted by nothing: there is no comparison
  # group, so the only honest output is no finding at all.
  df_Data <- data.frame(x = c(NA, 1:9), all_missing = NA_real_, complete = 1:10)
  result <- DetectDataAnomalies(df_Data, methods = "missingness")
  expect_true(all(result$Results$Status == "not applicable"))
  expect_equal(nrow(result$CellResults), 0L)
  expect_true("missingness_unexplained" %in% result$Diagnostics$Issue)
  expect_equal(nrow(CreateAnomalyReviewQueue(result)$Issues), 0L)
  expect_error(DetectDataAnomalies(df_Data, methods = "missingness", method_args = list(missingness = list(unknown = 1))), "Unknown")
  expect_error(DetectDataAnomalies(df_Data, methods = "missingness", method_args = list(missingness = list(min_reference = 1))), "integer")
  expect_error(DetectDataAnomalies(df_Data, methods = "missingness", thresholds = c(missingness = 2)), "between")
})

test_that("missingness respects strata, exact raw missingness, types and high cardinality", {
  # A three-variable form skipped by rows 2:21 of the large site; row 1 carries
  # two of the three, which is the finding.
  df_Data <- data.frame(Site = c(rep("large", 50), rep("small", 10), NA),
    x = c(NA, rep(0, 60)), y = c(0, rep(NA, 20), rep(0, 40)), z = c(0, rep(NA, 20), rep(0, 40)),
    category = paste0("unique", 1:61),
    When = as.Date("2020-01-01") + 1:61, elapsed = as.difftime(1:61, units = "mins"))
  df_Data$x[2:21] <- NA
  result <- DetectDataAnomalies(df_Data, stratify_by = "Site", methods = "missingness")
  expect_true(result$Results$Flag[1])
  expect_true("missingness_high_cardinality" %in% result$Diagnostics$Issue)
  expect_true("missingness_blocks_found" %in% result$Diagnostics$Issue)
  first <- subset(result$CellResults, .Row == 1 & Variable == "x")
  expect_equal(first$GroupStatistics[[1]]$block_size, 3L)
  expect_equal(first$GroupStatistics[[1]]$observed_in_block, 2L)
  expect_match(first$Explanation, "2 of the 3 variables in a group that go missing together")
  # The twenty records missing the whole form are not findings.
  expect_false(any(subset(result$CellResults, .Row %in% 2:21)$Flag %in% TRUE))
  expect_identical(result$ReviewData$elapsed, df_Data$elapsed)
})

test_that("a group with too few carriers leaves the blank unavailable, not clean", {
  df_Data <- data.frame(f1 = c(0, rep(NA, 25), rep(0, 14)), f2 = c(NA, rep(NA, 25), rep(0, 14)),
    f3 = c(0, rep(NA, 25), rep(0, 14)))
  result <- DetectDataAnomalies(df_Data, methods = "missingness")
  first <- subset(result$CellResults, .Row == 1)
  expect_true(all(is.na(first$Flag)))
  expect_equal(unique(first$Status), "unavailable")
  expect_match(first$Explanation[1], "too few comparable records")
})

test_that("variable summaries count source rows once and stay outside queue identity", {
  df_Data <- data.frame(a = c(NA, rep(1, 49)), a2 = c(1, rep(NA, 20), rep(1, 29)),
    a3 = c(1, rep(NA, 20), rep(1, 29)), b = rep(c("yes", "no"), each = 25))
  df_Data$a[2:21] <- NA
  result <- DetectDataAnomalies(df_Data, methods = "missingness", codebook = data.frame(Variable = "a", Label = "A <measure>", Category = "Domain"))
  queue <- CreateAnomalyReviewQueue(list(One = result, Two = result))
  saved <- queue
  display <- ScidrReviewNarratives(queue)
  summary <- ScidrReviewVariableSummary(queue, display)
  summary$Table <- summary$Table[summary$Table$Variable == "a", , drop = FALSE]
  expect_equal(nrow(summary$Table), 2L)
  expect_equal(summary$Table$AffectedRecords, c(1L, 1L))
  expect_equal(summary$Table$Variable, c("a", "a"))
  # The group names itself, which is the whole point of reporting the pattern.
  expect_match(summary$Table$StrongestFinding[1], "2 of the 3 variables in a group that go missing together")
  expect_equal(summary$Previews[[1]]$kind, "Missingness")
  expect_identical(queue, saved)
  html <- as.character(ExploreDataAnomalies(queue, IncludePlots = FALSE))
  expect_match(html, "Variables to review")
  expect_match(html, "Anomaly type")
  expect_false(grepl("Record combinations to review", html, fixed = TRUE))
  expect_match(html, "scidrVariable")
  expect_match(html, "Show all findings")
  expect_match(html, "A &lt;measure&gt;")
})

test_that("previews use numeric, categorical, date and circular representations", {
  df_Data <- data.frame(x = c(1:49, 1000), category = c(rep("common", 49), "rare"),
    When = as.Date("2020-01-01") + 1:50, clock = as.difftime(c(1:49, 86399), units = "secs"))
  result <- DetectDataAnomalies(df_Data, methods = "dates", time_types = c(clock = "clock"))
  expect_equal(ScidrReviewMiniDistribution(result, "x", 1L, 50L)$kind, "Distribution")
  expect_equal(ScidrReviewMiniDistribution(result, "category", 1L, 50L)$kind, "Categories")
  expect_true(tail(ScidrReviewMiniDistribution(result, "x", 1L, 50L)$highlighted, 1))
  expect_equal(ScidrReviewMiniDistribution(result, "clock", 1L, 50L)$kind, "Clock")
  expect_equal(sum(ScidrReviewMiniDistribution(result, "When", 1L, 50L)$counts), 50)
})

test_that("Excel uses short wording and preserves immutable evidence", {
  skip_if_not_installed("openxlsx")
  df_Data <- data.frame(x = c(NA, rep(1, 49)), group = rep("A", 50))
  queue <- CreateAnomalyReviewQueue(DetectDataAnomalies(df_Data, methods = "missingness"))
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file))
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  imported <- ImportAnomalyReviewDecisions(queue, file)
  expect_identical(imported$Issues, queue$Issues)
  expect_identical(imported$Metadata$QueueFingerprint, queue$Metadata$QueueFingerprint)
})

test_that("miniature SVG JavaScript is safe and executable", {
  skip_if(Sys.which("node") == "", "Node unavailable")
  script <- tempfile(fileext = ".js")
  on.exit(unlink(script))
  writeLines(c(paste0("const render=", ScidrReviewMiniJS(), ";"),
    "const x=render('a',{a:{kind:'Categories',labels:['<script>'],counts:[5],highlighted:[true]}});",
    "if(!x.includes('<svg') || !x.includes('&lt;script&gt;') || x.includes('<script>')) throw Error('Unsafe SVG');",
    "console.log('PASS');"), script)
  output <- system2(Sys.which("node"), shQuote(script), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  expect_true("PASS" %in% output)
})

test_that("previous long workbook wording remains accepted", {
  skip_if_not_installed("openxlsx")
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file))
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  previous <- ScidrReviewNarratives(queue, concise = FALSE)
  for (sheet in c("Review", "Evidence")) {
    old <- previous[[if (sheet == "Review") "Issues" else "Evidence"]]
    EditReviewWorkbook(file, function(df_Sheet) {
      df_Sheet$`Why flagged` <- old$WhyFlagged
      df_Sheet$`Comparison details` <- old$ComparisonDetails
      df_Sheet
    }, sheet)
  }
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
})

test_that("variable order uses strongest evidence and distinct records, excluding whole-record leads", {
  queue <- CreateAnomalyReviewQueue(ReviewQueueFixture())
  queue$Sources$Analysis$ReviewData$Value <- seq_len(nrow(queue$Sources$Analysis$ReviewData))
  queue$Sources$Analysis$Metadata$Labels["Value"] <- "Example measure"
  display <- ScidrReviewNarratives(queue)
  e <- display$Evidence[rep(1, 6), ]
  e$IssueID <- c("value", "value", "date1", "date2", "date3", "whole")
  e$EvidenceID <- letters[1:6]
  e$Method <- "dates"
  e$SourceRow <- c(1L, 1L, 1L, 2L, 3L, 1L)
  e$OtherRow <- NA_integer_
  e$Variable <- c("Value", "Value", "Visit", "Visit", "Visit", NA)
  e$OtherVariable <- NA_character_
  e$ReviewTier <- c(1L, 1L, 1L, 1L, 1L, 3L)
  e$ScorePercentile <- c(.99, .99, .9, .8, .7, 1)
  display$Evidence <- e
  summary <- ScidrReviewVariableSummary(queue, display)
  expect_equal(summary$Table$Variable, c("Value", "Visit"))
  expect_equal(summary$Table$AffectedRecords, c(1L, 3L))
  expect_false("whole" %in% unlist(lapply(summary$Table$IssueIDs, jsonlite::fromJSON)))
  expect_identical(summary, ScidrReviewVariableSummary(queue, display))
})

test_that("participant summaries do not double-count missing cells", {
  # Five variables on one form: row 1 carries three of them, so two blanks are
  # findings; rows 2:21 skipped the form entirely and are not.
  df_Data <- data.frame(ID = sprintf("%03d", 1:50), x = c(NA, rep(1, 49)), y = c(NA, rep(2, 49)),
    z = c(1, rep(2, 49)), v = c(1, rep(2, 49)), w = c(1, rep(2, 49)))
  for (column in c("x", "y", "z", "v", "w")) df_Data[[column]][2:21] <- NA
  result <- DetectDataAnomalies(df_Data, id_var = "ID", methods = "missingness")
  expect_equal(nrow(subset(result$CellResults, .Row == 1 & Flag %in% TRUE)), 2L)
  participant <- subset(result$ParticipantSummary, ID == "001")
  expect_equal(participant$FlaggedRecords, 1L)
  expect_equal(participant$MethodsFlagged, 1L)
  expect_equal(result$Summary$Flagged, 1L)
})

test_that("missingness is judged against who the record is, not overall completeness", {
  # Menopause status is collected from women only, so half the column is blank by
  # design. Only the women with no value are findings.
  df_Data <- data.frame(Sex = rep(c("F", "M"), each = 60),
    Age = rep(c(41:70, 41:70), times = 2),
    Meno = c(rep("post", 57), NA, NA, NA, rep(NA_character_, 60)), stringsAsFactors = FALSE)
  result <- DetectDataAnomalies(df_Data, methods = "missingness")
  cells <- subset(result$CellResults, Variable == "Meno")
  expect_setequal(cells$.Row[which(cells$Flag)], 58:60)
  expect_false(any(cells$Flag[cells$.Row > 60], na.rm = TRUE))
  flagged <- cells[cells$.Row == 58, ]
  expect_equal(flagged$Conditions[[1]][[1]]$column, "Sex")
  expect_match(flagged$Explanation, "records with Sex = F have a value")
})

test_that("missingness finds conditions a one-field comparison cannot see", {
  # Expected blank when Male OR under 18; neither field alone explains the pattern.
  set.seed(7)
  df_Data <- data.frame(Sex = rep(c("F", "M"), each = 80), Age = rep(c(10:49, 10:49), times = 2),
    stringsAsFactors = FALSE)
  exempt <- df_Data$Sex == "M" | df_Data$Age < 18
  df_Data$Measure <- ifelse(exempt, NA_real_, 1)
  # Well inside the predictable region: at the 18 boundary minbucket cannot resolve
  # the split finely enough, and the comparator correctly declines to judge.
  violator <- which(!exempt & df_Data$Age == 30)[1]
  df_Data$Measure[violator] <- NA_real_
  result <- DetectDataAnomalies(df_Data, methods = "missingness")
  cells <- subset(result$CellResults, Variable == "Measure")
  expect_true(isTRUE(cells$Flag[cells$.Row == violator]))
  expect_false(any(cells$Flag[cells$.Row %in% which(exempt)], na.rm = TRUE))
  columns <- vapply(cells$Conditions[cells$.Row == violator][[1]], function(x) x$column, character(1))
  expect_true(all(c("Sex", "Age") %in% columns))
})

test_that("missingness that nothing predicts is reported as unpredictable, not flagged", {
  set.seed(11)
  df_Data <- data.frame(Group = rep(c("a", "b"), each = 100), Other = rnorm(200))
  df_Data$Measure <- ifelse(runif(200) < 0.3, NA_real_, 1)
  result <- DetectDataAnomalies(df_Data, methods = "missingness")
  cells <- subset(result$CellResults, Variable == "Measure")
  expect_true(nrow(cells) > 0)
  expect_true(all(cells$Status %in% c("evaluated", "unavailable")))
  expect_false(any(cells$Flag, na.rm = TRUE))
})

test_that("every finding carries one anomaly type from the canonical list", {
  df_Types <- ScidrAnomalyTypes()
  expect_equal(anyDuplicated(df_Types$Key), 0L)
  expect_equal(anyDuplicated(df_Types$Type), 0L)
  combinations <- expand.grid(method = ScidrAnomalyMethods()$Key,
    kind = c("cell", "record", "ranked", "interval", "date_column", "date_interval"),
    field = c(TRUE, FALSE), conditioned = c(TRUE, FALSE), stringsAsFactors = FALSE)
  types <- mapply(ScidrAnomalyType, combinations$method, combinations$kind,
    combinations$field, combinations$conditioned)
  expect_true(all(types %in% df_Types$Type))
  expect_false(any(is.na(types)))
  expect_equal(ScidrAnomalyType("conditional", "cell", TRUE, TRUE, coded = TRUE), "Rare coded value")
})

test_that("the completion message reports one elapsed time, not one per record", {
  set.seed(7)
  n <- 120
  df_Data <- data.frame(ID = sprintf("P%03d", seq_len(n)), Age = round(stats::rnorm(n, 50, 10)),
    wbc = stats::rnorm(n, 7), rbc = stats::rnorm(n, 5), hgb = stats::rnorm(n, 14), hct = stats::rnorm(n, 42))
  df_Data[1:35, c("wbc", "rbc", "hgb", "hct")] <- NA
  df_Data$hgb[100] <- NA
  # The co-missing-group branch assigns a per-record logical vector to a local
  # name that once shadowed the screen's timer, so proc.time() - started gave one
  # elapsed time per record and the message pasted them all together.
  messages <- capture.output(
    DetectDataAnomalies(df_Data, id_var = "ID", methods = "missingness"),
    type = "message")
  completion <- grep("Missingness screening completed", messages, value = TRUE)
  expect_length(completion, 1L)
  expect_match(completion, "^Missingness screening completed in [0-9]+\\.?[0-9]* seconds\\.$")
})
