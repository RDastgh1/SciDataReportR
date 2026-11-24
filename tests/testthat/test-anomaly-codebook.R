CategoryFixture <- function(codebook = NULL, Relabel = TRUE) {
  df_Data <- data.frame(ID = sprintf("%03d", 1:40), DateA = as.Date("2020-01-01") + 1:40,
    DateB = as.Date("2020-01-01") + 11:50, value = 1:40)
  df_Data$DateA[1] <- as.Date("2099-01-01")
  attr(df_Data$value, "label") <- "Attached score"
  DetectDataAnomalies(df_Data, id_var = "ID", methods = "dates", Relabel = Relabel, codebook = codebook)
}

CategoryCodebook <- function() data.frame(Variable = c("DateA", "DateB", "value", "unused"),
  Label = c("Repeated label", "Repeated label", NA, "Unused label"),
  Category = c("Dates", "Visits", "Scores", "Other"))

test_that("codebooks supply display metadata without changing detection or source values", {
  codebook <- CategoryCodebook()
  before <- codebook
  rng <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
  original <- CategoryFixture()
  labeled <- CategoryFixture(codebook)
  if (!is.null(rng)) expect_identical(.Random.seed, rng)
  expect_identical(codebook, before)
  for (component in c("Results", "ReviewData", "CellResults", "DateResults", "IntervalResults"))
    expect_identical(original[[component]], labeled[[component]])
  expect_identical(ScidrReviewSnapshot(original), ScidrReviewSnapshot(labeled))
  queue <- CreateAnomalyReviewQueue(labeled)
  unchanged <- CreateAnomalyReviewQueue(original)
  expect_identical(queue$Issues, unchanged$Issues)
  expect_identical(queue$Metadata$QueueFingerprint, unchanged$Metadata$QueueFingerprint)
  metadata <- ScidrReviewCatalog(queue, "Analysis")
  expect_equal(metadata$Label[metadata$Variable == "DateA"], "Repeated label")
  expect_equal(metadata$Label[metadata$Variable == "value"], "Attached score")
  expect_true(all(c("codebook_unmatched_variables", "codebook_unused_entries") %in% labeled$Diagnostics$Issue))
  expect_false(any(labeled$Diagnostics$Issue == "missing_label" & labeled$Diagnostics$Detail == "DateA"))
  raw <- CreateAnomalyReviewQueue(CategoryFixture(codebook, Relabel = FALSE))
  expect_equal(ScidrReviewCatalog(raw, "Analysis")$Label, ScidrReviewCatalog(raw, "Analysis")$Variable)
  expect_equal(with(ScidrReviewCatalog(raw, "Analysis"), Category[Variable == "DateA"]), "Dates")
})

test_that("codebook validation handles duplicates, blanks, missing metadata and exact names", {
  codebook <- CategoryCodebook()
  expect_equal(CategoryFixture(rbind(codebook, codebook[1, ]))$Metadata$VariableMetadata,
    CategoryFixture(codebook)$Metadata$VariableMetadata)
  conflicting <- rbind(codebook, transform(codebook[1, ], Category = "Different"))
  expect_error(CategoryFixture(conflicting), "Conflicting duplicate")
  expect_error(CategoryFixture(data.frame(Label = "x")), "Variable")
  expect_error(CategoryFixture(data.frame(Variable = "")), "nonmissing, nonblank")
  partial <- CategoryFixture(data.frame(Variable = c("DateA", "dateb"), Label = c(NA, NA), Category = c(NA, NA)))
  expect_true(all(partial$Metadata$VariableMetadata$Category == "Uncategorized"))
  expect_true(all(c("codebook_missing_labels", "codebook_missing_categories", "codebook_unused_entries") %in% partial$Diagnostics$Issue))
  expect_equal(CategoryFixture(data.frame(Variable = "DateA"))$Metadata$VariableMetadata$Category[2], "Uncategorized")
  expect_true(all(ScidrReviewNarratives(CreateAnomalyReviewQueue(CategoryFixture()))$Issues$CodebookCategory %in% "Uncategorized"))
})

test_that("date pairs keep both categories and condition categories stay separate", {
  queue <- CreateAnomalyReviewQueue(CategoryFixture(CategoryCodebook()))
  display <- ScidrReviewNarratives(queue)
  pair <- display$Issues[!is.na(display$Issues$OtherVariable), ]
  expect_setequal(jsonlite::fromJSON(pair$CodebookCategoryKeys[1]), c("Dates", "Visits"))
  expect_match(pair$CategoryVariables[1], "Repeated label (DateA)", fixed = TRUE)
  expect_match(pair$CategoryVariables[1], "Repeated label (DateB)", fixed = TRUE)
  e <- queue$Evidence[queue$Evidence$Kind == "date_column", ][1, ]
  e$Method <- "conditional"
  e$Conditions <- list(list(list(column = "value", comparison = ">", value_comp = 4, value_this = 5)))
  decorated <- ScidrReviewCategoryFields(queue, e)
  expect_equal(decorated$CodebookCategory, "Dates")
  expect_match(decorated$ConditionCategories, "Scores")
  expect_false(grepl("Scores", decorated$CodebookCategory, fixed = TRUE))
})

test_that("whole-record categories use five positive effects, including saved SHAP", {
  skip_if_not_installed("isotree")
  set.seed(18)
  df_Data <- as.data.frame(matrix(rnorm(40 * 8), nrow = 40))
  codebook <- data.frame(Variable = names(df_Data), Category = paste0("Category ", seq_len(8)))
  result <- DetectDataAnomalies(df_Data, methods = "isolation", method_args = list(isolation = list(ntrees = 20)),
    explanation_args = list(rows = 1), codebook = codebook)
  unlabelled <- DetectDataAnomalies(df_Data, methods = "isolation", method_args = list(isolation = list(ntrees = 20)),
    explanation_args = list(rows = 1))
  expect_identical(result$Results, unlabelled$Results)
  expect_identical(result$FeatureSensitivity, unlabelled$FeatureSensitivity)
  expect_identical(result$ExplanationMetadata$Fingerprint, unlabelled$ExplanationMetadata$Fingerprint)
  # Set deterministic saved explanation evidence to exercise grouping, not forest ranking.
  result$FeatureSensitivity$ScoreDecrease <- c(8, 7, 6, 5, 4, 3, 0, -1)
  result$ExplanationMetadata$Fingerprint <- ScidrAnomalyExplanationFingerprint(result)
  queue <- CreateAnomalyReviewQueue(result, top_n = 40)
  display <- ScidrReviewNarratives(queue)$Issues
  finding <- display[display$SourceRow == 1, ]
  expect_setequal(jsonlite::fromJSON(finding$CodebookCategoryKeys), paste0("Category ", 1:5))
  expect_equal(finding$CategoryBasis, "Suggested inspection categories")
  expect_false(grepl("V6", finding$CategoryVariables, fixed = TRUE))
  expect_true(any(display$CodebookCategory == "Whole record\u2014category unspecified"))
  old <- result
  old$FeatureContributions <- transform(as.data.frame(old$FeatureSensitivity), Contribution = ScoreDecrease)
  old$FeatureSensitivity <- NULL
  old$ExplanationMetadata$Method <- "shap"
  old$ExplanationMetadata$Fingerprint <- ScidrAnomalyExplanationFingerprint(old)
  expect_equal(ScidrReviewInspectionColumns(old, 1), paste0("V", 1:5))
  no_positive <- result
  no_positive$FeatureSensitivity$ScoreDecrease <- rep(-1, 8)
  expect_length(ScidrReviewInspectionColumns(no_positive, 1), 0)
})

test_that("analysis-specific catalogs and common overrides leave queues unchanged", {
  a <- CategoryFixture(CategoryCodebook())
  other <- CategoryCodebook()
  other$Category[1] <- "Other dates"
  b <- CategoryFixture(other)
  queue <- CreateAnomalyReviewQueue(list(A = a, B = b))
  expect_equal(with(ScidrReviewCatalog(queue, "A"), Category[Variable == "DateA"]), "Dates")
  expect_equal(with(ScidrReviewCatalog(queue, "B"), Category[Variable == "DateA"]), "Other dates")
  override <- data.frame(Variable = "DateB", Label = "Override label", Category = "Override category")
  replaced <- CreateAnomalyReviewQueue(list(A = a, B = b), codebook = override)
  expect_identical(queue$Sources, replaced$Sources)
  expect_identical(queue$Issues, replaced$Issues)
  expect_identical(queue$Metadata$QueueFingerprint, replaced$Metadata$QueueFingerprint)
  for (analysis in c("A", "B")) expect_equal(with(ScidrReviewCatalog(replaced, analysis), Category[Variable == "DateB"]), "Override category")
  before <- queue
  html <- ExploreDataAnomalies(queue, codebook = override, IncludePlots = FALSE)
  expect_identical(queue, before)
  expect_match(htmltools::renderTags(html)$html, "Override category", fixed = TRUE)
})

test_that("longitudinal results retain categories without recoding variables", {
  df_Data <- expand.grid(ID = 1:12, Time = 0:3)
  df_Data$value <- sin(df_Data$ID) + df_Data$Time
  codebook <- data.frame(Variable = "value", Label = "Observed measure", Category = "Measurements", Type = "categorical")
  raw <- DetectLongitudinalAnomalies(df_Data, "value", "ID", "Time", methods = "change")
  result <- DetectLongitudinalAnomalies(df_Data, "value", "ID", "Time", methods = "change", codebook = codebook)
  expect_identical(result$Results, raw$Results)
  expect_identical(result$IntervalResults, raw$IntervalResults)
  expect_true(is.numeric(result$ReviewData$value))
  queue <- CreateAnomalyReviewQueue(result)
  expect_equal(ScidrReviewCatalog(queue, "Analysis")$Category[1], "Measurements")
})

test_that("Excel categories and relabeling accept legacy workbooks and preserve decisions", {
  skip_if_not_installed("openxlsx")
  original <- CreateAnomalyReviewQueue(CategoryFixture())
  queue <- CreateAnomalyReviewQueue(original$Sources$Analysis, codebook = CategoryCodebook())
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  ExportAnomalyReviewQueue(queue, file, layout = "full")
  expect_identical(ImportAnomalyReviewDecisions(queue, file), queue)
  data <- openxlsx::read.xlsx(file, "Review", check.names = FALSE, sep.names = " ")
  expect_true("Codebook category" %in% names(data))
  expect_true(any(grepl("Dates", data$`Codebook category`, fixed = TRUE)))
  ExportAnomalyReviewQueue(original, file, overwrite = TRUE, layout = "full")
  # Emulate a matching pre-category workbook, including its original labels/reasons.
  for (sheet in c("Review", "Evidence", "Contributions")) EditReviewWorkbook(file, function(d)
    d[setdiff(names(d), ScidrReviewHeaders(c("CodebookCategory", "CodebookCategoryKeys", "CategoryBasis", "CategoryVariables", "ConditionCategories")))], sheet)
  edited_id <- openxlsx::read.xlsx(file, "Review", check.names = FALSE, sep.names = " ")$`Issue ID`[1]
  EditReviewWorkbook(file, function(d) { d$`Review status`[1] <- "Correction needed"; d })
  updated <- ImportAnomalyReviewDecisions(queue, file)
  expect_equal(updated$Issues$ReviewStatus[match(edited_id, updated$Issues$IssueID)], "Correction needed")
  expect_identical(updated$Metadata$QueueFingerprint, queue$Metadata$QueueFingerprint)
  expect_identical(updated$Sources, queue$Sources)
  ExportAnomalyReviewQueue(queue, file, overwrite = TRUE, layout = "full")
  EditReviewWorkbook(file, function(d) { d$`Recorded value`[1] <- "Changed"; d }, "Evidence")
  expect_error(ImportAnomalyReviewDecisions(queue, file), "evidence.*changed")
})

test_that("explorer removes decisions from queue columns and explains methods", {
  skip_if_not_installed("DT")
  queue <- CreateAnomalyReviewQueue(CategoryFixture(CategoryCodebook()))
  before <- queue
  html <- ExploreDataAnomalies(queue)
  tables <- ExplorerTables(html)
  issue <- Filter(function(x) "CodebookCategoryKeys" %in% x$x$options$scidrMeta$columns, tables)[[1]]
  columns <- issue$x$options$scidrMeta$columns
  hidden <- issue$x$options$columnDefs[[2]]$targets + 1
  expect_false("ReviewStatus" %in% columns[setdiff(seq_along(columns), hidden)])
  expect_true("CodebookCategory" %in% columns[setdiff(seq_along(columns), hidden)])
  text <- htmltools::renderTags(html)$html
  expect_match(as.character(issue$x$callback), "Anomaly type", fixed = TRUE)
  expect_match(text, "Codebook category", fixed = TRUE)
  expect_match(text, "What the methods mean", fixed = TRUE)
  expect_match(text, "Isolation forest", fixed = TRUE)
  expect_false(grepl("Cell evidence", text, fixed = TRUE))
  expect_identical(queue, before)
})

test_that("revalued inputs and escaped codebook text remain presentation-only", {
  df_Data <- data.frame(value = seq_len(30), response = rep(c(0, 1), 15))
  codebook <- data.frame(Variable = c("value", "response"), Label = c("Observed value", "Observed response"),
    Type = c("continuous", "categorical"), Code = c(NA, "0=No;1=Yes"), Category = c("A", "AA"))
  revalued <- RevalueData(df_Data, codebook)$RevaluedData
  before <- revalued
  result <- DetectDataAnomalies(revalued, methods = "dates", codebook = codebook)
  expect_identical(revalued, before)
  expect_identical(result$ReviewData$response, revalued$response)
  expect_true(is.factor(result$ReviewData$response))
  queue <- CreateAnomalyReviewQueue(CategoryFixture())
  hostile <- data.frame(Variable = "DateA", Label = "</script><img src=x onerror=alert(1)>", Category = "A < B + (C); AA")
  view <- ScidrReviewCodebook(queue, ScidrAnomalyCodebook(hostile))
  display <- ScidrReviewNarratives(view)$Issues
  date <- display[is.na(display$OtherVariable), ]
  expect_equal(jsonlite::fromJSON(date$CodebookCategoryKeys[1]), hostile$Category)
  html <- htmltools::renderTags(ExploreDataAnomalies(queue, codebook = hostile, IncludePlots = FALSE))$html
  expect_true(grepl("&lt;img", html, fixed = TRUE))
  skip_if_not_installed("xml2")
  document <- xml2::read_html(html)
  expect_length(xml2::xml_find_all(document, ".//img[@src='x']"), 0)
})

test_that("duplicate codebook rows that disagree only outside the read fields are accepted", {
  # A data dictionary carries value ranges, counts and notes. Two rows for one
  # variable that differ only there say nothing contradictory about the variable.
  codebook <- data.frame(
    Variable = c("age", "age", "sex"),
    Label = c("Age at visit", "Age at visit", "Sex"),
    Category = c("Demographics", "Demographics", "Demographics"),
    NMissing = c(3L, 7L, 0L),
    Values = c("18-90", "18-85", "F/M"))
  resolved <- SciDataReportR:::ScidrAnomalyCodebook(codebook)
  expect_equal(nrow(resolved), 2L)
  expect_equal(resolved$Variable, c("age", "sex"))
  expect_equal(resolved$Label[resolved$Variable == "age"], "Age at visit")
  expect_false(any(c("NMissing", "Values") %in% names(resolved)))
})

test_that("duplicate rows are merged when one of them simply leaves a field blank", {
  codebook <- data.frame(
    Variable = c("age", "age", "age"),
    Label = c("Age at visit", " Age at visit ", NA),
    Category = c(NA, "", "Demographics"))
  # Whitespace is not a disagreement, and a blank states nothing rather than
  # contradicting the row that does state something.
  resolved <- SciDataReportR:::ScidrAnomalyCodebook(codebook)
  expect_equal(nrow(resolved), 1L)
  expect_equal(resolved$Label, "Age at visit")
  expect_equal(resolved$Category, "Demographics")
})

test_that("a real codebook conflict names the field and the competing values", {
  codebook <- data.frame(
    Variable = c("age", "age", "sex", "sex"),
    Label = c("Age at visit", "Age at enrolment", "Sex", "Sex"),
    Category = c("Demographics", "Demographics", "Demographics", "Biology"))
  expect_error(SciDataReportR:::ScidrAnomalyCodebook(codebook),
    "Label: Age at visit vs Age at enrolment")
  expect_error(SciDataReportR:::ScidrAnomalyCodebook(codebook),
    "Category: Demographics vs Biology")
  # Only one field disagrees here, so only that one is named.
  one <- data.frame(Variable = c("age", "age"), Label = c("Age at visit", "Age at enrolment"),
    Category = c("Demographics", NA))
  expect_error(SciDataReportR:::ScidrAnomalyCodebook(one), "age \\(Label: ")
  expect_false(grepl("Category", tryCatch(SciDataReportR:::ScidrAnomalyCodebook(one),
    error = function(e) conditionMessage(e)), fixed = TRUE))
})

test_that("a long list of conflicts is summarised rather than printed in full", {
  codebook <- data.frame(Variable = rep(paste0("v", 1:8), each = 2),
    Label = paste0("label", seq_len(16)), Category = "Clinical")
  message_text <- tryCatch(SciDataReportR:::ScidrAnomalyCodebook(codebook),
    error = function(e) conditionMessage(e))
  expect_match(message_text, "and 3 more")
  expect_false(grepl("v8", message_text, fixed = TRUE))
})
