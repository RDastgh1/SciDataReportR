ScidrReviewSnapshot <- function(object) {
  evidence <- list(ReviewData = object$ReviewData, Results = object$Results,
    CellResults = object$CellResults, DateResults = object$DateResults,
    IntervalResults = object$IntervalResults, Variables = object$Metadata$Variables,
    Methods = object$Metadata$Methods, Thresholds = object$Metadata$Thresholds)
  temporal <- object$Metadata$TimeConversions
  if (isTRUE(nrow(temporal) > 0)) {
    # Only new temporal snapshots include temporal modeling semantics. Legacy hashes
    # and presentation-only label/category overrides retain their original identity.
    temporal <- temporal[c("Variable", "TimeType", "TimeUnit", "TimeFormat")]
    for (i in seq_len(nrow(temporal))) {
      x <- object$ReviewData[[temporal$Variable[i]]]
      if (inherits(x, "difftime")) temporal$TimeUnit[i] <- as.character(attr(x, "units"))
      if (!is.character(x)) temporal$TimeFormat[i] <- NA_character_
    }
    evidence$TemporalModel <- temporal
  }
  ScidrReviewHash(evidence)
}

ScidrReviewCheckDecisions <- function(df_Decisions) {
  for (column in c("ReviewStatus", "Priority")) {
    values <- as.character(df_Decisions[[column]])
    allowed <- if (column == "ReviewStatus") ScidrReviewStatuses() else c("Unassigned", "High", "Medium", "Low")
    if (any(!is.na(values) & nzchar(trimws(values)) & !values %in% allowed)) stop("Invalid ", column, "; use the documented dropdown choices.")
  }
  dates <- ScidrReviewFormat(df_Decisions$ReviewedOn)
  dates <- dates[!is.na(dates) & nzchar(trimws(dates))]
  parsed <- suppressWarnings(as.Date(dates, format = "%Y-%m-%d"))
  if (length(dates) && any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", dates) | is.na(parsed) | as.character(parsed) != dates)) stop("ReviewedOn must use valid YYYY-MM-DD dates.")
}

ScidrReviewValidate <- function(object) {
  if (!inherits(object, "AnomalyReviewQueue") || !all(c("Issues", "Evidence", "Summary", "Metadata", "Sources", "Diagnostics") %in% names(object))) stop("Provide an AnomalyReviewQueue from CreateAnomalyReviewQueue().")
  if (!identical(object$Metadata$FormatVersion, "1")) stop("Unsupported review queue format version.")
  if (!all(ScidrReviewFields() %in% names(object$Issues)) || anyNA(object$Issues$IssueID) || anyDuplicated(object$Issues$IssueID)) stop("Invalid queue issue IDs or review columns.")
  if (!identical(names(object$Sources), names(object$Metadata$Snapshots))) stop("Queue snapshot manifest is invalid.")
  for (analysis in names(object$Sources)) if (!identical(ScidrReviewSnapshot(object$Sources[[analysis]]), unname(object$Metadata$Snapshots[analysis]))) stop("Source snapshot changed; rebuild the queue before review.")
  fingerprint <- ScidrReviewHash(list(object$Issues[setdiff(names(object$Issues), ScidrReviewFields())], object$Evidence,
    as.list(object$Metadata$Snapshots)))
  if (!identical(fingerprint, object$Metadata$QueueFingerprint)) stop("Queue evidence changed; rebuild the queue rather than editing original findings.")
  ScidrReviewCheckDecisions(object$Issues)
  invisible(TRUE)
}

ScidrReviewFlat <- function(df_Data) {
  ScidrReviewRequire("jsonlite")
  df_Flat <- as.data.frame(df_Data)
  df_Flat[] <- lapply(df_Flat, function(x) {
    values <- if (is.list(x)) vapply(x, function(value) as.character(jsonlite::toJSON(value,
      auto_unbox = TRUE, na = "null", null = "null", digits = NA, POSIXt = "ISO8601")), character(1)) else ScidrReviewFormat(x)
    values[is.na(values)] <- ""
    unname(values)
  })
  rownames(df_Flat) <- NULL
  df_Flat
}

ScidrReviewHeaders <- function(columns) {
  labels <- c(IssueID = "Issue ID", EvidenceID = "Evidence ID", SourceRow = "Source row", OtherRow = "Other source row",
    StudyID = "Study ID", OtherStudyID = "Other study ID", Variable = "Variable key", Label = "Measure",
    OtherVariable = "Other variable key", OtherLabel = "Compared measure", Observed = "Recorded value",
    OtherObserved = "Compared value", ReviewCategory = "Review category",
    ExtremeValue = "Extreme value", ExtremeRatio = "Extreme ratio", WeakSupport = "Weak support", Explanation = "Why review",
    SuggestedCheck = "Suggested check", EvidenceCount = "Evidence count", ReviewStatus = "Review status",
    ProposedCorrection = "Proposed correction", SourceChecked = "Source checked", ReviewedOn = "Reviewed on",
    ReferenceN = "Reference observations", ReferenceLower = "Reference lower limit", ReferenceUpper = "Reference upper limit",
    ElapsedUnit = "Elapsed time unit", ActualIntervalDays = "Actual interval (days)", GroupStatistics = "Group statistics", NativeEvidence = "Original evidence",
    ExplanationMethod = "Explanation method", ScoreDecrease = "Score decrease after replacement",
    ReplacementScoreMean = "Mean replacement score", ReplacementScoreMin = "Minimum replacement score",
    ReplacementScoreMax = "Maximum replacement score", ReplacementN = "Replacement count",
    .Row = "Source row", ContributionID = "Contribution ID", Contribution = "Mean signed contribution",
    SimulationSD = "Simulation SD", SimulationMin = "Simulation minimum", SimulationMax = "Simulation maximum",
    PositiveRank = "Positive contribution rank", PositiveSensitivityRank = "Positive screening rank", AdditionalEvidence = "Additional evidence issue IDs (not confirmation)",
    Baseline = "Background baseline", SumContributions = "Sum of contributions", ReconstructionResidual = "Reconstruction residual",
    CodebookCategory = "Codebook category", CodebookCategoryKeys = "Category filter keys",
    CategoryBasis = "Category assignment", CategoryVariables = "Columns supplying categories", ConditionCategories = "Comparison-variable categories",
    WhyFlagged = "Why flagged", ComparisonDetails = "Comparison details",
    SuggestedReviewOrder = "Suggested review order", ReviewTier = "Review tier", OrderRationale = "Order rationale",
    ScorePercentile = "Within-method percentile (0-1)", ReviewReferenceN = "Review reference count", FindingType = "Anomaly type", Methods = "Supporting methods")
  out <- columns
  matched <- columns %in% names(labels)
  out[matched] <- unname(labels[columns[matched]])
  out
}

ScidrReviewManifest <- function(object) data.frame(FormatVersion = object$Metadata$FormatVersion,
  Analysis = names(object$Metadata$Snapshots), Snapshot = unname(object$Metadata$Snapshots),
  QueueFingerprint = object$Metadata$QueueFingerprint, stringsAsFactors = FALSE)

ScidrReviewReadme <- function(object) {
  df_Readme <- data.frame(Topic = c("Purpose", "Identity", "Editable columns", "Decisions",
  "Priority", "Proposed correction", "Dates", "Evidence", "Coverage", "Import", "Blank fields", "Limits"),
  Guidance = c("A worklist of review candidates, not confirmed errors. No source data are changed.",
    "Issue IDs and original row references identify findings within this detection snapshot. Study IDs are optional. Never move decisions to a changed or reordered dataset by ID alone.",
    paste(ScidrReviewHeaders(ScidrReviewFields()), collapse = ", "),
    paste(ScidrReviewStatuses(), collapse = " / "),
    "Reviewer-assigned: Unassigned / High / Medium / Low. Statistical rarity does not set clinical urgency.",
    "Free text. Specify each affected variable and proposed value if multiple fields need correction. Nothing is automatically applied.",
    "Enter Reviewed on as YYYY-MM-DD. Original midnight-only date values display without fractional seconds; original timestamps/time zones remain in the R source object.",
    "Evidence links each issue to methods, conditions, counts and limits. Different date pairs remain distinct. Agreement is not independent confirmation. Variable keys are visible alongside labels. Tracking IDs, original backend explanations and structured evidence are hidden columns; unhide them to inspect without editing. Conditional pattern means unusual within a learned comparison group, not a clinical rule or confirmed error.",
    "Coverage by analysis and method is listed below. Zero flags without a cutoff is not proof of normality.",
    "Use ImportAnomalyReviewDecisions(queue, file). Keep all columns; you may sort rows and remove issue rows to return a subset. Original evidence and the hidden snapshot manifest must remain unchanged.",
    "Blank imported reviewer fields leave existing decisions unchanged. Explicit Unreviewed resets status. Omitted issue rows keep their decisions.",
    "Conditional detection currently exposes at most one finding per record. Diagnostics are coverage information, not correction instructions. The workbook is not encrypted or a security audit signature."),
  stringsAsFactors = FALSE)

  df_Readme <- rbind(df_Readme, data.frame(Topic = "Contributions", Guidance = "Available isolation explanations are linked by analysis, source row and issue IDs. The Explanation method distinguishes sensitivity from SHAP. For sensitivity, score decrease is original score minus mean score after replacing one column with reference values; positive suggests a column to inspect. Effects are not additive, replacement-score ranges are not confidence intervals, and replacements can create unrealistic combinations. For SHAP: Positive signed contributions raise the model score relative to the sampled background; negative ones lower it. Simulation SD/range describe Monte Carlo variability, not error probability. Correlated/derived columns can share attribution. These explain the score, not errors or corrections. Contributions and their separate explanation fingerprint are verified on import; older workbooks may omit this sheet."))
  df_Readme <- rbind(df_Readme, data.frame(Topic = "Codebook categories", Guidance = "Codebook categories describe affected fields. Suggested inspection categories for isolation records use the five largest positive screening effects or saved SHAP contributions, not identified problem fields. Labels and categories are presentation metadata and do not change source values, models or issue identities. Review decisions are edited here in Excel; the explorer displays imported decisions in finding details."))
  if (nrow(object$Summary)) {
    df_Coverage <- data.frame(Topic = paste("Coverage:", object$Summary$Analysis, object$Summary$Method, sep = " / "),
      Guidance = paste0(object$Summary$Flagged, " flagged of ", object$Summary$Evaluated, " evaluated records; ",
        object$Summary$Records, " source records; ", object$Summary$Unavailable, " unavailable record checks."))
    df_Readme <- rbind(df_Readme, df_Coverage)
  }
  df_Readme
}

#' Export an anomaly review queue for collaborator decisions
#'
#' Write Review, Evidence, Read me, Contributions and a hidden manifest.
#' The optional Contributions sheet has a separately validated explanation fingerprint.
#' Original values and evidence are read-only by convention and verified on import;
#' yellow review columns are editable. No source data are changed.
#' @param ReviewObj An AnomalyReviewQueue from [CreateAnomalyReviewQueue()].
#' @param file Destination .xlsx path in an existing directory.
#' @param overwrite Allow replacement of an existing workbook? Default FALSE.
#' @param layout `"compact"` (the default) or `"full"`. Compact is built for
#'   scanning a queue: **one sheet and nothing else**, one line per finding, and
#'   only the columns a reviewer reads (order, record, measure, value, anomaly
#'   type, tier, codebook category, extreme value, why flagged) followed by the
#'   editable decision fields. A column empty for every finding is dropped, and
#'   the codebook category is dropped unless a codebook assigned a category to
#'   at least one screened column -- supply `codebook` to
#'   [DetectDataAnomalies()] or [CreateAnomalyReviewQueue()] to sort and filter
#'   the sheet by category. The long
#'   explanations, the Evidence sheet and the snapshot manifest are **not
#'   written**, so a compact workbook cannot be read back by
#'   [ImportAnomalyReviewDecisions()] -- it is a deliverable, not a round trip.
#'   Use `"full"` when decisions must be imported: every sheet, every column,
#'   cells wrapped over as many lines as their text needs.
#' @return The destination path, invisibly. Requires optional openxlsx and jsonlite.
#' @details Reviewer fields are free text except dropdown status and priority.
#'   IDs, observed values and proposed corrections are exported as text to preserve
#'   leading zeros and prevent automatic date/number reinterpretation. The workbook
#'   is not encrypted. Keep the R queue object for importing decisions.
#' @export
ExportAnomalyReviewQueue <- function(ReviewObj, file, overwrite = FALSE,
    layout = c("compact", "full")) {
  layout <- match.arg(layout)
  compact <- identical(layout, "compact")
  ScidrReviewRequire("openxlsx")
  ScidrReviewRequire("jsonlite")
  ScidrReviewValidate(ReviewObj)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file) || !grepl("\\.xlsx$", file, ignore.case = TRUE)) stop("file must be one .xlsx path.")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) stop("overwrite must be TRUE or FALSE.")
  if (file.exists(file) && !overwrite) stop("Workbook already exists; choose another file or set overwrite = TRUE.")
  if (!dir.exists(dirname(file))) stop("Destination directory does not exist.")
  contributions <- ScidrReviewContributions(ReviewObj)
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = 11, fontName = "Arial")
  display <- ScidrReviewNarratives(ReviewObj)
  sheets <- list(Review = ScidrReviewFlat(display$Issues), Evidence = ScidrReviewFlat(display$Evidence),
    `Read me` = ScidrReviewReadme(ReviewObj), Manifest = ScidrReviewManifest(ReviewObj),
    Contributions = ScidrReviewFlat(contributions))
  sheets$Manifest$ExplanationFingerprint <- vapply(names(ReviewObj$Sources), function(analysis) {
    source <- ReviewObj$Sources[[analysis]]
    if (is.null(source$ExplanationMetadata)) "" else ScidrAnomalyExplanationFingerprint(source)
  }, character(1))
  header <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#263F53", textDecoration = "bold", wrapText = TRUE, valign = "center")
  # Unwrapped text keeps every row one line high; the full value is still in the
  # cell and shows in the formula bar when a row is selected.
  text <- openxlsx::createStyle(numFmt = "@", wrapText = !compact, valign = if (compact) "center" else "top")
  editable <- openxlsx::createStyle(numFmt = "@", fgFill = "#FFF2CC", wrapText = !compact, valign = if (compact) "center" else "top")
  if (compact) sheets <- sheets["Review"]
  # A codebook that named no category for any screened column leaves nothing to
  # sort or filter on, so the compact sheet omits the column entirely.
  categorized <- any(vapply(names(ReviewObj$Sources), function(analysis)
    any(ScidrReviewCatalog(ReviewObj, analysis)$CategorySource %in% c("codebook", "saved")), logical(1)))
  for (sheet in names(sheets)) {
    df_Data <- sheets[[sheet]]
    if (compact) {
      names(df_Data) <- ScidrReviewHeaders(names(df_Data))
      keep <- intersect(ScidrReviewCompactColumns(), names(df_Data))
      # Only three columns are dropped, and only because they do not apply to
      # this queue: no ID variable was supplied, no finding was extreme, or no
      # codebook assigned a category. The category is decided from the resolved
      # catalogs rather than from the written cells, because an unsupplied
      # codebook still fills those cells -- "Uncategorized" for a column and a
      # whole-record placeholder otherwise -- which a reviewer can neither sort
      # nor filter on. Every other column keeps its place even when this queue
      # happens not to fill it, so the sheet has the same shape from one run to
      # the next.
      optional <- ScidrReviewHeaders(c("StudyID", "ExtremeValue"))
      if (!categorized) keep <- setdiff(keep, ScidrReviewHeaders("CodebookCategory"))
      used <- vapply(keep, function(column) !column %in% optional ||
        any(!is.na(df_Data[[column]]) & nzchar(as.character(df_Data[[column]]))), logical(1))
      df_Data <- df_Data[keep[used]]
    }
    if (any(vapply(df_Data, function(x) any(nchar(as.character(x), type = "chars") > 32767, na.rm = TRUE), logical(1)))) stop("Workbook cell exceeds Excel's 32767-character limit; reduce context/evidence size before export.")
    if (!compact && sheet %in% c("Review", "Evidence", "Contributions")) {
      preferred <- c("SuggestedReviewOrder", "ReviewTier", "OrderRationale", "IssueID", "EvidenceID", "Analysis", "StudyID", "SourceRow", ".Row", "Label", "Variable", "OtherLabel", "OtherVariable", "Observed", "OtherObserved", "CodebookCategory", "WhyFlagged", "ComparisonDetails", "ReviewStatus", "Priority", "SuggestedCheck")
      df_Data <- df_Data[c(intersect(preferred, names(df_Data)), setdiff(names(df_Data), preferred))]
      names(df_Data) <- ScidrReviewHeaders(names(df_Data))
    }
    openxlsx::addWorksheet(wb, sheet, visible = compact || sheet != "Manifest", gridLines = FALSE)
    if (nrow(df_Data)) openxlsx::writeDataTable(wb, sheet, df_Data, tableStyle = "TableStyleMedium2") else openxlsx::writeData(wb, sheet, df_Data)
    openxlsx::addStyle(wb, sheet, header, rows = 1, cols = seq_along(df_Data), gridExpand = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 32)
    if (nrow(df_Data)) {
      openxlsx::addStyle(wb, sheet, text, rows = seq_len(nrow(df_Data)) + 1, cols = seq_along(df_Data), gridExpand = TRUE)
    }
    openxlsx::freezePane(wb, sheet, firstActiveRow = 2,
      firstActiveCol = if (sheet != "Review") 3 else if (compact) 4 else 6)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(df_Data), widths = if (sheet == "Read me") c(24, 110) else 23)
    if (sheet %in% c("Review", "Evidence", "Contributions")) {
      wide <- which(names(df_Data) %in% c("Measure", "Compared measure", "Why review", "Why flagged", "Comparison details", "Suggested check", "Proposed correction", "Notes"))
      if (length(wide)) openxlsx::setColWidths(wb, sheet, cols = wide, widths = 48)
      hidden <- which(names(df_Data) %in% c("Issue ID", "Evidence ID", "Original evidence", "Group statistics", "Conditions", "Why review", "Contribution ID", "Category filter keys"))
      if (sheet == "Review") hidden <- unique(c(hidden, which(names(df_Data) %in% c("Study ID", "Other study ID", "Other source row", "Context") &
        vapply(df_Data, function(x) all(is.na(x) | x == ""), logical(1)))))
      if (length(hidden)) openxlsx::setColWidths(wb, sheet, cols = hidden, widths = 24, hidden = TRUE)
    }
    if (nrow(df_Data) && sheet != "Manifest") {
      widths <- if (sheet == "Read me") c(24, 110) else rep(23, ncol(df_Data))
      hidden_cols <- integer()
      if (sheet %in% c("Review", "Evidence", "Contributions")) {
        widths[wide] <- 48
        hidden_cols <- hidden
        short <- which(names(df_Data) %in% c("Source row", "Other source row", "Study ID", "Other study ID", "Priority"))
        if (length(short)) {
          widths[short] <- 14
          openxlsx::setColWidths(wb, sheet, cols = short, widths = 14, hidden = short %in% hidden_cols)
        }
        decision <- which(names(df_Data) == "Review status")
        if (length(decision)) {
          widths[decision] <- 27
          openxlsx::setColWidths(wb, sheet, cols = decision, widths = 27)
        }
      }
      visible_cols <- setdiff(seq_len(ncol(df_Data)), hidden_cols)
      heights <- vapply(seq_len(nrow(df_Data)), function(i) {
        lines <- vapply(visible_cols, function(j) {
          value <- as.character(df_Data[[j]][i])
          if (is.na(value)) return(1)
          sum(pmax(1, ceiling(nchar(strsplit(value, "\n", fixed = TRUE)[[1]]) / max(widths[j] - 3, 1))))
        }, numeric(1))
        if (compact) 22 else min(400, max(44, max(lines) * 14 + 12))
      }, numeric(1))
      openxlsx::setRowHeights(wb, sheet, rows = seq_len(nrow(df_Data)) + 1, heights = heights)
    }
    if (sheet == "Review" && nrow(df_Data)) {
      cols <- match(ScidrReviewHeaders(ScidrReviewFields()), names(df_Data))
      openxlsx::addStyle(wb, sheet, editable, rows = seq_len(nrow(df_Data)) + 1, cols = cols, gridExpand = TRUE, stack = TRUE)
      for (column in c("ReviewStatus", "Priority")) {
        choices <- if (column == "ReviewStatus") ScidrReviewStatuses() else c("Unassigned", "High", "Medium", "Low")
        openxlsx::dataValidation(wb, sheet, cols = match(ScidrReviewHeaders(column), names(df_Data)),
          rows = seq_len(nrow(df_Data)) + 1, type = "list", value = paste0('"', paste(choices, collapse = ","), '"'), allowBlank = TRUE, showErrorMsg = TRUE)
      }
    }
  }
  openxlsx::saveWorkbook(wb, file, overwrite = overwrite)
  invisible(normalizePath(file, mustWork = TRUE))
}

#' Import collaborator decisions into an existing anomaly review snapshot
#'
#' Validate a review workbook and copy only editable review fields into a queue.
#' Original source data and detection evidence are never changed.
#' @inheritParams ExportAnomalyReviewQueue
#' @return An updated AnomalyReviewQueue. The input object is unchanged. Omitted
#'   issues and blank reviewer cells preserve existing decisions; an explicit
#'   Unreviewed resets status. ReviewedOn must be a valid YYYY-MM-DD string.
#' @details Import requires the matching snapshot manifest, unique known issue IDs,
#'   all original columns, and unchanged original findings and evidence. Sorting
#'   rows is allowed. Review may contain a subset of issues; Evidence must remain
#'   complete. Contributions, when present, must match their separate explanation
#'   fingerprint and contents. Older four-sheet workbooks remain supported.
#'   This checks accidental mismatches, not malicious tampering.
#' @export
ImportAnomalyReviewDecisions <- function(ReviewObj, file) {
  ScidrReviewRequire("openxlsx")
  ScidrReviewRequire("jsonlite")
  ScidrReviewValidate(ReviewObj)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !file.exists(file)) stop("file must name an existing review workbook.")
  wb <- openxlsx::loadWorkbook(file)
  if (!all(c("Review", "Evidence", "Read me", "Manifest") %in% names(wb)))
    stop("Workbook is missing required review sheets. A compact workbook carries only the Review sheet and cannot be imported; export with layout = \"full\" when decisions need to be read back.")
  legacy_display <- NULL
  previous_wording <- NULL
  previous_legacy_wording <- NULL
  LegacyDisplay <- function() {
    if (is.null(legacy_display)) {
      legacy <- ReviewObj
      legacy$Metadata$VariableMetadata <- NULL
      for (analysis in names(legacy$Sources)) legacy$Sources[[analysis]]$Metadata$VariableMetadata <- NULL
      legacy_display <<- ScidrReviewNarratives(legacy)
    }
    legacy_display
  }
  Read <- function(sheet, expected) {
    df_Data <- openxlsx::read.xlsx(wb, sheet, check.names = FALSE, sep.names = " ", skipEmptyRows = FALSE, skipEmptyCols = FALSE,
      detectDates = FALSE, na.strings = character())
    if (sheet == "Manifest" && "Contributions" %in% names(wb) && !"ExplanationFingerprint" %in% names(df_Data)) stop("Missing explanation fingerprint.")
      optional <- c("ExplanationFingerprint", ScidrReviewHeaders(ScidrReviewPresentationFields()), ScidrReviewHeaders(c("WhyFlagged", "ComparisonDetails", "SuggestedReviewOrder", "ReviewTier", "OrderRationale", "ScorePercentile", "ReviewReferenceN", "FindingType", "Methods", "ExtremeValue", "ExtremeRatio", "WeakSupport")),
      # Workbooks written before these columns were named carry the raw field name.
      "ExtremeValue", "ExtremeRatio", "WeakSupport")
    if (anyDuplicated(names(df_Data)) || !all(setdiff(names(expected), optional) %in% names(df_Data)) || any(!names(df_Data) %in% c(names(expected), optional))) stop(sheet, " columns changed; retain original headers and use Notes for annotations.")
    df_Data <- df_Data[intersect(names(expected), names(df_Data))]
    # Labels and categories are display metadata, never imported as source evidence.
    presentation <- setdiff(ScidrReviewPresentationFields(), c("WhyFlagged", "ComparisonDetails"))
    if (sheet %in% c("Review", "Evidence")) {
      key <- if (sheet == "Review") "Issue ID" else "Evidence ID"
      index <- match(df_Data[[key]], expected[[key]])
      for (column in intersect(c("Why flagged", "Comparison details"), names(df_Data))) {
        target <- expected[[column]][index]
        differs <- !is.na(index) & !is.na(df_Data[[column]]) & df_Data[[column]] != target
        differs[is.na(differs)] <- FALSE
        if (any(differs)) {
          if (is.null(previous_wording)) previous_wording <<- ScidrReviewNarratives(ReviewObj, concise = FALSE)
          if (is.null(previous_legacy_wording)) {
            legacy <- ReviewObj
            legacy$Metadata$VariableMetadata <- NULL
            for (analysis in names(legacy$Sources)) legacy$Sources[[analysis]]$Metadata$VariableMetadata <- NULL
            previous_legacy_wording <<- ScidrReviewNarratives(legacy, concise = FALSE)
          }
          for (version in list(previous_wording, LegacyDisplay(), previous_legacy_wording)) {
            old <- version[[if (sheet == "Review") "Issues" else "Evidence"]]
            names(old) <- ScidrReviewHeaders(names(old))
            old_index <- match(df_Data[[key]], old[[key]])
            compatible <- differs & !is.na(old_index) & df_Data[[column]] == old[[column]][old_index]
            compatible[is.na(compatible)] <- FALSE
            df_Data[[column]][compatible] <- target[compatible]
          }
        }
      }
    }
    df_Data <- df_Data[setdiff(names(df_Data), ScidrReviewHeaders(presentation))]
    ScidrReviewFlat(df_Data)
  }
  Compare <- function(actual, expected, key, sheet) {
    if (anyDuplicated(actual[[key]]) || any(!nzchar(actual[[key]])) || !setequal(actual[[key]], expected[[key]])) stop(sheet, " IDs do not match the snapshot.")
    expected <- expected[match(actual[[key]], expected[[key]]), names(actual), drop = FALSE]
    rownames(expected) <- NULL
    if (!identical(actual, expected)) stop(sheet, " original evidence or snapshot manifest changed.")
  }
  manifest <- ScidrReviewFlat(ScidrReviewManifest(ReviewObj))
  if ("Contributions" %in% names(wb)) {
    expected_contributions <- ScidrReviewFlat(ScidrReviewContributions(ReviewObj))
    names(expected_contributions) <- ScidrReviewHeaders(names(expected_contributions))
    Compare(Read("Contributions", expected_contributions), expected_contributions, "Contribution ID", "Contributions")
    manifest$ExplanationFingerprint <- vapply(names(ReviewObj$Sources), function(analysis) {
      source <- ReviewObj$Sources[[analysis]]
      if (is.null(source$ExplanationMetadata)) "" else ScidrAnomalyExplanationFingerprint(source)
    }, character(1))
  } else {
    # A legacy workbook may have no explanation sheet or fingerprint.
    original_read <- Read
    Read <- function(sheet, expected) {
      if (sheet == "Manifest") {
        df_Manifest <- openxlsx::read.xlsx(wb, sheet, check.names = FALSE, sep.names = " ", na.strings = character())
        df_Manifest$ExplanationFingerprint <- NULL
        if (!setequal(names(df_Manifest), names(expected))) stop("Manifest columns changed.")
        return(ScidrReviewFlat(df_Manifest[names(expected)]))
      }
      original_read(sheet, expected)
    }
  }
  Compare(Read("Manifest", manifest), manifest, "Analysis", "Manifest")
  display <- ScidrReviewNarratives(ReviewObj)
  evidence <- ScidrReviewFlat(display$Evidence)
  names(evidence) <- ScidrReviewHeaders(names(evidence))
  Compare(Read("Evidence", evidence), evidence, "Evidence ID", "Evidence")
  expected <- ScidrReviewFlat(display$Issues)
  names(expected) <- ScidrReviewHeaders(names(expected))
  review <- Read("Review", expected)
  if (any(!nzchar(review$`Issue ID`)) || anyDuplicated(review$`Issue ID`) || any(!review$`Issue ID` %in% expected$`Issue ID`)) stop("Review contains duplicate, blank or unknown issue IDs.")
  index <- match(review$`Issue ID`, expected$`Issue ID`)
  immutable <- setdiff(names(review), ScidrReviewHeaders(ScidrReviewFields()))
  expected_subset <- expected[index, immutable, drop = FALSE]
  rownames(expected_subset) <- NULL
  if (!identical(review[immutable], expected_subset)) stop("Review original evidence changed; edit only reviewer decision columns.")
  fields <- c("IssueID", ScidrReviewFields())
  review <- review[ScidrReviewHeaders(fields)]
  names(review) <- fields
  ScidrReviewCheckDecisions(review)
  index <- match(review$IssueID, ReviewObj$Issues$IssueID)
  updated <- ReviewObj
  for (column in ScidrReviewFields()) {
    values <- review[[column]]
    nonblank <- !is.na(values) & nzchar(trimws(values))
    updated$Issues[[column]][index[nonblank]] <- values[nonblank]
  }
  updated
}

#' Summarize collaborator review decisions
#'
#' Report what reviewers decided, after [ImportAnomalyReviewDecisions()] has
#' copied a returned workbook into a queue. Read-only: source data is never
#' modified and no corrections are applied, so acting on the decisions remains
#' the analyst's judgement.
#' @param ReviewObj An `AnomalyReviewQueue`, normally one returned by
#'   [ImportAnomalyReviewDecisions()].
#' @param by Grouping for the counts: `status` (default), `variable`,
#'   `analysis`, `priority`, or `reviewer`.
#' @return A list with `Counts` (a tibble of issue counts by the requested
#'   grouping and review status), `Corrections` (one row per issue with a
#'   proposed correction, with its record, field, recorded value, proposed
#'   value, reviewer and notes), and `Progress` (issues reviewed, outstanding
#'   and the reviewed fraction). Unreviewed issues are counted, never dropped.
#' @family data anomalies
#' @seealso [ImportAnomalyReviewDecisions()], [CreateAnomalyReviewQueue()]
#' @examples
#' set.seed(24)
#' df_Example <- data.frame(ID = 1:60, Visit = as.Date("2024-01-01") + c(0:58, 900))
#' result <- DetectDataAnomalies(df_Example, "Visit", id_var = "ID", methods = "dates")
#' queue <- CreateAnomalyReviewQueue(result)
#' SummarizeAnomalyDecisions(queue)$Progress
#' @export
SummarizeAnomalyDecisions <- function(ReviewObj, by = c("status", "variable", "analysis", "priority", "reviewer")) {
  if (!inherits(ReviewObj, "AnomalyReviewQueue")) stop("ReviewObj must be an AnomalyReviewQueue.", call. = FALSE)
  by <- match.arg(by)
  df_Issues <- ReviewObj$Issues
  statuses <- ScidrReviewStatuses()
  status <- if (!nrow(df_Issues)) character() else
    ifelse(is.na(df_Issues$ReviewStatus) | !nzchar(df_Issues$ReviewStatus), "Unreviewed", df_Issues$ReviewStatus)
  reviewed <- status != "Unreviewed"
  group <- switch(by,
    status = status,
    variable = if (!nrow(df_Issues)) character() else ifelse(is.na(df_Issues$Variable),
      "No attributed field", df_Issues$Variable),
    analysis = df_Issues$Analysis,
    priority = if (!nrow(df_Issues)) character() else ifelse(is.na(df_Issues$Priority) |
      !nzchar(df_Issues$Priority), "No priority assigned", df_Issues$Priority),
    reviewer = if (!nrow(df_Issues)) character() else ifelse(is.na(df_Issues$Reviewer) |
      !nzchar(df_Issues$Reviewer), "No reviewer recorded", df_Issues$Reviewer))
  df_Counts <- if (!nrow(df_Issues)) tibble::tibble(Group = character(), Status = character(), Issues = integer()) else
    tibble::as_tibble(stats::aggregate(list(Issues = rep(1L, length(group))),
      by = list(Group = group, Status = factor(status, levels = unique(c(statuses, status)))), FUN = sum))
  if (nrow(df_Counts)) {
    df_Counts$Status <- as.character(df_Counts$Status)
    df_Counts <- df_Counts[order(df_Counts$Group, df_Counts$Status), ]
  }
  names(df_Counts)[names(df_Counts) == "Group"] <- switch(by, status = "ReviewStatus",
    variable = "Variable", analysis = "Analysis", priority = "Priority", reviewer = "Reviewer")
  if (by == "status") df_Counts$Status <- NULL
  proposed <- if (!nrow(df_Issues)) logical() else
    !is.na(df_Issues$ProposedCorrection) & nzchar(df_Issues$ProposedCorrection)
  df_Corrections <- tibble::as_tibble(df_Issues[proposed, intersect(c("IssueID", "Analysis", "SourceRow",
    "StudyID", "Variable", "Label", "Observed", "ProposedCorrection", "ReviewStatus", "Priority",
    "Reviewer", "ReviewedOn", "Notes"), names(df_Issues)), drop = FALSE])
  list(Counts = df_Counts, Corrections = df_Corrections,
    Progress = tibble::tibble(Issues = nrow(df_Issues), Reviewed = sum(reviewed),
      Outstanding = sum(!reviewed), WithProposedCorrection = sum(proposed),
      ReviewedFraction = if (nrow(df_Issues)) sum(reviewed) / nrow(df_Issues) else NA_real_))
}

# The columns a reviewer actually scans: who, what, what it says, and the
# decision fields they fill in. Everything else stays in the sheet for the
# importer and for anyone who unhides it.
ScidrReviewCompactColumns <- function() c(
  ScidrReviewHeaders(c("SuggestedReviewOrder", "StudyID", "SourceRow", "Label", "Variable",
    "Observed", "FindingType", "ReviewTier", "CodebookCategory", "ExtremeValue", "WhyFlagged",
    ScidrReviewFields())))
