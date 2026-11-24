# A review finding is one row with explicit names and already validated scalar fields.
# Avoid repeated tibble expression/name processing when building large queues.
ScidrReviewOneRow <- function(...) structure(list(...), class = "data.frame", row.names = 1L)

ScidrReviewRequire <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install package '", package, "' for anomaly review queues.", call. = FALSE)
}

ScidrReviewHash <- function(x) {
  digest::digest(x, algo = "sha256", serializeVersion = 2)
}

ScidrReviewFields <- function() c("ReviewStatus", "Priority", "ProposedCorrection",
  "SourceChecked", "Reviewer", "ReviewedOn", "Notes")

ScidrReviewStatuses <- function() c("Unreviewed", "Confirmed correct", "Correction needed", "Definition clarification")

# A CHARSXP is R's internal string type, and it is never meant to surface in
# R-level code. The conditional backend builds its condition values in C++ and
# can hand one back, and then nothing in R can touch the value: inherits(),
# format(), paste() and as.character() all abort on it with "cannot have
# attributes on a CHARSXP". One such value in one condition of one finding would
# otherwise take down every page, plot and export built from the whole queue.
# deparse() is the only accessor that reads a CHARSXP, and it returns the text
# verbatim -- unquoted and unescaped, even for embedded quotes and backslashes --
# so the value is restored to the character string it should always have been.
ScidrRecoverValues <- function(x) {
  if (typeof(x) == "char") return(deparse(x))
  if (typeof(x) == "list") {
    out <- lapply(x, ScidrRecoverValues)
    attributes(out) <- attributes(x)
    return(out)
  }
  x
}

ScidrReviewFormat <- function(x) {
  if (typeof(x) == "char") x <- deparse(x)
  if (inherits(x, "Date")) return(as.character(x))
  if (inherits(x, "POSIXct")) {
    if (all(is.na(x) | format(x, "%H:%M:%OS6") == "00:00:00.000000")) return(format(x, "%Y-%m-%d"))
    return(format(x, "%Y-%m-%d %H:%M:%OS6 %Z"))
  }
  as.character(x)
}

ScidrReviewField <- function(df_Data, name, i, default = NA) {
  if (!name %in% names(df_Data)) return(default)
  if (is.list(df_Data[[name]])) df_Data[[name]][[i]] else df_Data[[name]][i]
}

ScidrReviewEmptyIssues <- function() tibble::tibble(IssueID = character(), Analysis = character(),
  SourceRow = integer(), OtherRow = integer(), StudyID = character(), OtherStudyID = character(),
  Context = character(), Variable = character(), Label = character(), OtherVariable = character(),
  OtherLabel = character(), Observed = character(), OtherObserved = character(),
  ReviewCategory = character(), Explanation = character(), SuggestedCheck = character(),
  EvidenceCount = integer(), ReviewStatus = character(), Priority = character(),
  ProposedCorrection = character(), SourceChecked = character(), Reviewer = character(),
  ReviewedOn = character(), Notes = character())

ScidrReviewEmptyEvidence <- function() tibble::tibble(IssueID = character(), EvidenceID = character(),
  Analysis = character(), SourceRow = integer(), OtherRow = integer(), Variable = character(),
  Label = character(), OtherVariable = character(), OtherLabel = character(), Observed = character(),
  OtherObserved = character(), Method = character(), Kind = character(), Score = double(), Rank = double(),
  Flag = logical(), Explanation = character(), SuggestedCheck = character(), ReferenceN = double(),
  ReferenceLower = character(), ReferenceUpper = character(), Expected = character(), Elapsed = double(), ElapsedUnit = character(), ActualIntervalDays = double(),
  Conditions = list(), GroupStatistics = list(), NativeEvidence = list())

ScidrReviewSummarize <- function(sources) {
  dplyr::bind_rows(lapply(names(sources), function(analysis) {
    df_Results <- sources[[analysis]]$Results
    df_Results %>% dplyr::group_by(.data$Method) %>% dplyr::summarise(
      Records = dplyr::n_distinct(.data$.Row),
      Evaluated = sum(.data$Status == "evaluated"),
      Flagged = sum(.data$Flag %in% TRUE & .data$Status == "evaluated"),
      Unavailable = sum(.data$Status == "unavailable"),
      .groups = "drop") %>% dplyr::mutate(Analysis = analysis,
        FlagFraction = ifelse(.data$Evaluated > 0, .data$Flagged / .data$Evaluated, NA_real_))
  }))
}

#' Build a collaborator review queue from anomaly evidence
#'
#' Convert flagged findings into specific, traceable issues without changing source
#' data. Study IDs are optional; source rows identify records within a detection
#' snapshot. Statistical rarity is kept separate from reviewer-assigned priority.
#' @param AnomalyObj A DataAnomalyResult or a uniquely named list of such results.
#'   Different analyses retain independent row namespaces.
#' @param id_var Optional retained ID column name. NULL uses each detector's ID
#'   metadata when available. An override can add retained IDs without refitting.
#' @param context_vars Optional character vector of additional retained columns
#'   to display. Recorded time and stratification columns are included automatically.
#' @param include_ranked Include high-scoring unflagged records in addition to
#'   native and cutoff-based flags? Default TRUE.
#' @param codebook Optional presentation codebook as in [DetectDataAnomalies()].
#'   Overrides matching saved labels/categories across analyses. Source detection
#'   objects, issue identities and evidence remain unchanged. The resolved catalog
#'   is retained in queue metadata for exploration and Excel.
#' @param top_n Positive integer limit per method and stratum for unclassified
#'   ranked candidates, only when include_ranked is TRUE.
#' @param variables Optional character vector of screened column names. NULL
#'   keeps every finding. Otherwise only findings attributed to a listed column
#'   are kept, together with whole-record findings, which are raised about no one
#'   field. Detection is unchanged: the models were fitted on the full selection,
#'   so a kept finding is still conditioned on every screened column, and coverage
#'   in Summary still describes the whole screen.
#' @return An AnomalyReviewQueue list with Issues, Evidence, Summary, Diagnostics,
#'   Metadata, and Sources (unaltered original detection objects). Issues contain
#'   editable ReviewStatus, Priority, ProposedCorrection, SourceChecked, Reviewer,
#'   ReviewedOn and Notes. No corrections are applied automatically.
#' @details Issue IDs are deterministic within a fingerprinted snapshot. Rows,
#'   selected values or detection evidence changing create a different snapshot;
#'   study IDs do not authorize transferring decisions between snapshots.
#'   Related methods on the same field set and rows are consolidated. Different
#'   date pairs stay separate. Method agreement is not independent confirmation.
#'   The explorer and Excel derive an unvalidated suggested review order from
#'   evidence tiers, within-method/stratum midrank score percentiles, reference
#'   counts and issue IDs. These presentation fields do not alter fingerprints,
#'   reviewer priority or error probabilities. Existing queues retain membership.
#'   Older conditional findings retain their original explanation; missing
#'   structured group statistics remain unavailable.
#' @examples
#' df_Example <- data.frame(Visit = as.Date("2020-01-01") + seq_len(40))
#' df_Example$Visit[1] <- as.Date("2099-01-01")
#' result <- DetectDataAnomalies(df_Example, "Visit", methods = "dates")
#' queue <- CreateAnomalyReviewQueue(result)
#' queue$Issues
#' @family data anomalies
#' @seealso [ExportAnomalyReviewQueue()], [SummarizeAnomalyDecisions()]
#' @export
CreateAnomalyReviewQueue <- function(AnomalyObj, id_var = NULL, context_vars = NULL,
    include_ranked = TRUE, top_n = 10, codebook = NULL, variables = NULL) {
  codebook <- ScidrAnomalyCodebook(codebook)
  if (!is.logical(include_ranked) || length(include_ranked) != 1L || is.na(include_ranked)) stop("include_ranked must be TRUE or FALSE.")
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != floor(top_n)) stop("top_n must be a positive integer.")
  for (arg in c("id_var", "context_vars")) {
    value <- get(arg)
    if (!is.null(value) && (!is.character(value) || !length(value) || anyNA(value) || any(!nzchar(value)) || anyDuplicated(value))) stop(arg, " must be NULL or unique column names.")
  }
  if (length(id_var) > 1L) stop("id_var must name one retained column.")
  sources <- if (inherits(AnomalyObj, "DataAnomalyResult")) list(Analysis = AnomalyObj) else AnomalyObj
  if (!is.list(sources) || !length(sources) || is.null(names(sources)) || anyNA(names(sources)) || any(!nzchar(names(sources))) || anyDuplicated(names(sources)) ||
      !all(vapply(sources, inherits, logical(1), "DataAnomalyResult"))) stop("Provide a detection result or a uniquely named list of detection results.")
  evidence <- list(ScidrReviewEmptyEvidence())
  snapshots <- ids <- contexts <- list()
  diagnostics <- list()
  for (analysis in names(sources)) {
    object <- sources[[analysis]]
    if (!all(c("ReviewData", "Results", "Metadata") %in% names(object))) stop("Incomplete detection result: ", analysis)
    df_Data <- object$ReviewData
    if (!".Row" %in% names(df_Data) || anyNA(df_Data$.Row) || anyDuplicated(df_Data$.Row) ||
        !is.numeric(df_Data$.Row) || any(!is.finite(df_Data$.Row)) || any(df_Data$.Row < 1 | df_Data$.Row != floor(df_Data$.Row))) stop("ReviewData requires unique positive original row keys.")
    selected_id <- if (is.null(id_var)) object$Metadata$IDVariable else id_var
    if (is.null(selected_id)) selected_id <- character()
    if (any(!c(selected_id, context_vars) %in% names(df_Data))) stop("ID/context columns must be retained in ReviewData: ", analysis)
    ids[[analysis]] <- selected_id
    contexts[[analysis]] <- unique(c(object$Metadata$TimeVariable, object$Metadata$StrataVariables, context_vars))
    snapshots[[analysis]] <- ScidrReviewSnapshot(object)
    AddDiagnostic <- function(issue, detail) {
      diagnostics[[length(diagnostics) + 1L]] <<- tibble::tibble(Analysis = analysis, Issue = issue, Detail = detail)
    }
    if (!length(selected_id)) AddDiagnostic("row_identity", "No study ID supplied; row references are valid only within this snapshot.")
    else {
      if (anyNA(df_Data[[selected_id]])) AddDiagnostic("missing_id", "Missing IDs retained; row keys identify records.")
      if (anyDuplicated(df_Data[[selected_id]][!is.na(df_Data[[selected_id]])])) AddDiagnostic("duplicate_id", "Duplicate IDs retained as distinct source rows.")
    }
    if (isTRUE(nrow(object$Diagnostics) > 0)) diagnostics[[length(diagnostics) + 1L]] <- tibble::tibble(
      Analysis = analysis, Issue = object$Diagnostics$Issue, Detail = object$Diagnostics$Detail)
    Label <- function(variable) {
      if (is.na(variable)) return(NA_character_)
      value <- object$Metadata$Labels[variable]
      if (!length(value) || is.na(value) || !nzchar(value)) variable else unname(value)
    }
    Value <- function(row, variable) {
      if (is.na(variable)) return(NA_character_)
      if (!variable %in% names(df_Data)) stop("Evidence names a variable absent from ReviewData: ", variable)
      ScidrTimeObserved(df_Data, variable, match(row, df_Data$.Row), object$Metadata$TimeConversions)
    }
    covered <- character()
    Add <- function(row, method, kind, variable = NA_character_, other_row = NA_integer_,
        other_variable = NA_character_, score = NA_real_, rank = NA_real_, flag = TRUE,
        explanation = NA_character_, native = list(), conditions = list(), statistics = list(),
        reference_n = NA_real_, lower = NA_character_, upper = NA_character_, interval = NA_real_) {
      if (length(row) != 1L || is.na(row) || !row %in% df_Data$.Row ||
          (!is.na(other_row) && !other_row %in% df_Data$.Row)) stop("Evidence contains an unknown source row in ", analysis)
      if (!is.na(other_variable) && is.na(other_row)) other_row <- row
      positions <- data.frame(Row = c(row, other_row), Variable = c(variable, other_variable), stringsAsFactors = FALSE)
      positions <- positions[!is.na(positions$Row), , drop = FALSE]
      positions <- positions[order(positions$Row, positions$Variable, na.last = TRUE), , drop = FALSE]
      rownames(positions) <- NULL
      issue_id <- ScidrReviewHash(list(analysis, snapshots[[analysis]], positions))
      suggested <- switch(kind,
        date_interval = "Verify both dates and the events they represent; check encounter assignment and source records.",
        date_column = "Verify the recorded date against the source; an unusual calendar date can still be valid.",
        interval = "Check both endpoint values, elapsed time and visit assignment. A negative change can be legitimate.",
        cell = "Check valid coding, variable definition and source entry; derived changes can legitimately be negative.",
        "Inspect retained source values to identify which fields, if any, need verification.")
      if (method == "missingness") suggested <- "Check collection requirements, skip logic, study design and source records. A blank can be legitimate; do not fill it automatically."
      temporal <- object$Metadata$TimeConversions
      time_kind <- temporal$TimeType[match(variable, temporal$Variable)]
      if (method != "missingness" && length(time_kind) && !is.na(time_kind)) suggested <- if (time_kind == "clock")
        "Verify the original clock time and field definition, including its date/context. A clock alone does not establish elapsed time." else
        "Verify the recorded duration, input units and field definition. Negative or long durations are not automatically incorrect."
      if (is.na(explanation) || !nzchar(explanation)) explanation <- if (kind == "interval")
        paste(Label(variable), "has an unusual transition or elapsed-time-adjusted change; inspect both endpoints.") else
          paste(if (is.na(variable)) "Record combination" else Label(variable), "selected by", method, "for review.")
      evidence[[length(evidence) + 1L]] <<- ScidrReviewOneRow(IssueID = issue_id,
        EvidenceID = ScidrReviewHash(list(issue_id, method, kind, native, score, flag)),
        Analysis = analysis, SourceRow = as.integer(row), OtherRow = as.integer(other_row),
        Variable = variable, Label = Label(variable), OtherVariable = other_variable, OtherLabel = Label(other_variable),
        Observed = Value(row, variable), OtherObserved = if (!is.na(other_row)) Value(other_row, other_variable) else NA_character_,
        Method = method, Kind = kind, Score = as.numeric(score), Rank = as.numeric(rank), Flag = as.logical(flag),
        Explanation = as.character(explanation), SuggestedCheck = suggested, ReferenceN = as.numeric(reference_n),
        ReferenceLower = as.character(lower), ReferenceUpper = as.character(upper),
        Expected = as.character(if (!is.null(native$ExpectedDisplay)) native$ExpectedDisplay else if (!is.null(native$Expected)) native$Expected else NA_character_),
        Elapsed = as.numeric(if (!is.null(native$Elapsed)) native$Elapsed else NA_real_),
        ElapsedUnit = as.character(if (!is.null(object$Metadata$Settings[[method]]$TimeUnits)) object$Metadata$Settings[[method]]$TimeUnits else NA_character_),
        ActualIntervalDays = as.numeric(interval),
        # Evidence from a C++ backend is repaired on the way in, so a value R
        # cannot touch never reaches the fingerprint, the pages, the plots or the
        # workbook. Repairing it here also means a result object detected before
        # this fix only needs a new queue, not a new detection run.
        Conditions = list(ScidrRecoverValues(conditions)),
        GroupStatistics = list(ScidrRecoverValues(statistics)),
        NativeEvidence = list(ScidrRecoverValues(native)))
      covered <<- unique(c(covered, paste(method, c(row, other_row[!is.na(other_row)]))))
    }
    df_Records <- object$Results
    if (!all(c(".Row", "Method", "Status", "Score", "Flag") %in% names(df_Records)) ||
        any(!df_Records$.Row %in% df_Data$.Row)) stop("Invalid record results in ", analysis)
    if (anyDuplicated(df_Records[c(".Row", "Method")])) stop("Duplicate record/method results in ", analysis)
    IsEvaluated <- function(row, method) any(df_Records$.Row == row & df_Records$Method == method & df_Records$Status == "evaluated")
    for (table_name in c("DateResults", "CellResults", "IntervalResults")) {
      df_Table <- object[[table_name]]
      if (is.null(df_Table) || !nrow(df_Table)) next
      if (!all(c(".Row", "Method", "Flag") %in% names(df_Table))) stop("Incomplete evidence table: ", table_name)
      for (i in which(df_Table$Flag %in% TRUE)) {
        row <- df_Table$.Row[i]
        method <- df_Table$Method[i]
        if (!IsEvaluated(row, method) || ("Status" %in% names(df_Table) && !is.na(df_Table$Status[i]) && df_Table$Status[i] != "evaluated")) next
        variable <- as.character(ScidrReviewField(df_Table, "Variable", i, NA_character_))
        other_row <- NA_integer_
        other_variable <- NA_character_
        statistics <- ScidrReviewField(df_Table, "GroupStatistics", i, list())
        if (!is.list(statistics)) statistics <- list()
        lower <- upper <- NA_character_
        reference_n <- NA_real_
        interval <- NA_real_
        if (table_name == "DateResults") {
          pair <- identical(df_Table$CheckType[i], "interval")
          kind <- if (pair) "date_interval" else "date_column"
          if (pair) { other_variable <- df_Table$OtherVariable[i]; other_row <- row }
          reference_n <- df_Table$N[i]
          lower <- if (pair) as.character(df_Table$LowerFenceDays[i]) else df_Table$LowerDate[i]
          upper <- if (pair) as.character(df_Table$UpperFenceDays[i]) else df_Table$UpperDate[i]
          if (pair) interval <- df_Table$ValueDays[i]
        } else if (table_name == "IntervalResults") {
          kind <- "interval"
          other_row <- df_Table$PreviousRow[i]
          other_variable <- variable
          # IntervalResults stores current row in .Row, preceding row in PreviousRow.
          reference_n <- ScidrReviewField(df_Table, "ConditioningCount", i, NA_real_)
        } else {
          if (method == "conditional" && !length(statistics)) AddDiagnostic("legacy_conditional_statistics", "Structured group statistics are unavailable in older findings; original explanations are retained. Rerun detection only if structured references are needed.")
          kind <- if (ScidrAnomalyIsDate(df_Data[[variable]])) "date_column" else "cell"
          reference_n <- if (!is.null(statistics$n_obs)) statistics$n_obs else NA_real_
          FormatBound <- function(value) {
            if (is.null(value)) return(NA_character_)
            conversion <- object$Metadata$DateConversions[[variable]]
            declaration <- object$Metadata$TimeConversions
            kind <- declaration$TimeType[match(variable, declaration$Variable)]
            if (!is.null(conversion)) ScidrAnomalyFormatDate(value, conversion) else if (length(kind) && !is.na(kind) && kind == "duration") ScidrTimeFormat(value, "duration") else as.character(value)
          }
          lower <- FormatBound(statistics$lower_thr)
          upper <- FormatBound(statistics$upper_thr)
        }
        Add(row, method, kind, variable, other_row, other_variable,
          score = ScidrReviewField(df_Table, "Score", i, NA_real_),
          explanation = ScidrReviewField(df_Table, "Explanation", i, NA_character_),
          native = as.list(df_Table[i, , drop = FALSE]),
          conditions = ScidrReviewField(df_Table, "Conditions", i, list()), statistics = statistics,
          reference_n = reference_n, lower = lower, upper = upper, interval = interval)
      }
    }
    eligible <- which(df_Records$Status == "evaluated" & df_Records$Flag %in% TRUE)
    if (include_ranked) {
      # A default percentile cutoff must not shrink the candidate list below the top_n
      # highest scores a reviewer would otherwise see, so those methods still top up.
      # Methods with a native or user cutoff judged their unflagged records; leave them.
      budgeted <- setdiff(object$Metadata$FlagPercentileMethods, names(object$Metadata$Thresholds))
      unclassified <- which(df_Records$Status == "evaluated" & is.finite(df_Records$Score) &
        !df_Records$Flag %in% TRUE & (is.na(df_Records$Flag) | df_Records$Method %in% budgeted))
      groups <- interaction(df_Records$Method[unclassified],
        if ("Stratum" %in% names(df_Records)) df_Records$Stratum[unclassified] else rep(1, length(unclassified)), drop = TRUE)
      for (rows in split(unclassified, groups)) eligible <- c(eligible, head(rows[order(-df_Records$Score[rows])], top_n))
    }
    for (i in unique(eligible)) {
      row <- df_Records$.Row[i]
      method <- df_Records$Method[i]
      if (paste(method, row) %in% covered) next
      Add(row, method, if (is.na(df_Records$Flag[i])) "ranked" else "record",
        score = df_Records$Score[i], rank = ScidrReviewField(df_Records, "Rank", i, NA_real_),
        flag = df_Records$Flag[i], explanation = ScidrReviewField(df_Records, "Explanation", i, NA_character_),
        native = as.list(df_Records[i, , drop = FALSE]))
    }
  }
  df_Evidence <- dplyr::bind_rows(evidence)
  df_Evidence <- df_Evidence[!duplicated(df_Evidence$EvidenceID), , drop = FALSE]
  issues <- list(ScidrReviewEmptyIssues())
  for (issue_id in unique(df_Evidence$IssueID)) {
    df_Current <- df_Evidence[df_Evidence$IssueID == issue_id, , drop = FALSE]
    first <- df_Current[1, , drop = FALSE]
    analysis <- first$Analysis
    data <- sources[[analysis]]$ReviewData
    row_index <- match(first$SourceRow, data$.Row)
    other_index <- match(first$OtherRow, data$.Row)
    selected_id <- ids[[analysis]]
    kinds <- unique(df_Current$Kind)
    # Provisional; ScidrReviewOrder() recomputes it from the strongest evidence.
    category <- ScidrAnomalyType(df_Current$Method[1], if ("date_interval" %in% kinds) "date_interval" else kinds[1],
      !is.na(first$Variable) && nzchar(first$Variable),
      length(first$Conditions[[1]]) > 0L)
    vars_Context <- contexts[[analysis]]
    context <- if (length(vars_Context)) paste(vapply(vars_Context, function(variable) {
      label <- sources[[analysis]]$Metadata$Labels[variable]
      if (!length(label) || is.na(label) || !nzchar(label)) label <- variable
      paste0(label, ": ", ScidrReviewFormat(data[[variable]][row_index]),
        if (!is.na(other_index) && other_index != row_index) paste0("; other row: ", ScidrReviewFormat(data[[variable]][other_index])) else "")
    }, character(1)), collapse = " | ") else NA_character_
    issues[[length(issues) + 1L]] <- ScidrReviewOneRow(IssueID = issue_id, Analysis = analysis,
      SourceRow = first$SourceRow, OtherRow = first$OtherRow,
      StudyID = if (length(selected_id)) ScidrReviewFormat(data[[selected_id]][row_index]) else NA_character_,
      OtherStudyID = if (length(selected_id) && !is.na(other_index)) ScidrReviewFormat(data[[selected_id]][other_index]) else NA_character_,
      Context = context, Variable = first$Variable, Label = if (is.na(first$Variable)) "Whole record" else first$Label,
      OtherVariable = first$OtherVariable, OtherLabel = first$OtherLabel,
      Observed = first$Observed, OtherObserved = first$OtherObserved, ReviewCategory = category,
      Explanation = paste(unique(stats::na.omit(df_Current$Explanation)), collapse = "\n"),
      SuggestedCheck = paste(unique(df_Current$SuggestedCheck), collapse = "\n"), EvidenceCount = nrow(df_Current),
      ReviewStatus = "Unreviewed", Priority = "Unassigned", ProposedCorrection = NA_character_,
      SourceChecked = NA_character_, Reviewer = NA_character_, ReviewedOn = NA_character_, Notes = NA_character_)
  }
  df_Issues <- dplyr::bind_rows(issues)
  front <- c("IssueID", "Analysis", "StudyID", "SourceRow", "ReviewCategory", "ReviewStatus", "Priority",
    "Label", "Observed", "OtherLabel", "OtherObserved", "Explanation", "SuggestedCheck",
    "ProposedCorrection", "SourceChecked", "Reviewer", "ReviewedOn", "Notes")
  df_Issues <- df_Issues[c(front, setdiff(names(df_Issues), front))]
  df_Issues <- df_Issues[order(match(df_Issues$Analysis, names(sources)), df_Issues$SourceRow,
    ifelse(is.na(df_Issues$OtherRow), -1, df_Issues$OtherRow), df_Issues$IssueID), , drop = FALSE]
  queue <- structure(list(Issues = df_Issues, Evidence = df_Evidence, Summary = ScidrReviewSummarize(sources),
    Diagnostics = dplyr::distinct(dplyr::bind_rows(c(list(tibble::tibble(Analysis = character(), Issue = character(), Detail = character())), diagnostics))),
    Metadata = list(FormatVersion = "1", Snapshots = unlist(snapshots), IDVariables = ids, ContextVariables = contexts,
      IncludeRanked = include_ranked, TopN = top_n,
      Notice = "Flags are review candidates, not an estimated error rate. Conditional detection currently returns at most one finding per record.",
      QueueFingerprint = ScidrReviewHash(list(df_Issues[setdiff(names(df_Issues), ScidrReviewFields())], df_Evidence, snapshots))),
    Sources = sources), class = c("AnomalyReviewQueue", "list"))
  queue <- ScidrReviewSubsetVariables(queue, variables)
  ScidrReviewCodebook(queue, codebook)
}

# Narrow a queue to the variables a reader cares about, without refitting. The
# models stay as they were fitted on the full selection, so a kept finding is
# still conditioned on every screened column; only what is listed changes.
# Findings with no attributed field concern the whole record and are kept: they
# were not raised about any one variable, so no variable list excludes them.
ScidrReviewSubsetVariables <- function(queue, variables) {
  if (is.null(variables)) return(queue)
  if (!is.character(variables) || !length(variables) || anyNA(variables) || any(!nzchar(variables)) || anyDuplicated(variables))
    stop("variables must be NULL or unique column names.")
  known <- unique(unlist(lapply(queue$Sources, function(source) source$Metadata$Variables)))
  unknown <- setdiff(variables, known)
  if (length(unknown)) stop("variables names columns that were not screened: ", paste(unknown, collapse = ", "))
  df_Evidence <- queue$Evidence
  keep <- is.na(df_Evidence$Variable) | df_Evidence$Variable %in% variables |
    (!is.na(df_Evidence$OtherVariable) & df_Evidence$OtherVariable %in% variables)
  queue$Evidence <- df_Evidence[keep, , drop = FALSE]
  df_Issues <- queue$Issues[queue$Issues$IssueID %in% queue$Evidence$IssueID, , drop = FALSE]
  df_Issues$EvidenceCount <- as.integer(table(factor(queue$Evidence$IssueID, levels = df_Issues$IssueID))[df_Issues$IssueID])
  queue$Issues <- df_Issues
  queue$Metadata$Variables <- variables
  queue$Metadata$QueueFingerprint <- ScidrReviewHash(list(df_Issues[setdiff(names(df_Issues), ScidrReviewFields())],
    queue$Evidence, as.list(queue$Metadata$Snapshots)))
  queue
}
