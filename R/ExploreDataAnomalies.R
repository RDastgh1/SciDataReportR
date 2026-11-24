#' Explore data anomaly results interactively
#'
#' Display coverage, scores, flags, explanations, and longitudinal histories.
#' @param AnomalyObj A result from [DetectDataAnomalies()] or
#'   [DetectLongitudinalAnomalies()], or an AnomalyReviewQueue from
#'   [CreateAnomalyReviewQueue()]. Detection results are converted automatically;
#'   the issue worklist appears before the original technical views.
#' @param codebook Optional presentation codebook as in [DetectDataAnomalies()].
#'   Overrides matching saved labels/categories for this view only, including
#'   older results. Does not refit models or change decisions and fingerprints.
#' @param Title Dashboard title.
#' @param TopN Positive integer retained for compatibility with earlier history
#'   previews. Interactive plots are selected individually; all rows are accessible.
#' @param TableHeight Positive scrolling table height in pixels, default 500. NULL
#'   also uses 500 pixels for compatibility.
#' @param Detail `Compact` for collapsed detail sections, or `Full` to open them.
#' @param IncludePlots Logical; enable interactive distributions, date intervals and
#'   histories, drawn only when an issue is selected. Requires optional plotly;
#'   default TRUE selects and plots the first issue automatically when plotly is
#'   installed. Set FALSE for smaller HTML. Use
#'   [PlotDataAnomalyTrajectories()] separately for selected histories.
#' @details Codebook category is distinct from issue type. Field findings retain
#'   categories of affected variables; isolation records use the five largest
#'   positive screening effects or saved SHAP contributions as suggested inspection
#'   categories. Without positive explanations the category is unspecified.
#'   Variable-level evidence supports findings and is not independent confirmation.
#' @details The first review table groups findings by analysis, stratum, and
#'   original variable. It orders variables by their strongest supported evidence,
#'   then distinct affected-record count and variable name. Selecting a variable
#'   filters the individual findings; Show all findings resets that selection.
#'   Both tables carry a dropdown per categorical column, covering extreme
#'   value, anomaly type, severity tier, codebook category, method and analysis
#'   group; a column with fewer than two distinct values is omitted.
#'   Whole-record combinations share the findings table and do not count
#'   isolation screening leads as affected fields; filter by anomaly type to
#'   list them on their own. Lightweight, lazily rendered distribution
#'   previews highlight review values against retained within-stratum data;
#'   conditional subgroup evidence remains in selected details. Missingness uses
#'   observed/missing bars and clock times use circular previews. Short reasons
#'   are shared with Excel; full explanations remain in comparison details.
#'   Review decisions are shown in selected details, not the main or printed queue.
#' @details Every finding carries one anomaly type from a
#'   fixed list, shown as a sortable column on both tables with a matching
#'   filter, and one severity tier: the rank of that kind of finding, also a
#'   sortable column, which drives the suggested review order and is defined on
#'   the page. A value at least ten times further from its column median than
#'   any other record is marked as an extreme value and listed ahead of every
#'   tier. An anomaly map above the tables draws one cell per modelled
#'   variable and record, grey where nothing was reported and coloured by the
#'   kind of anomaly, with the record, value and reason on hover. A cell takes
#'   its most serious finding; findings with no attributed field, and fields
#'   outside the modelled set, cannot be placed and are counted beneath it. Developer-facing views are deliberately absent from the page:
#'   per-record scores, per-cell evidence and screening diagnostics remain on
#'   the result and queue objects and in the Excel workbook.
#' @details Row details are generated in the browser only when selected. The
#'   record grid contains scores and context; retained source values are stored
#'   once and looked up by original row key in a separate evidence panel. Grid searches
#'   search grid columns, not source values inside unexpanded details. Full opens
#'   section accordions but does not eagerly display every record. Tables use
#'   deferred rendering and virtual scrolling; TopN never limits table access.
#' @param variables Optional character vector of screened column names. NULL
#'   shows every finding. Otherwise the page lists only findings attributed to a
#'   listed column, together with whole-record findings, which are raised about
#'   no one field. This narrows the display only: detection is unchanged, the
#'   models were fitted on the full selection, so a shown finding is still
#'   conditioned on every screened column, and the coverage summary still
#'   describes the whole screen.
#' @return An `htmltools::tagList()` for Quarto, RStudio Viewer, or standalone HTML.
#'   No files are saved. Requires DT and htmltools; interactive plots use
#'   the locally bundled Plotly JavaScript, with no external services. Requires digest for queue identity. Scores retain their native scales and are not error probabilities.
#' @family data anomalies
#' @seealso [SaveDataAnomalyExplorer()], [CreateDataAnomalyReport()]
#' @export
ExploreDataAnomalies <- function(AnomalyObj, Title = "Data anomaly explorer",
    TopN = 10, TableHeight = 500, Detail = c("Compact", "Full"), IncludePlots = TRUE, codebook = NULL,
    variables = NULL) {
  codebook <- ScidrAnomalyCodebook(codebook)
  # DT advises server-side processing above 1.5 MB per table; a standalone,
  # emailable file cannot use it, so the advice is not actionable here.
  previous_options <- options(DT.warn.size = FALSE)
  on.exit(options(previous_options), add = TRUE)
  Detail <- match.arg(Detail)
  ReviewObj <- if (inherits(AnomalyObj, "AnomalyReviewQueue")) AnomalyObj else CreateAnomalyReviewQueue(AnomalyObj)
  ScidrReviewValidate(ReviewObj)
  ReviewObj <- ScidrReviewSubsetVariables(ReviewObj, variables)
  ReviewObj <- ScidrReviewCodebook(ReviewObj, codebook)
  sources <- stats::setNames(lapply(names(ReviewObj$Sources), function(analysis) ScidrReviewDisplaySource(ReviewObj, analysis)), names(ReviewObj$Sources))
  AnomalyObj <- sources[[1]]
  if (!is.character(Title) || length(Title) != 1L || is.na(Title)) stop("Title must be one string.")
  if (!is.numeric(TopN) || length(TopN) != 1L || !is.finite(TopN) || TopN < 1 || TopN != floor(TopN)) stop("TopN must be a positive integer.")
  if (!is.null(TableHeight) && (!is.numeric(TableHeight) || length(TableHeight) != 1L || !is.finite(TableHeight) || TableHeight <= 0)) stop("TableHeight must be positive.")
  if (!is.logical(IncludePlots) || length(IncludePlots) != 1L || is.na(IncludePlots)) stop("IncludePlots must be TRUE or FALSE.")
  for (package in c("DT", "htmltools")) if (!requireNamespace(package, quietly = TRUE)) stop("Install package '", package, "' to use the explorer.")
  plot_dependencies <- NULL
  if (IncludePlots) {
    # Render directly with the optional package's bundled JS, without loading R plotting engines.
    fn_Plotly <- system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly")
    if (!nzchar(fn_Plotly)) {
      message("Install package 'plotly' to enable interactive anomaly figures.")
      IncludePlots <- FALSE
    }
    if (IncludePlots)
    plot_dependencies <- htmltools::htmlDependency("scidr-plotly", as.character(utils::packageVersion("plotly")),
      src = dirname(fn_Plotly), script = basename(fn_Plotly), all_files = FALSE)
  }
  display <- ScidrReviewNarratives(ReviewObj)
  variable_review <- ScidrReviewVariableSummary(ReviewObj, display)
  df_Coverage <- ReviewObj$Summary
  df_Coverage$ScoredOnly <- vapply(seq_len(nrow(df_Coverage)), function(i) {
    records <- sources[[df_Coverage$Analysis[i]]]$Results
    sum(records$Method == df_Coverage$Method[i] & records$Status == "evaluated" & is.na(records$Flag) & is.finite(records$Score))
  }, integer(1))
  df_Coverage$CandidatesShown <- vapply(seq_len(nrow(df_Coverage)), function(i)
    length(unique(ReviewObj$Evidence$IssueID[ReviewObj$Evidence$Analysis == df_Coverage$Analysis[i] & ReviewObj$Evidence$Method == df_Coverage$Method[i]])), integer(1))
  df_Coverage$Classification <- vapply(seq_len(nrow(df_Coverage)), function(i) {
    object <- sources[[df_Coverage$Analysis[i]]]
    method <- df_Coverage$Method[i]
    scored <- object$Results[object$Results$Method == method & object$Results$Status == "evaluated", ]
    if (nrow(scored) && all(is.na(scored$Flag))) "No cutoff supplied; scores only"
    else if (method %in% object$Metadata$FlagPercentileMethods && !method %in% names(object$Metadata$Thresholds))
      paste0("Default top ", signif(100 * (1 - object$Metadata$FlagPercentile), 3), "% cutoff")
    else "Native or explicit cutoff flags"
  }, character(1))
  unclassified <- df_Coverage$Classification == "No cutoff supplied; scores only"
  df_Coverage$Flagged <- as.character(df_Coverage$Flagged)
  df_Coverage$Flagged[unclassified] <- "Not classified"
  df_Coverage$FlagFraction[unclassified] <- NA_real_
  df_Coverage$FlagFraction <- NULL
  if (length(sources) == 1L) df_Coverage$Analysis <- NULL
  names(df_Coverage)[names(df_Coverage) == "ScoredOnly"] <- "Scored only"
  names(df_Coverage)[names(df_Coverage) == "CandidatesShown"] <- "Candidates shown"
  df_Methods <- ScidrAnomalyMethods()
  method_names <- stats::setNames(df_Methods$Method, df_Methods$Key)
  requested <- unique(ReviewObj$Summary$Method)
  df_Methods <- df_Methods[order(match(df_Methods$Key, requested), match(df_Methods$Key, ScidrAnomalyMethods()$Key), na.last = TRUE), ]
  df_Methods$Method <- paste0(df_Methods$Method, " (", df_Methods$Key, ")")
  df_Methods$AppliesTo <- ifelse(df_Methods$Key %in% c("change", "trajectory", "transition"), "Longitudinal", "Cross-sectional")
  df_Coverage$Method <- unname(method_names[df_Coverage$Method])
  display_names <- c(SuggestedReviewOrder = "Suggested review order", ReviewTier = "Severity tier", ExtremeValue = "Extreme value", OrderRationale = "Order rationale", FindingType = "Anomaly type", MeasureValue = "Measure and recorded value", .Row = "Source row", SourceRow = "Source row", OtherRow = "Other source row",
    StudyID = "Study ID", OtherStudyID = "Other study ID", Variable = "Measure", OtherVariable = "Compared measure",
    Observed = "Recorded value", OtherObserved = "Compared value", ReviewCategory = "Anomaly type", CodebookCategory = "Codebook category",
    CategoryBasis = "Category assignment", CategoryVariables = "Columns supplying categories", ConditionCategories = "Comparison-variable categories",
    ReviewStatus = "Review decision", .ReviewStatus = "Review status", WhyFlagged = "Why flagged",
    ComparisonDetails = "Comparison details", Label = "Measure label", OtherLabel = "Compared measure label",
    EvaluatedRecords = "Evaluated records", FlaggedRecords = "Records flagged", MethodsFlagged = "Methods flagging",
    FlagFraction = "Fraction flagged", NotApplicable = "Not applicable", Stratum = "Analysis group", Explanation = "Original explanation", SuggestedCheck = "Suggested check",
    StratumLabel = "Analysis group", AffectedRecords = "Affected records", StrongestFinding = "Strongest finding", CommonReasons = "Common reasons", MiniDistribution = "Column distribution",
    ReferenceN = "Reference observations", ReferenceLower = "Reference lower boundary", ReferenceUpper = "Reference upper boundary")
  # The identifier column is the study's own, so it is headed with that column's
  # label -- or its name, where no label was supplied -- rather than a generic
  # "Study ID" that does not appear anywhere in the data. Analyses that used
  # different identifier columns have no single heading to give, so those keep the
  # generic one. The workbook's headers are deliberately left fixed: decisions are
  # read back in by header, so those must not move with the data.
  id_label <- unique(unlist(lapply(sources, function(source) {
    id <- source$Metadata$IDVariable
    if (!length(id) || is.na(id)) return(NULL)
    label <- unname(source$Metadata$Labels[id])
    if (length(label) && !is.na(label) && nzchar(label)) label else id
  })))
  if (length(id_label) == 1L)
    display_names[c("StudyID", "OtherStudyID")] <- c(id_label, paste("Other", id_label))
  # The print view has no column chooser, so it picks its own columns: the same
  # ones the interactive table shows, minus any that are constant here.
  PrintVariableReview <- function(df_Table) {
    columns <- c(if (length(sources) > 1L) "Analysis",
      if (length(unique(df_Table$Stratum)) > 1L) "Stratum",
      "Label", "Variable", "ExtremeValue", "ReviewTier", "ReviewCategory",
      if (length(unique(df_Table$CodebookCategory)) > 1L) "CodebookCategory",
      "AffectedRecords", "CommonReasons")
    columns <- intersect(columns, names(df_Table))
    df_Print <- as.data.frame(df_Table[columns])
    names(df_Print) <- vapply(columns, function(name)
      if (name %in% names(display_names)) display_names[[name]] else name, character(1))
    df_Print
  }
  PrintTable <- function(df_Data, queue = FALSE) {
    Escape <- function(x) htmltools::htmlEscape(ifelse(is.na(x), "", as.character(x)))
    cells <- lapply(df_Data, function(x) paste0("<td>", Escape(x), "</td>"))
    rows <- if (nrow(df_Data)) do.call(paste0, cells) else character()
    htmltools::HTML(paste0('<table data-quarto-disable-processing="true" class="',
      if (queue) "scidr-static-table scidr-print-issues" else "scidr-static-table", '"><thead><tr>',
      paste0("<th>", Escape(names(df_Data)), "</th>", collapse = ""), "</tr></thead><tbody>",
      paste0("<tr>", rows, "</tr>", collapse = ""), "</tbody></table>"))
  }
  df_Print <- if (!nrow(display$Issues)) data.frame(Record=character(), Measure=character(), Values=character(), Reason=character(), CodebookCategory=character()) else data.frame(
    SuggestedOrder = display$Issues$SuggestedReviewOrder,
    Record = paste0("Source row ", display$Issues$SourceRow, ifelse(is.na(display$Issues$StudyID), "", paste0("; ID ", display$Issues$StudyID))),
    Measure = ifelse(is.na(display$Issues$Variable), "Record combination (no attributed field)", paste0(display$Issues$Label, " (", display$Issues$Variable, ")",
      ifelse(is.na(display$Issues$OtherVariable), "", paste0("; compared with ", display$Issues$OtherLabel, " (", display$Issues$OtherVariable, ")")))),
    Values = ifelse(is.na(display$Issues$Variable), "See retained record",
      paste0(vapply(display$Issues$Observed, ScidrReviewRoundedValue, character(1), USE.NAMES = FALSE),
        ifelse(is.na(display$Issues$OtherObserved), "",
          paste0("; compared value: ", vapply(display$Issues$OtherObserved, ScidrReviewRoundedValue, character(1), USE.NAMES = FALSE))))),
    Reason = display$Issues$WhyFlagged,
    CodebookCategory = display$Issues$CodebookCategory)
  names(df_Print)[names(df_Print) == "CodebookCategory"] <- "Codebook category"
  names(df_Print)[names(df_Print) == "SuggestedOrder"] <- "Suggested review order"
  plot_sources_emitted <- FALSE
  Table <- function(df_Table, kind = "evidence", source = AnomalyObj, issue_evidence = NULL) {
    if (is.null(df_Table) || !nrow(df_Table)) return(htmltools::tags$p("No results available."))
    df_Table <- ScidrAnomalyDisplayFields(as.data.frame(df_Table), source$Metadata$VariableMetadata)
    if (kind != "queue" && "Variable" %in% names(df_Table) && !is.null(source$Metadata$VariableMetadata)) {
      catalog <- source$Metadata$VariableMetadata
      df_Table$CodebookCategory <- vapply(seq_len(nrow(df_Table)), function(i) {
        fields <- unique(c(df_Table$Variable[i], df_Table$OtherVariable[i]))
        categories <- catalog$Category[match(intersect(fields, catalog$Variable), catalog$Variable)]
        if (length(categories)) paste(unique(categories), collapse = "; ") else "Uncategorized"
      }, character(1))
    }
    if (kind == "summary" && "Flagged" %in% names(df_Table) && !is.null(source$Results)) {
      df_Table$Flagged <- as.character(df_Table$Flagged)
      for (i in seq_len(nrow(df_Table))) {
        records <- source$Results[source$Results$Method == df_Table$Method[i] & source$Results$Status == "evaluated", ]
        if (nrow(records) && all(is.na(records$Flag))) df_Table$Flagged[i] <- "Not classified; scores only"
      }
    }
    if (kind == "queue") {
      df_Table$MeasureValue <- vapply(seq_len(nrow(df_Table)), function(i) {
        if (is.na(df_Table$Variable[i])) return(paste("Record combination; no attributed field.", df_Table$Methods[i]))
        paste0(df_Table$Label[i], " (", df_Table$Variable[i], "): ", ScidrReviewRoundedValue(df_Table$Observed[i]),
          if (!is.na(df_Table$OtherVariable[i])) paste0("; ", df_Table$OtherLabel[i], " (", df_Table$OtherVariable[i], "): ", ScidrReviewRoundedValue(df_Table$OtherObserved[i])) else "",
          "\n", df_Table$Methods[i])
      }, character(1))
    }

    preferred <- c(if (kind == "variables") c("SuggestedReviewOrder", "Analysis", "StratumLabel", "Variable", "Label", "ExtremeValue", "ReviewTier", "ReviewCategory", "CodebookCategory", "AffectedRecords", "CommonReasons", "StrongestFinding", "MiniDistribution"), "SuggestedReviewOrder", "StudyID", "SourceRow", "Analysis", "MeasureValue", ".Row", "Method", "Variable", "OtherVariable", "Observed", "OtherObserved", "WhyFlagged", "ReviewCategory", "ReviewStatus", "Priority")
    df_Table <- df_Table[unique(c(intersect(preferred, names(df_Table)), setdiff(names(df_Table), preferred)))]
    if ("Flag" %in% names(df_Table)) df_Table <- df_Table[order(!df_Table$Flag, na.last = TRUE), , drop = FALSE]
    # Keep structured technical evidence as JSON cells; format only when displayed.
    df_Table[] <- lapply(df_Table, function(x) if (is.list(x)) I(x) else ScidrReviewFormat(x))
    if ("Flag" %in% names(df_Table)) df_Table$.ReviewStatus <- ifelse(is.na(df_Table$Flag), "Not classified", ifelse(df_Table$Flag == "TRUE", "Review", "Not flagged"))
    id_var <- source$Metadata$IDVariable
    visible <- switch(kind,
      queue = c("SuggestedReviewOrder", if (any(!is.na(df_Table$StudyID) & nzchar(df_Table$StudyID))) "StudyID", "SourceRow",
        if (length(unique(df_Table$Analysis)) > 1L) "Analysis", "ExtremeValue", "ReviewTier", "ReviewCategory", "CodebookCategory", "MeasureValue", "WhyFlagged"),
      variables = c("SuggestedReviewOrder", if (length(sources) > 1L) "Analysis", if (length(unique(df_Table$StratumLabel)) > 1L) "StratumLabel", "Variable", "ExtremeValue", "ReviewTier", "ReviewCategory", "CodebookCategory", "AffectedRecords", "CommonReasons", "MiniDistribution"),
      summary = names(df_Table),
      records = c(id_var, "Method", "Status", ".ReviewStatus", "Rank", "Score", "Explanation"),
      dates = c(id_var, "CodebookCategory", "Variable", "OtherVariable", "Observed", "OtherObserved", ".ReviewStatus", "Explanation"),
      cells = c(id_var, "CodebookCategory", "Variable", "Observed", if ("ExpectedDisplay" %in% names(df_Table)) "ExpectedDisplay" else "Expected", ".ReviewStatus", "Score"),
      names(df_Table))
    if (!length(id_var)) visible <- c(".Row", visible)
    empty_compared <- intersect(c("OtherVariable", "OtherObserved"), names(df_Table))
    empty_compared <- empty_compared[vapply(df_Table[empty_compared], function(x) all(is.na(x) | x == ""), logical(1))]
    visible <- setdiff(visible, empty_compared)
    labels <- utils::modifyList(as.list(display_names), as.list(source$Metadata$Labels))
    meta <- list(methodNames = as.list(method_names), multipleAnalyses = length(sources) > 1L, columns = names(df_Table), labels = labels, evidence = list())
    if (kind == "variables") meta$previews <- variable_review$Previews
    if (!is.null(issue_evidence)) issue_evidence <- issue_evidence[issue_evidence$IssueID %in% df_Table$IssueID, , drop = FALSE]
    if (!is.null(issue_evidence) && nrow(issue_evidence)) {
      meta$evidence <- lapply(split(issue_evidence, issue_evidence$IssueID), function(df_Evidence) {
        lapply(seq_len(nrow(df_Evidence)), function(i) as.list(ScidrReviewFlat(df_Evidence[i, , drop = FALSE])))
      })
      meta$explanations <- ScidrReviewExplanationMeta(ReviewObj)
      meta$explanations <- meta$explanations[intersect(names(meta$explanations), df_Table$IssueID)]
      meta$sourceLabels <- lapply(sources, function(object) as.list(object$Metadata$Labels))
      if (IncludePlots) {
        if (!plot_sources_emitted) {
          meta$plotSources <- ScidrReviewPlotSources(ReviewObj)
          plot_sources_emitted <<- TRUE
        } else meta$useSharedPlotSources <- TRUE
        meta$plotConditions <- lapply(split(issue_evidence, issue_evidence$IssueID), function(df_Evidence) {
          conditional <- which(df_Evidence$Method == "conditional")
          if (!length(conditional)) return(list())
          i <- conditional[1]
          conditions <- df_Evidence$Conditions[[i]]
          object <- sources[[df_Evidence$Analysis[i]]]
          lapply(conditions, function(condition) {
            conversion <- object$Metadata$DateConversions[[condition$column]]
            if (!is.null(conversion) && is.numeric(condition$value_comp)) {
              condition$value_comp <- as.numeric(as.POSIXct(ScidrAnomalyFormatDate(condition$value_comp, conversion), tz = "UTC")) / 86400
            }
            condition
          })
        })
        # The learned boundary the finding is stated against, so the panel can draw
        # it beside the record's own value instead of leaving the reader to guess
        # where "above the subgroup boundary" sits.
        meta$plotBounds <- lapply(split(issue_evidence, issue_evidence$IssueID), function(df_Evidence) {
          i <- which(df_Evidence$Method == "conditional")
          if (!length(i)) return(NULL)
          Bound <- function(x) {
            x <- suppressWarnings(as.numeric(x[[1]]))
            if (length(x) == 1L && is.finite(x)) x else NULL
          }
          bounds <- list(lower = Bound(df_Evidence$ReferenceLower[i[1]]), upper = Bound(df_Evidence$ReferenceUpper[i[1]]))
          if (is.null(bounds$lower) && is.null(bounds$upper)) NULL else bounds
        })
      }
    }
    render <- DT::JS("function(data, type, row, meta) {
      const esc = function(x) { return String(x == null ? '' : x).replace(/[&<>\"']/g, function(c) {return {'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;', \"'\":'&#39;'}[c];}); };
      const opts = meta.settings.oInit.scidrMeta, name = opts.columns[meta.col];
      let text = data == null ? '' : (typeof data === 'object' ? JSON.stringify(data) : String(data)), label = '';
      if (name === 'MiniDistribution' && type === 'display') {
        const preview = opts.previews[text];
        if (!opts.previewCache) opts.previewCache = {};
        if (!opts.previewCache[text]) opts.previewCache[text] = opts.renderMini(text, opts.previews);
        return opts.previewCache[text];
      }
      if (name === 'MeasureValue' && type === 'display') {
        const value = function(key) {const i=opts.columns.indexOf(key);return i < 0 || row[i] == null ? '' : String(row[i]);};
        const variable=value('Variable');
        let html = variable ? esc(value('Label') || variable) + '<br><code>' + esc(variable) + '</code> = ' + esc(value('Observed') || (value('Methods').indexOf('missingness') >= 0 ? '(missing)' : '')) : 'Record combination; no attributed field';
        if(value('OtherVariable')) html += '<br>' + esc(value('OtherLabel')) + ' <code>' + esc(value('OtherVariable')) + '</code> = ' + esc(value('OtherObserved'));
        html += '<br><small>' + esc(value('Methods').split(', ').map(function(key) {return opts.methodNames[key] || key;}).join(', ')) + '</small>';
        return '<div class=\"scidr-cell scidr-measure\">' + html + '</div>';
      }
      if (name === 'Variable' || name === 'OtherVariable') {
        const index = opts.columns.indexOf(name === 'Variable' ? 'Label' : 'OtherLabel');
        label = index >= 0 && row[index] ? row[index] : text;
        if (type !== 'display') return label + ' ' + text;
        return '<div class=\"scidr-cell\"><span class=\"scidr-field-label\">' + esc(label) + '</span><code class=\"scidr-field-name\" title=\"' + esc(text) + '\">' + esc(text) + '</code></div>';
      }
      if (type === 'display' && ['WhyFlagged','Explanation','Detail'].indexOf(name) >= 0) return '<div class=\"scidr-cell scidr-reason\">' + esc(text) + '</div>';
      if (name === 'Method' || name === 'Methods') text = text.split(', ').map(function(key) {return opts.methodNames[key] || key;}).join(', ');
      if (name === 'FlagFraction' && text !== '' && !isNaN(Number(text))) text = Number((100 * Number(text)).toPrecision(3)) + '%';
      return type === 'display' ? '<div class=\"scidr-cell\">' + esc(text) + '</div>' : text;
    }")
    if (kind == "variables") meta$renderMini <- DT::JS(ScidrReviewMiniJS())
    definitions <- list(list(targets = seq_along(df_Table) - 1L, render = render),
      list(targets = which(!names(df_Table) %in% visible) - 1L, visible = FALSE),
      list(targets = which(names(df_Table) %in% c("WhyFlagged", "Explanation", "Detail")) - 1L, width = "320px"))
    headers <- names(df_Table)
    matched <- headers %in% names(display_names)
    headers[matched] <- unname(display_names[headers[matched]])
    headers <- htmltools::htmlEscape(headers)
    callback <- DT::JS(paste0("const renderDetails = ", ScidrAnomalyDetailJS(), ";
      const renderContributions = ", ScidrAnomalyContributionsJS(), ";
      const renderPlot = ", ScidrReviewPlotJS(), ";
      const root = $(table.table().container()).closest('.scidr-anomaly-explorer');
      const sharedMeta=table.settings()[0].oInit.scidrMeta;
      if(sharedMeta.plotSources) root.data('scidrPlotSources',sharedMeta.plotSources);
      else if(sharedMeta.useSharedPlotSources) sharedMeta.plotSources=root.data('scidrPlotSources');
      // One dropdown per named column: the distinct values become the choices, a
      // column with fewer than two of them is skipped, and 'blank' labels the
      // empty value where an empty cell is a meaningful answer rather than
      // missing information.
      const addFilters = function(meta, specs) {
        const filters = $('<div class=\"scidr-queue-filters\"></div>');
        let added = 0;
        specs.forEach(function(spec) {
          const index = meta.columns.indexOf(spec.name); if(index < 0) return;
          // An empty cell becomes spec.blank where emptiness is itself an answer
          // (no extreme value), which keeps the empty string free to mean All.
          const members = function(x) {
            const raw = spec.members ? spec.members(x) : [x == null ? '' : String(x)];
            return raw.map(function(v) {return v === '' && spec.blank ? spec.blank : v;});
          };
          const choices = new Set();
          table.column(index).data().toArray().forEach(function(x) {
            members(x).forEach(function(v) {if(v !== '') choices.add(v);});
          });
          if(choices.size < 2) return;
          const control = $('<select></select>').attr('aria-label', spec.title + ' filter');
          control.append($('<option></option>').val('').text('All'));
          Array.from(choices).sort().forEach(function(v) {
            control.append($('<option></option>').val(v).text(spec.label ? spec.label(v) : v));
          });
          let selected = '';
          $.fn.dataTable.ext.search.push(function(settings,data,dataIndex) {
            if(settings.nTable !== table.table().node() || !selected) return true;
            return members(table.row(dataIndex).data()[index]).indexOf(selected) >= 0;
          });
          control.on('change',function() {selected=this.value;table.draw();});
          filters.append($('<label></label>').text(spec.title + ' ').append(control));
          added++;
        });
        if(added) $(table.table().container()).before(filters);
      };
      const variableColumn = table.settings()[0].oInit.scidrMeta.columns.indexOf('MiniDistribution');
      if (variableColumn >= 0) {
        $(table.table().body()).on('click keydown', 'tr', function(event) {
          if(event.type === 'keydown' && event.key !== 'Enter') return;
          const values = table.row(this).data(); if (!values) return;
          $(table.rows().nodes()).removeClass('selected'); $(this).addClass('selected');
          const meta=table.settings()[0].oInit.scidrMeta;
          root.trigger('scidrVariable', [String(values[meta.columns.indexOf('Variable')] || '')]);
        });
        root.find('.scidr-variable-reset').on('click', function() {$(table.rows().nodes()).removeClass('selected');root.trigger('scidrVariable',['']);});
        // Every categorical column on this table gets its own dropdown, so a
        // reader can narrow to one anomaly type, one codebook category or one
        // analysis group without typing into a column search box.
        addFilters(sharedMeta, [
          {name: 'ExtremeValue', title: 'Extreme value', blank: 'No'},
          {name: 'ReviewCategory', title: 'Anomaly type'},
          {name: 'ReviewTier', title: 'Severity tier'},
          {name: 'CodebookCategory', title: 'Codebook category'},
          {name: 'StratumLabel', title: 'Analysis group'},
          {name: 'Analysis', title: 'Analysis'}]);
        return;
      }
      const panel = $('<section class=\"scidr-detail-panel\" aria-live=\"polite\"><p>Select a row to view its full evidence.</p></section>');
      $(table.table().container()).after(panel);
      const metaOptions = table.settings()[0].oInit.scidrMeta;
      if (metaOptions.evidence && Object.keys(metaOptions.evidence).length) {
        addFilters(metaOptions, [
          {name: 'ExtremeValue', title: 'Extreme value', blank: 'No'},
          {name: 'ReviewCategory', title: 'Anomaly type'},
          {name: 'ReviewTier', title: 'Severity tier'},
          {name: 'Methods', title: 'Method',
            members: function(x) {return String(x == null ? '' : x).split(', ');},
            label: function(v) {return metaOptions.methodNames[v] || v;}},
          {name: 'CodebookCategoryKeys', title: 'Codebook category',
            members: function(x) {try {return JSON.parse(String(x || '[]'));} catch(e) {return [];}}},
          {name: 'Analysis', title: 'Analysis'}]);
      }
      // Selecting a variable types its name into the measure search box instead of
      // applying a hidden filter, so the reader can see why rows disappeared and
      // clear it in place.
      if(metaOptions.columns.indexOf('IssueID') >= 0 && metaOptions.evidence) {
        const measureColumn = metaOptions.columns.indexOf('MeasureValue');
        root.on('scidrVariable', function(event, text) {
          if(measureColumn < 0) return;
          const value = text == null ? '' : String(text);
          const cell = $(table.table().header()).find('tr').last().find('td').eq(measureColumn);
          const input = cell.find('input');
          if(input.length) input.val(value);
          table.column(measureColumn).search(value).draw();
        });
      }
      let selectedIssue = null, selectionStarted = false;
      const select = function(tr) {
        const values = table.row(tr).data(); if (!values) return;
        $(table.rows().nodes()).removeClass('selected'); $(tr).addClass('selected');
        if(typeof Plotly !== 'undefined') panel.find('.js-plotly-plot').each(function() {Plotly.purge(this);});
        const meta = table.settings()[0].oInit.scidrMeta;
        selectedIssue = values[meta.columns.indexOf('IssueID')]; selectionStarted = true;
        panel.html(renderDetails(values, meta)); renderPlot(panel, values, meta); renderContributions(panel, values, meta);
      };
      table.on('draw.dt.scidrReview', function() {
        if(!selectionStarted || !metaOptions.evidence || !Object.keys(metaOptions.evidence).length) return;
        const issueColumn = metaOptions.columns.indexOf('IssueID');
        const filtered = table.rows({search:'applied'}).data().toArray();
        if(filtered.some(function(values) {return values[issueColumn] === selectedIssue;})) return;
        if(!filtered.length) {
          if(typeof Plotly !== 'undefined') panel.find('.js-plotly-plot').each(function() {Plotly.purge(this);});
          panel.html('<p>No findings match the current filters.</p>');selectedIssue=null;return;
        }
        const first=table.rows({search:'applied'}).indexes().toArray()[0];
        const node=table.row(first).node();if(node) select(node);
      });
      $(table.table().body()).on('click', 'tr', function() {select(this);});
      $(table.table().body()).on('keydown', 'tr', function(e) {if(e.key === 'Enter') {select(this);}});
      if (table.settings()[0].oInit.scidrMeta.evidence && Object.keys(table.settings()[0].oInit.scidrMeta.evidence).length && table.rows().count()) select(table.row(0).node());
      const resize = function() {table.columns.adjust(); if(table.settings()[0].oInit.scroller) table.scroller.measure();};
      $(table.table().container()).closest('.scidr-anomaly-explorer').find('details').on('toggle', resize);
    "))
    widget <- DT::datatable(df_Table, rownames = FALSE, colnames = headers, filter = "top", escape = FALSE,
      selection = "none", extensions = if (nrow(df_Table) > 50L) "Scroller" else character(), callback = callback,
      class = "stripe hover compact", options = list(
        deferRender = TRUE, scrollY = paste0(if (is.null(TableHeight)) 500 else TableHeight, "px"),
        scrollX = TRUE, scrollCollapse = TRUE, paging = nrow(df_Table) > 50L,
        scroller = if (nrow(df_Table) > 50L) list(displayBuffer = 5) else FALSE,
        pageLength = 50, dom = "fti", order = list(), autoWidth = FALSE,
        columnDefs = definitions, scidrMeta = meta,
        createdRow = DT::JS("function(row) {row.tabIndex = 0;}")))
    if (kind == "summary") widget$x$callback <- NULL
    widget
  }
  css <- htmltools::tags$style(htmltools::HTML("
    .scidr-anomaly-explorer { font-family:Arial,sans-serif; color:#263746; max-width:1400px; margin:24px auto; padding:0 20px; line-height:1.5; }
    .scidr-anomaly-explorer h2 { font-size:26px; margin:0 0 8px; }
    .scidr-anomaly-explorer h3 { font-size:18px; margin:28px 0 12px; }
    .scidr-anomaly-section { border:1px solid #d9e1e8; border-radius:6px; margin:14px 0; padding:12px; }
    .scidr-anomaly-section > summary { font-weight:600; cursor:pointer; padding:2px 0; }
    .scidr-anomaly-section[open] > summary { margin-bottom:12px; }
    .scidr-anomaly-explorer dt { font-weight:600; margin-top:10px; }
    .scidr-anomaly-explorer dd { margin:2px 0 10px; overflow-wrap:anywhere; }
    .scidr-anomaly-explorer .scidr-cell { height:76px; overflow:hidden; line-height:18px; min-width:100px; max-width:340px; white-space:normal; overflow-wrap:anywhere; }
    .scidr-measure { min-width:260px !important; }
    .scidr-field-label { display:block; max-height:36px; overflow:hidden; }
    .scidr-field-name { display:block; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
    .scidr-anomaly-explorer .scidr-reason { min-width:280px; }
    .scidr-anomaly-explorer .dataTables_wrapper { width:100%; }
    .scidr-anomaly-explorer table.dataTable thead { background:#edf2f6; }
    .scidr-anomaly-explorer table.dataTable tbody tr { cursor:pointer; }
    .scidr-detail-panel { border-left:3px solid #9aafbd; padding:10px 16px; margin:12px 0; background:#f7f9fb; }
    .scidr-anomaly-explorer code { font-size:11px; }
    .scidr-queue-filters { display:flex; gap:20px; margin:12px 0; }
    .scidr-overview { overflow-x:auto; }
    .scidr-raster { width:100%; }
    .scidr-anomaly-warning { border-left:4px solid #b2182b; background:#fdf2f2; padding:8px 12px; }
    .scidr-raster-note { font-size:13px; color:#5a6b7a; }
    .scidr-static-table { width:100%; border-collapse:collapse; font-size:13px; }
    .scidr-static-table th, .scidr-static-table td { padding:8px; border-bottom:1px solid #d9e1e8; text-align:left; }
    .scidr-anomaly-print-note { display:none; }
    .scidr-print-review { display:none; }
    @media print {
      .scidr-anomaly-explorer .datatables, .scidr-anomaly-explorer .scidr-detail-panel,
      .scidr-anomaly-explorer .scidr-anomaly-section, .scidr-queue-filters, .scidr-overview,
      .scidr-anomaly-explorer .scidr-raster { display:none !important; }
      .scidr-anomaly-explorer .scidr-print-review { display:block !important; }
      .scidr-overview { overflow-x:auto; }
    .scidr-static-table { width:100%; border-collapse:collapse; table-layout:fixed; font-size:8pt; }
      .scidr-static-table th, .scidr-static-table td { border:1px solid #bbb; padding:4px; overflow-wrap:anywhere; vertical-align:top; }
      .scidr-print-issues th:nth-child(5) { width:40%; }
      .scidr-static-table thead { display:table-header-group; }
      .scidr-static-table tr { break-inside:avoid; }

      .scidr-anomaly-explorer { margin:0; padding:0; font-size:10pt; }
      .scidr-anomaly-print-note { display:block; border:1px solid #999; padding:8px; }
      .scidr-anomaly-table .rt-search, .scidr-anomaly-table .rt-thead.-filters { display:none; }
      .scidr-anomaly-explorer .rt-table, .scidr-anomaly-explorer .rt-tbody { max-height:none !important; height:auto !important; overflow:visible !important; }
      .scidr-anomaly-explorer dl { max-height:none !important; overflow:visible !important; }
      .scidr-anomaly-section { break-inside:auto; }
    }
  "))
  # A method that stopped part-way reports fewer findings than it detected. The
  # page carries no technical diagnostics, so this would otherwise be invisible.
  df_Failures <- do.call(rbind, lapply(names(sources), function(analysis) {
    df_Diagnostics <- sources[[analysis]]$Diagnostics
    if (is.null(df_Diagnostics) || !nrow(df_Diagnostics)) return(NULL)
    df_Failed <- df_Diagnostics[df_Diagnostics$Issue == "backend_failure", , drop = FALSE]
    if (!nrow(df_Failed)) return(NULL)
    data.frame(Analysis = analysis,
      Method = unname(ifelse(is.na(df_Failed$Method), "unknown", method_names[df_Failed$Method])),
      Detail = df_Failed$Detail, stringsAsFactors = FALSE)
  }))
  neutral <- list(Metadata = list(Labels = character(), IDVariable = character()))
  Named <- function(source, variable) {
    label <- source$Metadata$Labels[variable]
    if (!length(label) || is.na(label) || !nzchar(label)) variable else unname(label)
  }
  df_Blocks <- do.call(rbind, lapply(names(sources), function(analysis) {
    df_Found <- ReviewObj$Sources[[analysis]]$Metadata$Settings$missingness$Blocks
    if (is.null(df_Found) || !nrow(df_Found)) return(NULL)
    data.frame(Analysis = analysis, Group = paste0("Group ", df_Found$Block), Variables = df_Found$Labels,
      Size = df_Found$Size, Started = df_Found$Started, Completed = df_Found$Completed,
      Completion = paste0(round(100 * df_Found$CompletionRate), "%"), stringsAsFactors = FALSE)
  }))
  df_Completeness <- do.call(rbind, lapply(names(sources), function(analysis) {
    source <- sources[[analysis]]
    variables <- intersect(source$Metadata$Variables, names(source$ReviewData))
    blanks <- vapply(variables, function(variable) sum(is.na(source$ReviewData[[variable]])), integer(1))
    if (!any(blanks > 0)) return(NULL)
    variables <- variables[blanks > 0]
    blanks <- blanks[blanks > 0]
    records <- nrow(source$ReviewData)
    grouped <- ReviewObj$Sources[[analysis]]$Metadata$Settings$missingness$Blocks
    members <- if (is.null(grouped) || !nrow(grouped)) character() else
      unlist(strsplit(grouped$Members, ", ", fixed = TRUE))
    reported <- unique(stats::na.omit(ReviewObj$Evidence$Variable[ReviewObj$Evidence$Analysis == analysis &
      ReviewObj$Evidence$Method == "missingness"]))
    df_Rows <- data.frame(Analysis = analysis, Measure = vapply(variables, function(v) Named(source, v), character(1)),
      Variable = variables, Blank = unname(blanks),
      Complete = paste0(round(100 * (records - blanks) / records), "%"),
      Group = ifelse(variables %in% members, "In a co-missing group", "Not grouped"),
      Reported = ifelse(variables %in% reported, "Yes", "No"), stringsAsFactors = FALSE)
    df_Rows[order(-df_Rows$Blank, df_Rows$Variable), , drop = FALSE]
  }))
  htmltools::tagList(plot_dependencies, css, htmltools::tags$div(`data-html-table-processing` = "none", class = "scidr-anomaly-explorer", htmltools::tags$h2(Title),
    htmltools::tags$details(class = "scidr-anomaly-section", open = if (Detail == "Full") NA else NULL,
      htmltools::tags$summary("How to read this page"),
      htmltools::tags$ol(
        htmltools::tags$li("Each row in \"Potential data issues\" is one thing worth a human look: a record, the field involved, its recorded value, its anomaly type, and why it was picked out."),
        htmltools::tags$li("Unexpected missingness means a field is blank where comparable records carry a value. Comparable is learned from the other columns, so a field that is deliberately skipped for part of the sample is not reported for those records."),
        htmltools::tags$li("A score says how unusual a record is relative to the other records in this dataset, on each method's own scale. Higher is more unusual. Scores from different methods are not comparable and are not probabilities."),
        htmltools::tags$li("A flag means the score passed a cutoff. Cutoffs are either learned from the data, set by the analyst, or a default top-percentile review budget - never a test of whether the value is wrong."),
        htmltools::tags$li("Nothing here is a confirmed error. Plenty of unusual values are correct, and real errors can look ordinary. Select a row to see the comparison group and supporting evidence before judging it."),
        htmltools::tags$li("To record decisions, work through the accompanying Excel workbook rather than this page; this page is read-only."))),
    htmltools::tags$p(ReviewObj$Metadata$Notice),
    if (!is.null(variables)) htmltools::tags$p(paste0(
      "This page lists findings for ", length(variables), " selected variable(s) only, plus findings that concern a whole record. ",
      "Detection ran on the full selection, so each finding shown is still conditioned on every screened column, and the coverage summary below describes the whole screen, not this subset.")),
    if (!is.null(df_Failures)) htmltools::tagList(
      htmltools::tags$p(class = "scidr-anomaly-warning",
        paste0("Incomplete screening: ", nrow(df_Failures),
          " method run(s) stopped before finishing, so findings from them may be missing. Rerun detection, and report the message below if it persists.")),
      PrintTable(data.frame(Analysis = df_Failures$Analysis, Method = df_Failures$Method,
        "What stopped it" = df_Failures$Detail, check.names = FALSE))),
    htmltools::tags$p(class = "scidr-anomaly-print-note", "Print view contains the complete issue queue and coverage summary. Technical scrolling tables are omitted; use Excel for full evidence and review decisions."),
    htmltools::tags$details(class = "scidr-anomaly-section", htmltools::tags$summary("What the methods mean"),
      htmltools::tags$p("Requested methods appear first. Methods generate findings; variable-level evidence supports them, and isolation column screening suggests inspection leads. None estimates the probability of an entry error. Cross-sectional defaults are Conditional pattern, Isolation forest, Date screening, and Unexpected missingness."),
      PrintTable(data.frame(Method = df_Methods$Method, Meaning = df_Methods$Meaning, "Applies to" = df_Methods$AppliesTo, check.names = FALSE)),
      htmltools::tags$p("Every finding also carries one anomaly type, in the \"Anomaly type\" column. Sort or filter either table by it to work through one kind of problem at a time."),
      PrintTable(data.frame("Anomaly type" = ScidrAnomalyTypes()$Type, Meaning = ScidrAnomalyTypes()$Definition, check.names = FALSE)),
      htmltools::tags$p("A value can be so far outside its column that no judgement is needed: the \"Extreme value\" column marks a record sitting at least ten times further from the column median than any other record does. Those are listed first, ahead of every tier, because the kind of finding is the same whether a value is slightly or absurdly out of range. The test needs at least ten other numeric values in the column and is not applied to dates, clock times or blanks."),
      htmltools::tags$p("Severity tier is the kind of finding, ranked. It drives the suggested review order and is a sortable column on both tables: tier 1 first, tier 8 last. Scores from different methods are not comparable, so the tier orders findings across methods and the within-method percentile only separates findings of the same tier. The ranking is an unvalidated review aid, not error probability or clinical urgency."),
      PrintTable(data.frame(Tier = ScidrAnomalySeverityTable()$Tier, "Anomaly types" = ScidrAnomalySeverityTable()$Types,
        "Why it ranks here" = ScidrAnomalySeverityTable()$Meaning, check.names = FALSE))),
    htmltools::tags$details(class = "scidr-anomaly-section",
      htmltools::tags$summary("Groups of variables that go missing together"),
      htmltools::tags$p("Unexpected missingness compares a blank against records that are comparable for a reason the data supplies, never against the column as a whole. There are two such reasons. A rule learned from another column's values (\"every record with Sex = Female has a value\") is a skip rule the record appears to break. A group of variables that go blank together is a form: the groups below were found by clustering variables on how similarly they go missing, with no knowledge of what they measure, and a blank is reported when a record carries most of a group but not all of it, while nearly every record carrying any of that group carries all of it."),
      htmltools::tags$p("Check the membership before trusting what the findings say. A group that mixes variables from unrelated forms is a sign the clustering has over-merged; narrow it with method_args = list(missingness = list(block_similarity = 0.9))."),
      if (is.null(df_Blocks)) htmltools::tags$p("No groups were found: no set of variables goes blank together often enough in this data. Missingness findings here, if any, rest on a learned value rule instead.")
        else PrintTable(data.frame("Analysis" = df_Blocks$Analysis, "Group" = df_Blocks$Group,
          "Variables" = df_Blocks$Variables, "Records carrying any" = df_Blocks$Started,
          "Carrying all" = df_Blocks$Completed, "Completion" = df_Blocks$Completion, check.names = FALSE)),
      htmltools::tags$p("A record missing most of a group is read as a form that was not done and is not reported. Blanks belonging to no group and explained by no value rule are not reported either; they are listed under \"Completeness\" below as a fact about the column, which is all that can honestly be said about them.")),
    if (!is.null(df_Completeness)) htmltools::tags$details(class = "scidr-anomaly-section",
      htmltools::tags$summary(paste0("Completeness (", nrow(df_Completeness), " columns with blanks)")),
      htmltools::tags$p("Every column carrying at least one blank, whether or not any blank was reported as a finding. \"Reported\" says whether this column produced a missingness finding; a No with many blanks is not a clean column, only one whose blanks nothing in the data explains."),
      PrintTable(data.frame("Analysis" = df_Completeness$Analysis, "Measure" = df_Completeness$Measure,
        "Variable" = df_Completeness$Variable, "Blank" = df_Completeness$Blank,
        "Complete" = df_Completeness$Complete, "Grouping" = df_Completeness$Group,
        "Reported" = df_Completeness$Reported, check.names = FALSE))),
    htmltools::tags$h3("Anomaly map"),
    htmltools::tags$p("One cell per modelled variable and record, in saved source-row order. Grey means no anomaly was reported for that value; a colour marks the kind of anomaly carried by the cell. Hover a cell for the record, the recorded value and the reason. A cell carrying several findings takes the colour of its most serious one. This is a map of what was reported, not of what is wrong."),
    if (IncludePlots) lapply(names(sources), function(analysis) ScidrAnomalyRasterTag(ReviewObj, display, analysis, length(sources) > 1L))
      else htmltools::tags$p("The map needs interactive plots. Use IncludePlots = TRUE and install plotly."),
    htmltools::tags$h3("Variables to review"),
    htmltools::tags$p("Ordered by the strongest supported finding, not error probability. Counts identify distinct affected source records, not repeated evidence. Red marks review values in within-stratum column distributions; these previews do not reproduce conditional subgroups. Select a variable to filter the findings below."),
    htmltools::tags$button(type = "button", class = "scidr-variable-reset", "Show all findings"),
    Table(variable_review$Table, kind = "variables", source = neutral),
    htmltools::tags$h3("Potential data issues"),
    if (!nrow(ReviewObj$Issues)) htmltools::tags$p("No flagged issues were added. Score-only methods may still have ranked candidates; use CreateAnomalyReviewQueue(include_ranked = TRUE) to add an explicitly exploratory queue."),
    if (!isTRUE(ReviewObj$Metadata$IncludeRanked) && any(df_Coverage$`Scored only` > 0)) htmltools::tags$p("This saved queue omits some score-only candidates. Its membership and decisions are unchanged. To create a new queue, use CreateAnomalyReviewQueue(AnomalyObj, include_ranked = TRUE); retain this queue for its existing decisions."),
    Table(display$Issues, kind = "queue", source = neutral, issue_evidence = display$Evidence),
    htmltools::tags$div(class = "scidr-print-review", htmltools::tags$h3("Coverage summary"), PrintTable(df_Coverage),
      htmltools::tags$h3("Variables to review"), PrintTable(PrintVariableReview(variable_review$Table)),
      htmltools::tags$h3("Complete review queue"), PrintTable(df_Print, queue = TRUE)),
    htmltools::tags$p(if (IncludePlots) "Interactive plots are enabled. The first issue is selected automatically; select another row for its plot. Exploratory records offer selectable comparison variables; plots do not attribute responsibility to those fields." else "Plots are off. Use IncludePlots = TRUE and install plotly to view interactive figures."),
    if (any(vapply(sources, function(object) inherits(object, "LongitudinalAnomalyResult"), logical(1))))
      htmltools::tags$p("Longitudinal findings carry a retained history; select one to view it. PlotDataAnomalyTrajectories() also returns ordinary ggplots for selected histories."),
    htmltools::tags$details(class = "scidr-anomaly-section", htmltools::tags$summary("Coverage"),
      htmltools::tags$div(class = "scidr-overview", PrintTable(df_Coverage)))))
}

#' Save an anomaly explorer as a standalone HTML file
#'
#' Write the [ExploreDataAnomalies()] dashboard to a single HTML file that a
#' collaborator can open in any browser without R, Quarto, or a network
#' connection. Interactive dependencies are inlined when pandoc is available.
#' @param AnomalyObj A detection result or an `AnomalyReviewQueue`, as in
#'   [ExploreDataAnomalies()].
#' @param file Destination `.html` path in an existing directory.
#' @param overwrite Allow replacement of an existing file? Default FALSE.
#' @param self_contained Inline scripts, styles and data into one file using
#'   pandoc? Default TRUE. Without pandoc, a sidecar `_files` directory of
#'   dependencies is kept beside the HTML and a message says so; both must then
#'   be sent together.
#' @param ... Passed to [ExploreDataAnomalies()], for example `Title`,
#'   `IncludePlots`, `Detail` or `codebook`.
#' @return The destination path, invisibly. No input data is modified.
#' @family data anomalies
#' @seealso [ExploreDataAnomalies()], [CreateDataAnomalyReport()],
#'   [ExportAnomalyReviewQueue()]
#' @examples
#' set.seed(24)
#' df_Example <- data.frame(ID = 1:60, Visit = as.Date("2024-01-01") + c(0:58, 900))
#' result <- DetectDataAnomalies(df_Example, "Visit", id_var = "ID", methods = "dates")
#' if (requireNamespace("DT", quietly = TRUE)) {
#'   SaveDataAnomalyExplorer(result, file.path(tempdir(), "AnomalyExplorer.html"))
#' }
#' @export
SaveDataAnomalyExplorer <- function(AnomalyObj, file, overwrite = FALSE,
    self_contained = TRUE, ...) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) stop("file must be one path.")
  if (!dir.exists(dirname(file))) stop("Directory does not exist: ", dirname(file))
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) stop("overwrite must be TRUE or FALSE.")
  if (!is.logical(self_contained) || length(self_contained) != 1L || is.na(self_contained)) stop("self_contained must be TRUE or FALSE.")
  if (file.exists(file) && !overwrite) stop("File exists; use overwrite = TRUE: ", file)
  ScidrAnomalyRequireDependencies(ScidrAnomalyDependencies(explorer = TRUE))
  page <- ExploreDataAnomalies(AnomalyObj, ...)
  libdir <- paste0(tools::file_path_sans_ext(basename(file)), "_files")
  htmltools::save_html(page, file = file, libdir = libdir)
  if (self_contained) {
    if (requireNamespace("rmarkdown", quietly = TRUE)) ScidrAnomalyFindPandoc()
    if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available("1.17")) {
      rmarkdown::pandoc_self_contained_html(file, file)
      unlink(file.path(dirname(file), libdir), recursive = TRUE)
    } else {
      message("pandoc is unavailable, so dependencies stay in '", libdir,
        "' beside the HTML file. Send the file and that directory together, or zip the pair.")
    }
  }
  message("Wrote explorer: ", file)
  invisible(file)
}
