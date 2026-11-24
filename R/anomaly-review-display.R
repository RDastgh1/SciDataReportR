# Recorded values are shown at full stored precision otherwise, so a reading like
# 4980.49871276928 crowds out the sentence around it. Round for display only;
# the stored evidence keeps every digit. Text, dates and categories pass through,
# and whole numbers are never reshaped into something that looks rounded.
ScidrReviewRoundedValue <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(as.character(x))
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value) || !is.finite(value)) return(as.character(x))
  # Outside this range plain notation is less readable than scientific, so leave
  # the default formatting to decide.
  plain <- abs(value) >= 1e-4 && abs(value) < 1e15
  if (value == round(value) && plain) return(format(value, trim = TRUE, scientific = FALSE))
  format(if (abs(value) < 1) signif(value, 3) else round(value, 2), trim = TRUE, scientific = !plain)
}

# Display projections are shared by HTML and Excel; immutable queue evidence stays intact.
ScidrReviewNarratives <- function(object, concise = TRUE) {
  display_sources <- stats::setNames(lapply(names(object$Sources), function(analysis) ScidrReviewDisplaySource(object, analysis)), names(object$Sources))
  df_Evidence <- ScidrReviewCategoryFields(object, object$Evidence)
  df_Evidence$WhyFlagged <- character(nrow(df_Evidence))
  df_Evidence$ComparisonDetails <- character(nrow(df_Evidence))
  Text <- function(x) paste(ScidrReviewFormat(x), collapse = ", ")
  Present <- function(x) length(x) && !all(is.na(x)) && any(nzchar(as.character(x)))
  Number <- function(x) if (!Present(x)) "unavailable" else format(signif(x, 5), trim = TRUE, scientific = FALSE)
  for (i in seq_len(nrow(df_Evidence))) {
    e <- lapply(df_Evidence, function(x) if (is.list(x)) x[[i]] else x[i])
    source <- display_sources[[e$Analysis]]
    Label <- function(variable) {
      label <- source$Metadata$Labels[variable]
      if (!Present(label)) variable else paste0(unname(label), " (", variable, ")")
    }
    statistics <- e$GroupStatistics
    native <- e$NativeEvidence
    conditions <- e$Conditions
    condition_text <- vapply(conditions, function(condition) {
      variable <- condition$column
      Format <- function(value) {
        conversion <- source$Metadata$DateConversions[[variable]]
        declaration <- source$Metadata$TimeConversions
        kind <- declaration$TimeType[match(variable, declaration$Variable)]
        if (is.numeric(value) && length(kind) && !is.na(kind) && kind == "duration") return(ScidrTimeFormat(value, "duration"))
        if (!is.null(conversion) && is.numeric(value)) ScidrAnomalyFormatDate(value, conversion) else if (is.numeric(value)) Number(value) else Text(value)
      }
      paste0(Label(variable), " ", Text(condition$comparison), " ", Format(condition$value_comp),
        if (Present(condition$value_this)) paste0("; recorded: ", Format(condition$value_this)) else "")
    }, character(1))
    group <- if (length(condition_text)) paste0("the learned subgroup defined by ", paste(condition_text, collapse = "; ")) else "the observed column within this analysis group"
    details <- character()
    reason <- e$Explanation
    if (e$Method == "conditional") {
      if (Present(e$ReferenceUpper) || Present(e$ReferenceLower)) {
        value <- suppressWarnings(as.numeric(e$Observed))
        upper <- suppressWarnings(as.numeric(e$ReferenceUpper))
        lower <- suppressWarnings(as.numeric(e$ReferenceLower))
        # Date bounds are ISO strings and can be compared in their original class.
        if (Present(e$Variable) && ScidrAnomalyIsDate(source$ReviewData[[e$Variable]])) {
          value <- as.numeric(as.POSIXct(e$Observed, tz = "UTC"))
          upper <- if (Present(e$ReferenceUpper)) as.numeric(as.POSIXct(e$ReferenceUpper, tz = "UTC")) else NA_real_
          lower <- if (Present(e$ReferenceLower)) as.numeric(as.POSIXct(e$ReferenceLower, tz = "UTC")) else NA_real_
        }
        if (is.finite(value) && is.finite(upper) && value > upper) reason <- paste0("Recorded value ", ScidrReviewRoundedValue(e$Observed), " is unusually high: above the learned comparison boundary of ", ScidrReviewRoundedValue(e$ReferenceUpper), ".")
        else if (is.finite(value) && is.finite(lower) && value < lower) reason <- paste0("Recorded value ", ScidrReviewRoundedValue(e$Observed), " is unusually low: below the learned comparison boundary of ", ScidrReviewRoundedValue(e$ReferenceLower), ".")
      } else if (Present(statistics$categs_common)) {
        reason <- paste0("Recorded response '", e$Observed, "' is uncommon in ", group, ". Common response(s): ", Text(statistics$categs_common),
          if (Present(statistics$pct_common)) paste0(" (", Number(100 * statistics$pct_common), "% of the reference group)") else "", ".")
      }
      details <- c(details, paste0("Comparison: ", group, "."))
      if (Present(e$ReferenceN)) details <- c(details, paste0("Reference observations: ", e$ReferenceN, "."))
      if (Present(statistics$mean)) {
        conversion <- source$Metadata$DateConversions[[e$Variable]]
        declaration <- source$Metadata$TimeConversions
        kind <- declaration$TimeType[match(e$Variable, declaration$Variable)]
        mean <- if (length(kind) && !is.na(kind) && kind == "duration") ScidrTimeFormat(statistics$mean, "duration") else if (is.null(conversion)) Number(statistics$mean) else ScidrAnomalyFormatDate(statistics$mean, conversion)
        details <- c(details, paste0("Reference group mean: ", mean, "."))
      }
      if (!length(statistics)) details <- c(details, "Structured comparison statistics were not retained in this older result; the original explanation is shown.")
    } else if (e$Method == "missingness") {
      reference <- if (length(condition_text)) group else "all eligible reference records in this stratum"
      details <- c(details, paste0("Comparison: ", reference, "."),
        paste0("Reference observations: ", Number(e$ReferenceN), "; observed: ", Number(statistics$n_observed),
          "; missing: ", Number(statistics$n_missing), "."),
        if (length(source$Metadata$IDVariable)) "Reference counts exclude this participant's records; missing IDs use row exclusion." else "Reference counts exclude the evaluated source row.",
        "The observed fraction describes a selected comparison group, not an error probability. Skip logic, study design and selective testing may explain a blank.")
    } else if (e$Method == "dates") {
      if (e$Kind == "date_interval") {
        reason <- paste0("The signed interval ", Label(e$OtherVariable), " minus ", Label(e$Variable), " is ", Number(e$ActualIntervalDays),
          " days, outside the learned screening boundaries [", e$ReferenceLower, ", ", e$ReferenceUpper, "] days.")
        if (identical(e$ReferenceLower, "0") && identical(e$ReferenceUpper, "0")) reason <- paste0(
          "These dates differ by ", Number(abs(e$ActualIntervalDays)), " days. The learned screening boundaries indicate matching dates. Verify whether these fields represent the same event and check both source dates.")
        details <- c(details, paste0("Direction: ", Label(e$OtherVariable), " minus ", Label(e$Variable), "; signed interval: ", Number(e$ActualIntervalDays), " days."))
        details <- c(details, paste0("Dates: ", e$Observed, " to ", e$OtherObserved, ". A negative interval is not automatically invalid."))
      } else reason <- paste0("Recorded date ", e$Observed, " is outside the learned calendar screening boundaries [", e$ReferenceLower, ", ", e$ReferenceUpper, "].")
      details <- c(details, paste0("Reference observations: ", e$ReferenceN, ". Boundaries describe the observed distribution, not clinical normality."))
    } else if (e$Method == "change") {
      reason <- paste0("Value changed from ", Number(as.numeric(e$OtherObserved)), " to ", Number(as.numeric(e$Observed)), " over ", Text(e$Elapsed), " ", e$ElapsedUnit,
        "; the change rate is ", Number(native$Rate), " per time unit. Its absolute deviation from the reference median is ", Number(e$Score), " MAD units.")
      details <- c(details, "This identifies an unusual pair of observations; either endpoint, both, or neither may need correction.")
    } else if (e$Method == "transition") {
      reason <- paste0("Transition '", e$OtherObserved, "' to '", e$Observed, "' occurred ", Text(native$TransitionCount), " times in ", Text(native$ConditioningCount),
        " reference transitions starting at '", e$OtherObserved, "', excluding this participant.")
      details <- c(details, paste0("Reference participants: ", Text(native$ReferenceParticipants), ". Score = 1 minus the observed conditional transition frequency; this is not an error probability."))
    } else if (e$Method == "trajectory") {
      reason <- paste0("Recorded value ", e$Observed, " differs from its fitted trajectory value ", e$Expected, "; absolute Pearson residual = ", Number(e$Score), ".")
      details <- c(details, "This is an in-sample random-intercept model diagnostic, not a prediction interval or error probability.")
    } else if (e$Method == "cellwise") {
      reason <- paste0("Recorded value ", e$Observed, " was flagged by the cellwise detector; absolute standardized residual = ", Number(e$Score),
        if (Present(e$Expected)) paste0("; expected value = ", e$Expected) else "", ".")
    } else if (e$Method %in% c("isolation", "local", "distance")) {
      reason <- paste0("The overall combination of selected values is unusual under the ", e$Method, " method (score ", Number(e$Score),
        "). This score does not identify a single responsible field.")
    }
    declaration <- source$Metadata$TimeConversions
    kind <- declaration$TimeType[match(e$Variable, declaration$Variable)]
    if (length(kind) && !is.na(kind) && kind == "clock" && e$Method == "conditional") {
      reason <- paste0("Recorded clock time ", e$Observed, " falls in an unusual clock region",
        if (Present(statistics$ClockRegion)) paste0(": ", statistics$ClockRegion) else " on the circular model representation", ".")
      details <- c(details, "Clock regions are learned comparison boundaries, not clinical normal ranges. Sine/cosine component statistics remain under original backend evidence.")
    }
    if (length(kind) && !is.na(kind) && kind == "duration" && e$Method == "conditional") {
      value <- ScidrTimeSeconds(source$ReviewData[[e$Variable]], declaration[match(e$Variable, declaration$Variable), ], strict = FALSE)[match(e$SourceRow, source$ReviewData$.Row)]
      if (Present(statistics$upper_thr) && is.finite(value) && value > statistics$upper_thr)
        reason <- paste0("Recorded duration ", e$Observed, " is unusually high: above the learned comparison boundary of ", ScidrTimeFormat(statistics$upper_thr, "duration"), ".")
      else if (Present(statistics$lower_thr) && is.finite(value) && value < statistics$lower_thr)
        reason <- paste0("Recorded duration ", e$Observed, " is unusually low: below the learned comparison boundary of ", ScidrTimeFormat(statistics$lower_thr, "duration"), ".")
    }
    if (length(kind) && !is.na(kind) && kind == "clock" && e$Method == "cellwise") {
      reason <- paste0("Recorded clock time ", e$Observed, " was flagged on its circular representation; strongest absolute standardized component residual = ", Number(e$Score), ". No replacement clock time is inferred.")
    }
    if (!Present(reason)) reason <- "Selected for review; a more detailed explanation was not retained."
    if (e$Method == "conditional" && !Present(statistics$categs_common) && length(statistics)) {
      reason <- paste(reason, paste0("Compared with ", group, if (Present(e$ReferenceN)) paste0("; reference n = ", e$ReferenceN) else "", "."))
    }
    if (e$Method == "change") {
      reference <- source$IntervalResults
      if (all(c("Variable", "Stratum", "Rate") %in% names(reference)) && Present(native$Stratum)) {
        rates <- reference$Rate[reference$Variable == e$Variable & reference$Stratum == native$Stratum & reference$Method == "change"]
        rates <- rates[is.finite(rates)]
        if (length(rates)) details <- c(details, paste0("Reference median rate: ", Number(stats::median(rates)),
          "; MAD: ", Number(stats::mad(rates)), "; intervals: ", length(rates), "."))
      }
    }
    if (e$Method == "conditional" && Present(e$Variable) && !e$Variable %in% source$Metadata$TimeConversions$Variable && is.numeric(source$ReviewData[[e$Variable]])) {
      distinct <- unique(stats::na.omit(source$ReviewData[[e$Variable]]))
      if (length(distinct) <= 3L && all(distinct %in% c(-1, 0, 1)) &&
          (Present(e$ReferenceLower) || Present(e$ReferenceUpper))) {
        boundary <- if (Present(e$ReferenceUpper)) e$ReferenceUpper else e$ReferenceLower
        reason <- paste0("Recorded value ", e$Observed, " differs from the learned common-value boundary of ", boundary,
          ". This column contains only ", paste(sort(distinct), collapse = ", "),
          "; confirm its coding and recorded value. A nonzero value or change can be valid.")
        details <- c(details, "Low-cardinality numeric data were screened numerically; no categorical recoding or clinical rule was inferred.")
      }
    }
    cutoff <- source$Metadata$Thresholds[e$Method]
    if (length(cutoff) && is.finite(cutoff)) details <- c(details, paste0("Explicit score cutoff: ", Number(cutoff), "; larger scores are more anomalous."))
    if (concise) details <- c(reason, details)
    df_Evidence$WhyFlagged[i] <- if (concise) ScidrReviewShortReason(e, source, reason) else reason
    df_Evidence$ComparisonDetails[i] <- paste(details, collapse = "\n")
  }
  df_Issues <- ScidrReviewCategoryFields(object, object$Issues)
  df_Issues$WhyFlagged <- vapply(df_Issues$IssueID, function(id) if (concise) df_Evidence$WhyFlagged[which(df_Evidence$IssueID == id)[1L]] else paste(unique(df_Evidence$WhyFlagged[df_Evidence$IssueID == id]), collapse = "\n"), character(1))
  df_Issues$ComparisonDetails <- vapply(df_Issues$IssueID, function(id) paste(unique(df_Evidence$ComparisonDetails[df_Evidence$IssueID == id]), collapse = "\n\n"), character(1))
  ScidrReviewOrder(object, list(Issues = df_Issues, Evidence = df_Evidence))
}

ScidrAnomalyDetailJS <- function() 'function(values, meta) {
  const esc = function(x) {return String(x == null ? "" : x).replace(/[&<>"\x27]/g, function(c) {
    return {"&":"&amp;", "<":"&lt;", ">":"&gt;", "\\\"":"&quot;", "\x27":"&#39;"}[c];
  });};
  const row = {}; meta.columns.forEach(function(name, i) {row[name] = values[i];});
  const present = function(x) {return x != null && String(x) !== "" && String(x) !== "NA" && String(x) !== "[]" && String(x) !== "{}";};
  const field = function(name, value, labels) {
    if (!present(value)) return "";
    return "<dt>" + esc((labels || meta.labels)[name] || name) + "</dt><dd style=\\"white-space:pre-wrap\\">" + esc(value) + "</dd>";
  };
  // Backend evidence is stored as JSON. Dumping the raw string on the page asks
  // the reader to parse it; the same object as named rows does not.
  const nativeNames = {".Row":"Source row", "CheckType":"Check", "ValueDays":"Interval (days)", "N":"Reference observations",
    "MinimumDays":"Smallest observed interval (days)", "MaximumDays":"Largest observed interval (days)",
    "LowerFenceDays":"Lower boundary (days)", "UpperFenceDays":"Upper boundary (days)",
    "MinimumDate":"Earliest reference date", "MaximumDate":"Latest reference date",
    "LowerDate":"Lower boundary date", "UpperDate":"Upper boundary date", "Rule":"Boundary rule",
    "Stratum":"Analysis group", "Status":"Evaluation status", "Flag":"Passed the cutoff",
    "NativeEvidence":"Backend evidence", "Conditions":"Comparison group", "GroupStatistics":"Comparison group statistics",
    "column":"Column", "comparison":"Comparison", "value_comp":"Compared against", "value_this":"This record"};
  // Fields already stated above the fold are not repeated inside the backend object.
  const hiddenNames = ["Method", "Label", "OtherLabel", "Explanation", "Variable", "OtherVariable", "Observed", "OtherObserved"];
  const spaced = function(name) {
    if (nativeNames[name]) return nativeNames[name];
    const label = meta.labels[name];
    if (label && label !== name) return label;
    return String(name).replace(/([a-z0-9])([A-Z])/g, "$1 $2");
  };
  const structured = function(name, value) {
    if (!present(value)) return "";
    let parsed = value;
    if (typeof value === "string") {try {parsed = JSON.parse(value);} catch (e) {return field(name, value);}}
    if (parsed === null || typeof parsed !== "object") return field(name, value);
    const entries = Array.isArray(parsed) ? parsed.map(function(x, i) {return ["Item " + (i + 1), x];}) :
      Object.keys(parsed).map(function(key) {return [key, parsed[key]];});
    const kept = entries.filter(function(entry) {return present(entry[1]) && hiddenNames.indexOf(entry[0]) < 0;});
    if (!kept.length) return "";
    let inner = "";
    kept.forEach(function(entry) {
      inner += (entry[1] !== null && typeof entry[1] === "object") ? structured(entry[0], entry[1]) :
        "<dt>" + esc(spaced(entry[0])) + "</dt><dd>" + esc(typeof entry[1] === "boolean" ? (entry[1] ? "Yes" : "No") : entry[1]) + "</dd>";
    });
    return "<dt>" + esc(spaced(name)) + "</dt><dd><dl>" + inner + "</dl></dd>";
  };
  let html = "<h4>Selected finding</h4><dl>";
  const primary = ["StudyID", "SourceRow", ".Row", "OtherRow", "Analysis", "SuggestedReviewOrder", "ReviewTier", "OrderRationale", "FindingType", "ReviewCategory", "CodebookCategory", "CategoryBasis", "CategoryVariables", "ConditionCategories", "ReviewStatus", "Priority", "Label", "Variable", "Observed", "OtherLabel", "OtherVariable", "OtherObserved", "WhyFlagged", "ComparisonDetails", "SuggestedCheck", "Context", "ProposedCorrection", "SourceChecked", "Reviewer", "ReviewedOn", "Notes"];
  primary.forEach(function(name) {if (name !== "Analysis" || meta.multipleAnalyses !== false) html += field(name, row[name]);});
  const evidence = (meta.evidence || {})[row.IssueID];
  if (evidence) {
    html += "</dl><details><summary>Method scores and technical evidence</summary><h4>Supporting evidence</h4>";
    evidence.forEach(function(e) {
      html += "<dl>";
      ["Method", "ReferenceN", "ReferenceLower", "ReferenceUpper", "Expected", "Elapsed", "ElapsedUnit", "Score", "Rank"].forEach(function(name) {
        html += field(name, name === "Method" ? (meta.methodNames[e[name]] || e[name]) : e[name]);
      });
      html += "</dl><details><summary>Original backend evidence</summary><dl>";
      ["ConditionCategories", "Explanation"].forEach(function(name) {html += field(name, e[name]);});
      ["Conditions", "GroupStatistics", "NativeEvidence"].forEach(function(name) {html += structured(name, e[name]);});
      html += "</dl></details>";
    });
    html += "</details><dl>";
  }
  if (meta.reviewRows) {
    const index = meta.reviewRows.indexOf(Number(row[".Row"]));
    if (index >= 0) {
      html += "</dl><h4>Recorded source values</h4><dl>";
      Object.keys(meta.reviewValues).forEach(function(name) {
        const label = meta.labels[name];
        const labels = {}; labels[name] = label && label !== name ? label + " (" + name + ")" : name;
        html += field(name, meta.reviewValues[name][index], labels);
      });
    }
  }
  html += "</dl><details><summary>Tracking details</summary><dl>";
  Object.keys(row).filter(function(name) {return name !== "CodebookCategoryKeys" && primary.indexOf(name) < 0;}).forEach(function(name) {html += field(name, row[name]);});
  return html + "</dl></details>";
}'
