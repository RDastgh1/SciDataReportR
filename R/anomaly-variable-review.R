# Interval findings read as raw numbers unless the comparison is spelled out:
# a signed day count plus "bounds 0 to 0" hides that the two dates normally fall
# on the same day, which is the whole point of the finding.
ScidrReviewIntervalReason <- function(variable, other, days, lower, upper) {
  Days <- function(x) format(abs(as.numeric(x)), trim = TRUE)
  days <- suppressWarnings(as.numeric(days))
  lower <- suppressWarnings(as.numeric(lower))
  upper <- suppressWarnings(as.numeric(upper))
  gap <- if (!length(days) || is.na(days)) paste0(other, " and ", variable, " differ by an unknown interval") else
    if (days > 0) paste0(other, " falls ", Days(days), " days after ", variable) else
    if (days < 0) paste0(other, " falls ", Days(days), " days before ", variable) else
    paste0(other, " and ", variable, " fall on the same day")
  expected <- if (!length(lower) || !length(upper) || is.na(lower) || is.na(upper))
    "no comparison range was available" else
    if (isTRUE(all.equal(lower, upper))) {
      if (isTRUE(all.equal(lower, 0))) "in comparable records the two dates fall on the same day" else
      if (lower > 0) paste0("comparable records place it ", Days(lower), " days later") else
        paste0("comparable records place it ", Days(lower), " days earlier")
    } else paste0("comparable records fall between ", format(lower, trim = TRUE), " and ",
      format(upper, trim = TRUE), " days apart")
  paste0(gap, "; ", expected, ".")
}

# Concise display text is deliberately separate from immutable evidence.
ScidrReviewShortReason <- function(e, source, fallback) {
  Present <- function(x) length(x) == 1L && !is.na(x) && nzchar(as.character(x))
  Number <- function(x) if (Present(x)) ScidrReviewRoundedValue(x) else "unavailable"
  Value <- function(x) if (Present(x)) ScidrReviewRoundedValue(x) else as.character(x)
  if (e$Method == "missingness") {
    statistics <- e$GroupStatistics
    if (length(statistics$block_members)) return(ScidrAnomalyMissingnessBlockReason(
      statistics$observed_in_block, statistics$block_size, statistics$n_observed, e$ReferenceN,
      statistics$block_members, source$Metadata$Labels))
    return(ScidrAnomalyMissingnessReason(statistics$n_observed, e$ReferenceN, statistics$comparison_text))
  }
  if (e$Method %in% c("isolation", "local", "distance")) return(paste0("Uncommon combination of values across the record (", unname(ScidrAnomalyMethods()$Method[match(e$Method, ScidrAnomalyMethods()$Key)]), ")."))
  if (e$Method == "dates" && e$Kind == "date_interval")
    return(ScidrReviewIntervalReason(e$Variable, e$OtherVariable, e$ActualIntervalDays, e$ReferenceLower, e$ReferenceUpper))
  if (e$Method == "dates") return(paste0(Value(e$Observed), " is outside the learned date boundaries."))
  if (e$Method == "change") return(paste0(Value(e$OtherObserved), " \u2192 ", Value(e$Observed), " over ", Number(e$Elapsed), " ", e$ElapsedUnit, "; unusual change rate."))
  if (e$Method == "transition") return(paste0("Uncommon transition: ", Value(e$OtherObserved), " \u2192 ", Value(e$Observed), "."))
  if (e$Method == "trajectory") return(paste0(Value(e$Observed), " differs from its in-sample fitted trajectory."))
  if (e$Method == "cellwise") return(paste0(Value(e$Observed), " is inconsistent with related numeric fields."))
  if (e$Method == "conditional") {
    declaration <- source$Metadata$TimeConversions
    kind <- declaration$TimeType[match(e$Variable, declaration$Variable)]
    if (length(kind) && !is.na(kind) && kind == "clock") return(paste0(Value(e$Observed), " falls in an unusual clock region."))
    if (!length(e$GroupStatistics)) return(fallback)
    group <- if (length(e$Conditions)) "subgroup" else "column"
    value <- suppressWarnings(as.numeric(e$Observed))
    statistics <- e$GroupStatistics
    if (length(kind) && !is.na(kind) && kind == "duration") {
      return(paste0(Value(e$Observed), " is unusual for its ", group, " distribution."))
    }
    upper <- suppressWarnings(as.numeric(e$ReferenceUpper))
    lower <- suppressWarnings(as.numeric(e$ReferenceLower))
    if (length(value) && is.finite(value) && length(upper) && is.finite(upper) && value > upper)
      return(paste0(Value(e$Observed), " is above the ", group, " boundary of ", Number(e$ReferenceUpper), "."))
    if (length(value) && is.finite(value) && length(lower) && is.finite(lower) && value < lower)
      return(paste0(Value(e$Observed), " is below the ", group, " boundary of ", Number(e$ReferenceLower), "."))
    return(paste0(if (Present(e$Observed)) Value(e$Observed) else "This value", " is uncommon in its ", group, " distribution."))
  }
  "Selected for review; see the retained evidence."
}

# One strongest finding says little about a variable with many findings. Summarise
# what kinds of problem it has and how often, with a short example of each.
ScidrReviewCommonReasons <- function(categories, examples, limit = 3L) {
  counts <- sort(table(categories), decreasing = TRUE)
  parts <- vapply(names(counts)[seq_len(min(limit, length(counts)))], function(category) {
    example <- examples[categories == category]
    example <- example[!is.na(example) & nzchar(example)]
    label <- paste0(category, if (counts[[category]] > 1L) paste0(" (", counts[[category]], " records)") else "")
    if (!length(example)) label else paste0(label, ": ", example[1])
  }, character(1))
  remainder <- length(counts) - min(limit, length(counts))
  paste0(paste(parts, collapse = " "), if (remainder > 0L) paste0(" Plus ", remainder, " other kind",
    if (remainder > 1L) "s" else "", " of finding.") else "")
}

ScidrReviewVariableSummary <- function(object, display) {
  records <- list()
  catalogs <- stats::setNames(lapply(names(object$Sources), function(analysis) ScidrReviewCatalog(object, analysis)), names(object$Sources))
  for (i in seq_len(nrow(display$Evidence))) {
    e <- display$Evidence[i, ]
    source <- object$Sources[[e$Analysis]]
    stratum <- source$Results$Stratum[match(paste(e$SourceRow, e$Method), paste(source$Results$.Row, source$Results$Method))]
    fields <- unique(stats::na.omit(c(e$Variable, e$OtherVariable)))
    for (variable in fields) {
      catalog <- catalogs[[e$Analysis]]
      rows <- c(if (identical(variable, e$Variable)) e$SourceRow, if (identical(variable, e$OtherVariable)) e$OtherRow)
      records[[length(records) + 1L]] <- ScidrReviewOneRow(Analysis = e$Analysis, Stratum = stratum,
        Variable = variable, Label = catalog$Label[match(variable, catalog$Variable)],
        CodebookCategory = catalog$Category[match(variable, catalog$Variable)],
        IssueID = e$IssueID, EvidenceID = e$EvidenceID, SourceRows = list(unique(stats::na.omit(rows))),
        ReviewTier = e$ReviewTier, ScorePercentile = e$ScorePercentile, ReferenceN = e$ReferenceN,
        ExtremeRatio = e$ExtremeRatio,
        ReviewCategory = e$FindingType,
        StrongestFinding = e$WhyFlagged, Missingness = e$Method == "missingness")
    }
  }
  empty <- tibble::tibble(Analysis = character(), Stratum = integer(), Variable = character(), Label = character(),
    CodebookCategory = character(), StratumLabel = character(), ReviewCategory = character(), AffectedRecords = integer(), StrongestFinding = character(), CommonReasons = character(), MiniDistribution = character(),
    IssueIDs = character(), ReviewTier = integer(), ScorePercentile = double(), ReferenceN = double(),
    ExtremeRatio = double(), ExtremeValue = character(), SuggestedReviewOrder = integer())
  if (!length(records)) return(list(Table = empty, Previews = list()))
  evidence <- dplyr::bind_rows(records)
  groups <- split(seq_len(nrow(evidence)), interaction(evidence$Analysis, evidence$Stratum, evidence$Variable, drop = TRUE, lex.order = TRUE))
  previews <- list()
  tables <- lapply(groups, function(indices) {
    current <- evidence[indices, ]
    chosen <- order(is.na(current$ExtremeRatio), -current$ExtremeRatio, current$ReviewTier, -current$ScorePercentile, -current$ReferenceN, current$EvidenceID, na.last = TRUE)[1L]
    best <- current[chosen, ]
    row <- best[c("Analysis", "Stratum", "Variable", "Label", "CodebookCategory", "ReviewCategory", "ReviewTier", "ScorePercentile", "ReferenceN", "StrongestFinding")]
    row$ExtremeRatio <- best$ExtremeRatio
    row$ExtremeValue <- if (is.finite(best$ExtremeRatio)) "Extreme" else ""
    affected <- unique(unlist(current$SourceRows, use.names = FALSE))
    source <- object$Sources[[best$Analysis]]
    group_row <- match(affected[1], source$ReviewData$.Row)
    row$StratumLabel <- if (length(source$Metadata$StrataVariables)) paste(vapply(source$Metadata$StrataVariables, function(v)
      paste0(catalogs[[best$Analysis]]$Label[match(v, catalogs[[best$Analysis]]$Variable)], ": ", source$ReviewData[[v]][group_row]), character(1)), collapse = "; ") else "All records"
    row$AffectedRecords <- length(affected)
    row$CommonReasons <- ScidrReviewCommonReasons(current$ReviewCategory, current$StrongestFinding)
    row$IssueIDs <- as.character(jsonlite::toJSON(unique(current$IssueID), auto_unbox = FALSE))
    key <- as.character(length(previews) + 1L)
    row$MiniDistribution <- key
    previews[[key]] <<- ScidrReviewMiniDistribution(object$Sources[[best$Analysis]], best$Variable, best$Stratum, affected,
      any(current$Missingness))
    row
  })
  table <- dplyr::bind_rows(tables)
  table <- table[order(is.na(table$ExtremeRatio), -table$ExtremeRatio, table$ReviewTier, -table$ScorePercentile, -table$ReferenceN, -table$AffectedRecords, table$Variable, table$Analysis, table$Stratum, na.last = TRUE), ]
  table$SuggestedReviewOrder <- seq_len(nrow(table))
  list(Table = table, Previews = previews)
}

ScidrReviewMiniDistribution <- function(source, variable, stratum, affected, missingness = FALSE) {
  rows <- unique(source$Results$.Row[source$Results$Stratum %in% stratum])
  index <- match(rows, source$ReviewData$.Row)
  x <- source$ReviewData[[variable]][index]
  highlight <- rows %in% affected
  if (missingness && any(!is.na(x[highlight]))) {
    preview <- ScidrReviewMiniDistribution(source, variable, stratum, affected, FALSE)
    preview$completion <- c(sum(!is.na(x)), sum(is.na(x)))
    return(preview)
  }
  if (missingness) return(list(kind = "Missingness", labels = c("Observed", "Missing"), counts = c(sum(!is.na(x)), sum(is.na(x))), highlighted = c(FALSE, TRUE)))
  temporal <- source$Metadata$TimeConversions
  kind <- temporal$TimeType[match(variable, temporal$Variable)]
  if (length(kind) && !is.na(kind)) x <- tryCatch(ScidrTimeSeconds(x, temporal[match(variable, temporal$Variable), ], strict = FALSE), error = function(e) rep(NA_real_, length(x)))
  clock <- length(kind) && !is.na(kind) && kind == "clock"
  if (clock) x[!is.finite(x) | x < 0 | x >= 86400] <- NA_real_
  if (is.numeric(x) || ScidrAnomalyIsDate(x)) {
    is_date <- inherits(x, "Date")
    is_stamp <- inherits(x, "POSIXct")
    zone <- attr(x, "tzone")
    Format <- function(values) {
      if (clock) ScidrTimeFormat(values, "clock")
      else if (length(kind) && !is.na(kind) && kind == "duration") ScidrTimeFormat(values, "duration")
      else if (is_date) as.character(as.Date(values, origin = "1970-01-01"))
      else if (is_stamp) format(as.POSIXct(values, origin = "1970-01-01", tz = if (length(zone)) zone[1] else "UTC"))
      else format(signif(values, 4), trim = TRUE)
    }
    x <- as.numeric(x)
    valid <- is.finite(x)
    if (!any(valid)) return(list(kind = "Unavailable", labels = character(), counts = numeric(), highlighted = logical()))
    breaks <- if (clock) seq(0, 86400, length.out = 25L) else {
      limits <- range(x[valid]); if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
      seq(limits[1], limits[2], length.out = 25L)
    }
    bins <- cut(x, breaks, include.lowest = TRUE, labels = FALSE)
    counts <- tabulate(bins, nbins = 24L)
    flagged <- seq_len(24L) %in% bins[highlight & valid]
    return(list(kind = if (clock) "Clock" else "Distribution", labels = paste(Format(head(breaks, -1)), "to", Format(tail(breaks, -1))), counts = counts, highlighted = flagged,
      range = if (clock) "00:00\u201324:00" else paste(ScidrTimeObserved(source$ReviewData, variable, index[which(valid)[c(which.min(x[valid]), which.max(x[valid]))]], source$Metadata$TimeConversions), collapse = " to ")))
  }
  counts <- sort(table(as.character(x), useNA = "no"), decreasing = TRUE)
  labels <- names(counts)
  if (length(counts) > 24L) {
    important <- which(labels %in% as.character(x[highlight & !is.na(x)]))
    keep <- head(unique(c(head(important, 12L), seq_along(counts))), 23L)
    rest <- setdiff(seq_along(counts), keep)
    return(list(kind = "Categories", labels = c(labels[keep], "Other categories (combined)"),
      counts = c(as.numeric(counts[keep]), sum(counts[rest])),
      highlighted = c(keep %in% important, any(rest %in% important))))
  }
  list(kind = "Categories", labels = labels, counts = as.numeric(counts), highlighted = labels %in% as.character(x[highlight & !is.na(x)]))
}

ScidrReviewMiniJS <- function() r"---{function(key, previews) {
  const p=previews[key]; if(!p || !p.counts.length) return "No plotting data";
  const esc=function(x) {return String(x).replace(/[&<>"']/g,function(c) {return {"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c];});};
  const counts=p.counts, labels=p.labels, flags=p.highlighted, max=Math.max.apply(null,counts) || 1;
  let body="", width=160;
  counts.forEach(function(n,i) {
    const color=flags[i] ? "#bc3f36" : "#7299ac";
    const title=esc(labels[i]+": "+n+(flags[i] ? " (review values)" : ""));
    if(p.kind === "Clock") {
      const a=2*Math.PI*i/counts.length-Math.PI/2, r=7+11*n/max;
      body+=`<line x1="${80+7*Math.cos(a)}" y1="${21+7*Math.sin(a)}" x2="${80+r*Math.cos(a)}" y2="${21+r*Math.sin(a)}" stroke="${color}" stroke-width="3"><title>${title}</title></line>`;
    } else {
      const w=width/counts.length, h=34*n/max;
      body+=`<rect x="${i*w}" y="${36-h}" width="${Math.max(0.2,w-1)}" height="${h}" fill="${color}"><title>${title}</title></rect>`;
    }
  });
  if(p.completion) {
    const total=p.completion[0]+p.completion[1], missing=160*p.completion[1]/total;
    body+=`<rect x="0" y="39" width="${160-missing}" height="3" fill="#7299ac"><title>Observed: ${p.completion[0]}</title></rect><rect x="${160-missing}" y="39" width="${missing}" height="3" fill="#bc3f36"><title>Missing: ${p.completion[1]}</title></rect>`;
  }
  const title=esc(p.kind+"; red marks review values. "+(p.range || ""));
  return `<svg role="img" aria-label="${title}" viewBox="0 0 160 42" width="160" height="42"><title>${title}</title>${body}</svg>`;
}}---"
