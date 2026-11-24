# Derived presentation fields never participate in queue or snapshot fingerprints.
ScidrReviewOrder <- function(object, display) {
  score_cache <- new.env(parent = emptyenv())
  df_Evidence <- display$Evidence
  df_Severity <- ScidrAnomalySeverity()
  df_Evidence$ReviewTier <- integer(nrow(df_Evidence))
  df_Evidence$WeakSupport <- logical(nrow(df_Evidence))
  df_Evidence$ScorePercentile <- rep(NA_real_, nrow(df_Evidence))
  df_Evidence$FindingType <- df_Evidence$SeverityKey <- character(nrow(df_Evidence))
  df_Evidence$ExtremeRatio <- rep(NA_real_, nrow(df_Evidence))
  for (i in seq_len(nrow(df_Evidence))) {
    e <- df_Evidence[i, ]
    source <- object$Sources[[e$Analysis]]
    field <- !is.na(e$Variable) && nzchar(e$Variable)
    coded <- FALSE
    if (field && e$Method != "missingness" && !e$Variable %in% source$Metadata$TimeConversions$Variable && is.numeric(source$ReviewData[[e$Variable]]) && !ScidrAnomalyIsDate(source$ReviewData[[e$Variable]])) {
      values <- unique(stats::na.omit(source$ReviewData[[e$Variable]]))
      coded <- length(values) > 0L && length(values) <= 3L && all(values %in% c(-1, 0, 1))
    }
    # Extreme magnitude is a property of the recorded value against its own
    # column, so it is read from the retained data rather than from the method
    # that happened to report it. Missingness has no value to judge.
    if (field && e$Method != "missingness" && !e$Variable %in% source$Metadata$TimeConversions$Variable &&
        is.numeric(source$ReviewData[[e$Variable]]) && !ScidrAnomalyIsDate(source$ReviewData[[e$Variable]])) {
      ratio <- ScidrAnomalyExtremeRatio(source$ReviewData[[e$Variable]], match(e$SourceRow, source$ReviewData$.Row))
      if (is.finite(ratio) && ratio >= ScidrAnomalyExtremeFactor()) df_Evidence$ExtremeRatio[i] <- ratio
    }
    conditioned <- length(e$Conditions[[1]]) > 0L &&
      any(vapply(e$Conditions[[1]], function(condition) !identical(condition$column, e$Variable), logical(1)))
    key <- ScidrAnomalyTypeKey(e$Method, e$Kind, field, conditioned, coded)
    df_Evidence$FindingType[i] <- ScidrAnomalyTypeFromKey(key)
    df_Evidence$ReviewTier[i] <- df_Severity$Severity[match(key, df_Severity$Key)]
    df_Evidence$SeverityKey[i] <- key
    # Thin reference support weakens a finding but never promotes it past a more
    # serious kind, so it sorts within the type rather than ahead of it. A
    # record-level finding is not weak for lacking a field; it is weak only when
    # the comparison it rests on is thin or absent.
    df_Evidence$WeakSupport[i] <- !is.finite(e$ReferenceN) || e$ReferenceN < 20
    if (is.finite(e$Score)) {
      table_name <- switch(e$Method, missingness = "CellResults", dates = "DateResults", cellwise = "CellResults", trajectory = "CellResults",
        change = "IntervalResults", transition = "IntervalResults", "Results")
      stratum <- e$NativeEvidence[[1]]$Stratum
      if (is.null(stratum)) stratum <- source$Results$Stratum[match(paste(e$SourceRow, e$Method), paste(source$Results$.Row, source$Results$Method))]
      cache_key <- paste(e$Analysis, e$Method, paste(stratum, collapse = ","), sep = "\r")
      scores <- score_cache[[cache_key]]
      if (is.null(scores)) {
        reference <- source[[table_name]]
        scores <- numeric()
        if (!is.null(reference) && all(c("Method", "Score") %in% names(reference))) {
          reference <- reference[reference$Method == e$Method, , drop = FALSE]
          if (length(stratum) == 1L && !is.na(stratum) && "Stratum" %in% names(reference)) reference <- reference[reference$Stratum == stratum & !is.na(reference$Stratum), , drop = FALSE]
          scores <- sort(reference$Score[is.finite(reference$Score)])
        }
        score_cache[[cache_key]] <- scores
      }
      if (length(scores)) df_Evidence$ScorePercentile[i] <-
        (findInterval(e$Score, scores, left.open = TRUE) + findInterval(e$Score, scores)) / (2 * length(scores))
    }
  }
  df_Issues <- display$Issues
  df_Issues$ReviewTier <- integer(nrow(df_Issues))
  df_Issues$ScorePercentile <- rep(NA_real_, nrow(df_Issues))
  df_Issues$ReviewReferenceN <- rep(NA_real_, nrow(df_Issues))
  df_Issues$FindingType <- df_Issues$OrderRationale <- df_Issues$Methods <- character(nrow(df_Issues))
  df_Issues$WeakSupport <- logical(nrow(df_Issues))
  df_Issues$ExtremeRatio <- rep(NA_real_, nrow(df_Issues))
  df_Issues$ExtremeValue <- character(nrow(df_Issues))
  for (i in seq_len(nrow(df_Issues))) {
    evidence <- df_Evidence[df_Evidence$IssueID == df_Issues$IssueID[i], ]
    chosen <- order(is.na(evidence$ExtremeRatio), -evidence$ExtremeRatio, evidence$ReviewTier, evidence$WeakSupport, -evidence$ScorePercentile, -evidence$ReferenceN, evidence$EvidenceID, na.last = TRUE)[1]
    e <- evidence[chosen, ]
    df_Issues$ReviewTier[i] <- e$ReviewTier
    df_Issues$ExtremeRatio[i] <- e$ExtremeRatio
    df_Issues$ExtremeValue[i] <- if (is.finite(e$ExtremeRatio)) "Extreme" else ""
    df_Issues$WeakSupport[i] <- e$WeakSupport
    df_Issues$ScorePercentile[i] <- e$ScorePercentile
    df_Issues$ReviewReferenceN[i] <- e$ReferenceN
    df_Issues$FindingType[i] <- df_Issues$ReviewCategory[i] <- e$FindingType
    df_Issues$Methods[i] <- paste(sort(unique(evidence$Method)), collapse = ", ")
    df_Issues$OrderRationale[i] <- paste0(
      if (is.finite(e$ExtremeRatio)) paste0("Extreme magnitude: this value sits ",
        signif(e$ExtremeRatio, 3), " times further from the column median than any other record does, so it is listed ahead of the tiers. ") else "",
      "Tier ", e$ReviewTier, ": ",
      df_Severity$Rationale[match(e$SeverityKey, df_Severity$Key)], ". ",
      if (isTRUE(e$WeakSupport)) "Thin reference support. " else "",
      if (is.finite(e$ScorePercentile)) paste0("Within-method/stratum score percentile: ", round(100 * e$ScorePercentile, 1), "%. ") else "Comparable score unavailable. ",
      if (is.finite(e$ReferenceN)) paste0("Reference n: ", e$ReferenceN, ". ") else "Reference count unavailable. ",
      "Sorted by extreme magnitude first, then kind of finding, then reference support, then descending available percentile within that method, then reference count and issue ID. Scores are not comparable across methods. Unvalidated review aid, not error probability or clinical urgency.")
  }
  ordering <- order(is.na(df_Issues$ExtremeRatio), -df_Issues$ExtremeRatio, df_Issues$ReviewTier, df_Issues$WeakSupport, -df_Issues$ScorePercentile, -df_Issues$ReviewReferenceN, df_Issues$IssueID, na.last = TRUE)
  df_Issues$SuggestedReviewOrder <- match(seq_len(nrow(df_Issues)), ordering)
  df_Issues <- df_Issues[ordering, ]
  df_Evidence$SeverityKey <- NULL
  list(Issues = df_Issues, Evidence = df_Evidence)
}
