# Internal date adapters. Model copies use elapsed days; originals stay untouched.
ScidrAnomalyIsDate <- function(x) inherits(x, "Date") || inherits(x, "POSIXct")

ScidrAnomalyDateDays <- function(x) {
  if (inherits(x, "POSIXct")) as.numeric(x) / 86400 else as.numeric(x)
}

ScidrAnomalyFormatDate <- function(days, conversion) {
  absolute <- days + conversion$ReferenceDays
  absolute[!is.finite(absolute)] <- NA_real_
  if ("Date" %in% conversion$InputClass) return(as.character(as.Date(absolute, origin = "1970-01-01")))
  zone <- conversion$TimeZone
  if (!length(zone) || is.na(zone[1]) || !nzchar(zone[1])) zone <- "UTC"
  format(as.POSIXct(absolute * 86400, origin = "1970-01-01", tz = zone[1]),
    format = "%Y-%m-%d %H:%M:%OS6 %Z", tz = zone[1])
}

ScidrAnomalyDateInterval <- function(first, second) {
  # Subtract seconds before dividing to avoid cancellation in nearly constant
  # sub-day intervals across large calendar-date values. Date values mean UTC midnight.
  first_seconds <- if (inherits(first, "Date")) as.numeric(first) * 86400 else as.numeric(first)
  second_seconds <- if (inherits(second, "Date")) as.numeric(second) * 86400 else as.numeric(second)
  (second_seconds - first_seconds) / 86400
}

ScidrAnomalyEmptyDates <- function() tibble::tibble(.Row = integer(), Method = character(),
  Stratum = integer(), CheckType = character(), Variable = character(), Label = character(),
  OtherVariable = character(), OtherLabel = character(), Observed = character(),
  OtherObserved = character(), ValueDays = double(), N = integer(),
  MinimumDays = double(), MaximumDays = double(), MinimumDate = character(), MaximumDate = character(), LowerFenceDays = double(),
  UpperFenceDays = double(), LowerDate = character(), UpperDate = character(),
  Rule = character(), Status = character(), Score = double(), Flag = logical(), Explanation = character())

ScidrAnomalyDateScreen <- function(setup, data, model_data, variables, pairs, args) {
  df_Result <- ScidrAnomalyRecords(nrow(data), "dates", setup$groups, NULL)
  evidence <- list(ScidrAnomalyEmptyDates())
  if (!length(variables)) df_Result$Status[!is.na(setup$groups)] <- "not applicable"
  # These are per-column fences and learned interval distributions, not clinical rules.
  for (stratum in sort(unique(stats::na.omit(setup$groups)))) {
    rows <- which(setup$groups == stratum)
    checks <- c(lapply(variables, function(v) v),
      if (nrow(pairs)) lapply(seq_len(nrow(pairs)), function(i) unname(pairs[i, ])))
    for (check in checks) {
      variable <- check[1]
      other <- if (length(check) == 2L) check[2] else NA_character_
      pair <- !is.na(other)
      label <- setup$labels[[variable]]
      other_label <- if (pair) setup$labels[[other]] else NA_character_
      description <- if (pair) paste(other_label, "minus", label, "(days)") else label
      values <- if (pair) ScidrAnomalyDateInterval(data[[variable]][rows], data[[other]][rows]) else model_data[[variable]][rows]
      usable <- is.finite(values)
      observed <- values[usable]
      n <- length(observed)
      lower <- upper <- NA_real_
      rule <- "unavailable"
      reason <- "Insufficient finite observations."
      if (n >= args$min_observations) {
        quartiles <- stats::quantile(observed, c(0.25, 0.75), names = FALSE, type = 7)
        spread <- diff(quartiles)
        if (is.finite(spread) && spread > 0) {
          lower <- quartiles[1] - args$iqr_multiplier * spread
          upper <- quartiles[2] + args$iqr_multiplier * spread
          rule <- "outer_tukey"
        } else {
          unique_values <- unique(observed)
          counts <- tabulate(match(observed, unique_values))
          mode <- unique_values[which.max(counts)]
          if (max(counts) / n >= args$majority_fraction) {
            lower <- upper <- mode
            rule <- "learned_majority"
            spread <- 0
          } else reason <- "Zero IQR without sufficient agreement on a single value."
        }
      }
      evaluated <- usable & rule != "unavailable"
      flags <- rep(NA, length(rows))
      scores <- rep(NA_real_, length(rows))
      if (rule != "unavailable") {
        flags[evaluated] <- values[evaluated] < lower | values[evaluated] > upper
        scores[evaluated] <- pmax(lower - values[evaluated], values[evaluated] - upper, 0) / max(spread, 1)
      } else setup$Add("dates", stratum, "date_check_unavailable", paste(description, reason), rows)
      if (any(!usable)) setup$Add("dates", stratum, "missing_or_nonfinite_dates", description, rows[!usable])
      explanation <- paste(description, if (rule == "unavailable") reason else if (rule == "learned_majority")
        "Compared with a learned near-constant value; departures warrant review." else
        "Compared with empirical outer Tukey fences; outside values warrant review.")
      evidence[[length(evidence) + 1L]] <- tibble::tibble(.Row = rows, Method = "dates", Stratum = stratum,
        CheckType = if (pair) "interval" else "column", Variable = variable, Label = label,
        OtherVariable = other, OtherLabel = other_label,
        Observed = ScidrAnomalyFormatDate(model_data[[variable]][rows], setup$date_conversions[[variable]]),
        OtherObserved = if (pair) ScidrAnomalyFormatDate(model_data[[other]][rows], setup$date_conversions[[other]]) else NA_character_,
        ValueDays = values, N = as.integer(n), MinimumDays = if (n) min(observed) else NA_real_,
        MaximumDays = if (n) max(observed) else NA_real_,
        MinimumDate = if (!pair && n) ScidrAnomalyFormatDate(min(observed), setup$date_conversions[[variable]]) else NA_character_,
        MaximumDate = if (!pair && n) ScidrAnomalyFormatDate(max(observed), setup$date_conversions[[variable]]) else NA_character_,
        LowerFenceDays = lower, UpperFenceDays = upper,
        LowerDate = if (!pair && is.finite(lower)) ScidrAnomalyFormatDate(lower, setup$date_conversions[[variable]]) else NA_character_,
        UpperDate = if (!pair && is.finite(upper)) ScidrAnomalyFormatDate(upper, setup$date_conversions[[variable]]) else NA_character_,
        Rule = rule, Status = ifelse(evaluated, "evaluated", "unavailable"), Score = scores,
        Flag = flags, Explanation = explanation)
    }
  }
  df_Dates <- dplyr::bind_rows(evidence)
  for (row in unique(df_Dates$.Row[df_Dates$Status == "evaluated"])) {
    current <- df_Dates[df_Dates$.Row == row & df_Dates$Status == "evaluated", , drop = FALSE]
    df_Result$Status[row] <- "evaluated"
    df_Result$Score[row] <- max(current$Score)
    df_Result$Flag[row] <- any(current$Flag)
    df_Result$Explanation[row] <- if (any(current$Flag)) paste(unique(current$Explanation[current$Flag]), collapse = "\n") else "No evaluated date check exceeded its empirical limits."
  }
  list(Results = df_Result, DateResults = df_Dates)
}

ScidrAnomalyDatePairs <- function(date_pairs, variables) {
  if (is.null(date_pairs)) {
    if (choose(length(variables), 2) > 100L) warning("Automatic date selection produced ", choose(length(variables), 2), " pairs; supply date_pairs to limit comparisons.", call. = FALSE)
    pairs <- if (length(variables) >= 2L) t(utils::combn(variables, 2)) else matrix(character(), ncol = 2)
  } else {
    if (!is.matrix(date_pairs) || !is.character(date_pairs) || ncol(date_pairs) != 2L ||
        anyNA(date_pairs) || any(!date_pairs %in% variables)) stop("date_pairs must be a two-column character matrix naming selected Date/POSIXct columns.")
    pairs <- date_pairs
    if (nrow(pairs)) {
      if (any(pairs[, 1] == pairs[, 2])) stop("date_pairs cannot compare a column with itself.")
      # Preserve the selected-column order regardless of matrix row orientation.
      pairs <- t(vapply(seq_len(nrow(pairs)), function(i) variables[sort(match(pairs[i, ], variables))], character(2)))
      if (anyDuplicated(as.data.frame(pairs))) stop("date_pairs contains duplicate unordered pairs.")
    }
  }
  colnames(pairs) <- c("First", "Second")
  pairs
}

# Format numerical model conditions using original date units for human review.
ScidrAnomalyDateExplanation <- function(finding, variable, row, setup, data) {
  FormatValue <- function(value, column) {
    conversion <- setup$date_conversions[[column]]
    if (!is.null(conversion) && is.numeric(value)) ScidrAnomalyFormatDate(value, conversion)
    else paste(value, collapse = ", ")
  }
  observed <- if (ScidrAnomalyIsDate(data[[variable]])) {
    conversion <- setup$date_conversions[[variable]]
    ScidrAnomalyFormatDate(ScidrAnomalyDateDays(data[[variable]][row]) - conversion$ReferenceDays, conversion)
  } else as.character(data[[variable]][row])
  text <- paste(setup$labels[[variable]], "=", observed, "is unusual within the learned subgroup.")
  bounds <- finding$group_statistics
  if (!is.null(bounds$lower_thr)) text <- c(text, paste("Learned lower limit:", FormatValue(bounds$lower_thr, variable)))
  if (!is.null(bounds$upper_thr)) text <- c(text, paste("Learned upper limit:", FormatValue(bounds$upper_thr, variable)))
  for (condition in finding$conditions) {
    column <- condition$column
    text <- c(text, paste("Given", setup$labels[[column]], condition$comparison,
      FormatValue(condition$value_comp, column), "(observed:", FormatValue(condition$value_this, column), ")"))
  }
  paste(text, collapse = "\n")
}
