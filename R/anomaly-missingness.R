# Missingness is screened as a prediction problem: a blank is only a finding when
# comparable records, identified by a learned rule, almost always carry a value.
ScidrAnomalyMissingnessArgs <- function(args) {
  defaults <- list(min_reference = 20L, min_observed_fraction = 0.95,
    cp = 0.01, maxdepth = 4L, max_levels = 20L, max_predictors = 200L,
    min_block_events = 5L, min_block_size = 3L, block_similarity = 0.8, max_block_gap = 0.5)
  args <- ScidrAnomalyArgs(args, names(defaults), defaults)
  for (name in names(args)) if (!is.numeric(args[[name]]) || length(args[[name]]) != 1L ||
      !is.finite(args[[name]])) stop("Missingness settings must be finite numeric scalars.")
  for (name in c("min_reference", "max_levels", "max_predictors", "min_block_events", "min_block_size"))
    if (args[[name]] < 2 || args[[name]] != floor(args[[name]])) stop(name, " must be an integer >= 2.")
  if (args$block_similarity <= 0 || args$block_similarity > 1) stop("block_similarity must be > 0 and <= 1.")
  if (args$max_block_gap <= 0 || args$max_block_gap > 1) stop("max_block_gap must be > 0 and <= 1.")
  if (args$maxdepth < 1 || args$maxdepth > 30 || args$maxdepth != floor(args$maxdepth))
    stop("maxdepth must be an integer between 1 and 30.")
  if (args$cp < 0 || args$cp >= 1) stop("cp must be >= 0 and < 1.")
  if (args$min_observed_fraction <= 0.5 || args$min_observed_fraction > 1)
    stop("min_observed_fraction must be > 0.5 and <= 1.")
  args
}

# Shared by the detection cells and the queue's short reason so the two cannot drift.
ScidrAnomalyMissingnessReason <- function(n_observed, n, group_text = NULL) {
  Number <- function(x) if (length(x) == 1L && !is.na(x) && is.finite(suppressWarnings(as.numeric(x))))
    format(as.numeric(x), trim = TRUE) else "an unknown number of"
  if (!length(n) || is.na(n) || !is.finite(suppressWarnings(as.numeric(n))))
    return("Blank; no comparable records supported a comparison.")
  complete <- length(n_observed) == 1L && !is.na(n_observed) &&
    isTRUE(all.equal(as.numeric(n_observed), as.numeric(n)))
  if (!length(group_text) || is.na(group_text) || !nzchar(group_text))
    return(if (complete)
      paste0("Blank, although every other record in this analysis group has a value (",
        Number(n), " of ", Number(n), ").") else
      paste0("Blank, although ", Number(n_observed), " of the ", Number(n),
        " other records in this analysis group have a value."))
  if (complete) return(paste0("Blank, although every record with ", group_text,
    " has a value (", Number(n), " of ", Number(n), ")."))
  paste0("Blank, although ", Number(n_observed), " of ", Number(n), " records with ",
    group_text, " have a value.")
}

# Synthetic predictor names keep rpart's printed split strings unambiguous whatever
# the real column names contain.
# Built once per stratum over every selected variable; each fit then drops the two
# columns belonging to its own target. Rebuilding it per target made the step
# quadratic in variable count, which is minutes of work on a 1000-column frame.
ScidrAnomalyMissingnessDesign <- function(setup, model_data, rows, args) {
  design <- list()
  lookup <- list()
  for (variable in setup$variables) {
    x <- model_data[[variable]][rows]
    key <- paste0("m", length(lookup) + 1L)
    design[[key]] <- factor(ifelse(is.na(x), "missing", "observed"), levels = c("missing", "observed"))
    lookup[[key]] <- list(column = variable, indicator = TRUE)
    if (is.numeric(x)) {
      if (length(unique(stats::na.omit(x))) < 2L) next
      key <- paste0("p", length(lookup) + 1L)
      design[[key]] <- as.numeric(x)
      lookup[[key]] <- list(column = variable, indicator = FALSE)
    } else {
      levels <- unique(as.character(stats::na.omit(x)))
      if (length(levels) < 2L) next
      if (length(levels) > args$max_levels) next
      key <- paste0("p", length(lookup) + 1L)
      design[[key]] <- factor(as.character(x), levels = sort(levels))
      lookup[[key]] <- list(column = variable, indicator = FALSE)
    }
  }
  usable <- vapply(design, function(x) length(unique(stats::na.omit(x))) > 1L, logical(1))
  list(Design = design[usable], Lookup = lookup[usable])
}

# rpart prints splits as "p3< 18.5" or "p7=a,b"; synthetic names make this parseable.
ScidrAnomalyMissingnessPath <- function(fit, node, lookup, data, setup, rows, row) {
  steps <- tryCatch(rpart::path.rpart(fit, nodes = node, print.it = FALSE)[[1]], error = function(e) character())
  steps <- steps[steps != "root"]
  conditions <- list()
  for (step in steps) {
    parts <- regmatches(step, regexec("^([pm][0-9]+)(>=|<=|<|>|=)(.*)$", step))[[1]]
    if (length(parts) != 4L) next
    entry <- lookup[[parts[2]]]
    if (is.null(entry)) next
    operator <- parts[3]
    value <- parts[4]
    if (isTRUE(entry$indicator)) {
      levels <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      if (length(levels) != 1L) next
      condition <- list(column = entry$column, comparison = if (levels == "missing") "is missing" else "is observed", value_comp = "")
    } else if (operator == "=") {
      levels <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      condition <- list(column = entry$column,
        comparison = if (length(levels) > 1L) "is one of" else "=",
        value_comp = paste(levels, collapse = ", "))
    } else {
      numeric_value <- suppressWarnings(as.numeric(trimws(value)))
      if (!is.finite(numeric_value)) next
      condition <- list(column = entry$column, comparison = operator, value_comp = numeric_value)
    }
    observed <- ScidrTimeObserved(data, entry$column, row, setup$time_metadata)
    condition$value_this <- if (is.na(observed)) "(missing)" else observed
    conditions[[length(conditions) + 1L]] <- condition
  }
  conditions
}

ScidrAnomalyMissingnessGroupText <- function(conditions, setup) {
  if (!length(conditions)) return("")
  paste(vapply(conditions, function(condition) {
    label <- unname(setup$labels[condition$column])
    if (!length(label) || is.na(label) || !nzchar(label)) label <- condition$column
    value <- condition$value_comp
    if (is.numeric(value)) {
      conversion <- setup$date_conversions[[condition$column]]
      kind <- setup$time_metadata$TimeType[match(condition$column, setup$time_metadata$Variable)]
      value <- if (!is.null(conversion)) ScidrAnomalyFormatDate(value, conversion)
        else if (length(kind) && !is.na(kind) && kind %in% c("clock", "duration")) ScidrTimeFormat(value, kind)
        else format(signif(value, 6), trim = TRUE)
    }
    trimws(paste(label, condition$comparison, value))
  }, character(1)), collapse = " and ")
}

ScidrAnomalyMissingness <- function(setup, data, model_data, args) {
  timer_started <- proc.time()[[3]]
  result <- ScidrAnomalyRecords(nrow(data), "missingness", setup$groups, NULL)
  cells <- list(ScidrAnomalyEmptyCells())
  block_records <- list()
  message("Screening unexpected missingness against learned comparison groups...")
  cutoff <- if ("missingness" %in% names(setup$thresholds)) setup$thresholds[["missingness"]] else args$min_observed_fraction
  if (cutoff < 0 || cutoff > 1) stop("The missingness cutoff must be between 0 and 1.")
  available <- requireNamespace("rpart", quietly = TRUE)
  if (!available) setup$Add("missingness", NA_integer_, "missingness_model_unavailable",
    "rpart is not installed; missingness was compared against whole analysis groups only.")
  for (stratum in sort(unique(stats::na.omit(setup$groups)))) {
    rows <- which(setup$groups == stratum)
    missing <- vapply(data[rows, setup$variables, drop = FALSE], is.na, logical(length(rows)))
    if (is.null(dim(missing))) missing <- matrix(missing, ncol = length(setup$variables))
    colnames(missing) <- setup$variables
    targets <- which(colSums(missing) > 0 & colSums(missing) < length(rows))
    result$Status[rows] <- "not applicable"
    all_missing <- setup$variables[colSums(missing) == length(rows)]
    if (length(all_missing)) setup$Add("missingness", stratum, "missingness_all_missing", paste(all_missing, collapse = ", "), rows)
    if (!length(targets)) {
      setup$Add("missingness", stratum, "no_mixed_missingness", "No selected variable has both observed and missing values.", rows)
      next
    }
    wide <- setup$variables[vapply(setup$variables, function(variable) {
      x <- model_data[[variable]][rows]
      !is.numeric(x) && length(unique(as.character(stats::na.omit(x)))) > args$max_levels
    }, logical(1))]
    if (length(wide)) setup$Add("missingness", stratum, "missingness_high_cardinality",
      paste(paste(wide, collapse = ", "), "value splits skipped; missing/observed indicators retained."))
    blocks <- ScidrAnomalyMissingnessBlocks(missing, args)
    block_of <- stats::setNames(rep(NA_integer_, length(setup$variables)), setup$variables)
    for (b in seq_along(blocks)) block_of[blocks[[b]]] <- b
    coverage <- lapply(blocks, function(block) ScidrAnomalyMissingnessCoverage(missing, block))
    if (length(blocks)) setup$Add("missingness", stratum, "missingness_blocks_found",
      paste0(length(blocks), " group(s) of variables go missing together and are used as comparison groups: ",
        paste(vapply(blocks, function(block) paste0("(", paste(block, collapse = ", "), ")"), character(1)), collapse = "; ")), rows)
    block_records[[length(block_records) + 1L]] <- if (!length(blocks)) NULL else dplyr::bind_rows(lapply(seq_along(blocks), function(b) {
      started <- sum(coverage[[b]]$Started)
      tibble::tibble(Stratum = stratum, Block = b, Size = length(blocks[[b]]),
        Members = paste(blocks[[b]], collapse = ", "),
        Labels = ScidrAnomalyMissingnessBlockText(blocks[[b]], setup$labels, limit = length(blocks[[b]])),
        Started = started, Completed = sum(coverage[[b]]$Complete),
        CompletionRate = if (started) sum(coverage[[b]]$Complete) / started else NA_real_)
    }))
    ids <- if (length(setup$id_var)) as.character(data[[setup$id_var]][rows]) else rep(NA_character_, length(rows))
    peers_by_row <- lapply(seq_along(rows), function(i) if (!is.na(ids[i])) which(!is.na(ids) & ids == ids[i]) else i)
    stratum_design <- ScidrAnomalyMissingnessDesign(setup, model_data, rows, args)
    stratum_frame <- as.data.frame(stratum_design$Design, stringsAsFactors = FALSE)
    stratum_columns <- vapply(stratum_design$Lookup, function(entry) entry$column, character(1))
    # Ranking matrix for the predictor screen below: factors enter as integer codes,
    # which is crude but only orders candidates -- rpart still does the real split
    # search on whichever survive.
    stratum_numeric <- if (length(stratum_design$Design) > args$max_predictors)
      vapply(stratum_design$Design, function(x) as.numeric(x), numeric(length(rows))) else NULL
    if (!is.null(stratum_numeric)) setup$Add("missingness", stratum, "missingness_predictors_screened",
      paste0(length(stratum_design$Design), " candidate predictors exceed max_predictors (", args$max_predictors,
        "); the strongest by marginal association were kept for each variable."), rows)
    report_every <- max(1L, as.integer(ceiling(length(targets) / 20)))
    for (position in seq_along(targets)) {
      target_index <- targets[position]
      variable <- setup$variables[target_index]
      if (length(targets) > 50L && (position %% report_every == 0L || position == length(targets)))
        message("Screening missingness: variable ", position, "/", length(targets), ".")
      blank <- missing[, target_index]
      keep <- stratum_columns != variable
      if (!is.null(stratum_numeric) && sum(keep) > args$max_predictors) {
        strength <- abs(suppressWarnings(stats::cor(as.numeric(blank),
          stratum_numeric[, keep, drop = FALSE], use = "pairwise.complete.obs")[1, ]))
        strength[!is.finite(strength)] <- 0
        cut_at <- sort(strength, decreasing = TRUE)[args$max_predictors]
        keep[keep][strength < cut_at] <- FALSE
      }
      design <- list(Design = stratum_design$Design[keep], Lookup = stratum_design$Lookup[keep])
      leaves <- rep(1L, length(rows))
      fit <- NULL
      block <- block_of[[variable]]
      if (available && length(design$Design)) {
        frame <- stratum_frame[, names(design$Design), drop = FALSE]
        frame$.missing <- factor(ifelse(blank, "missing", "observed"), levels = c("missing", "observed"))
        fit <- tryCatch(rpart::rpart(.missing ~ ., data = frame, method = "class",
          control = rpart::rpart.control(minbucket = args$min_reference, minsplit = 2L * args$min_reference,
            cp = args$cp, maxdepth = args$maxdepth, xval = 0L)), error = function(e) NULL)
        if (!is.null(fit) && nrow(fit$frame) > 1L) leaves <- as.integer(fit$where) else fit <- NULL
      }
      if (is.null(fit) && is.na(block)) setup$Add("missingness", stratum, "missingness_unexplained",
        paste(variable, "blanks are explained neither by another column's values nor by a group of variables that go missing together; no comparison group, so no finding."),
        rows[blank])
      # No learned rule and no co-missing group: nothing to compare against, so
      # this variable yields no findings. Status stays "not applicable".
      if (is.null(fit) && is.na(block)) next
      node_numbers <- if (is.null(fit)) NULL else as.integer(rownames(fit$frame))
      blank_rows <- which(blank)
      reference_n <- reference_observed <- scores <- rep(NA_real_, length(blank_rows))
      flags <- rep(NA, length(blank_rows))
      cell_conditions <- vector("list", length(blank_rows))
      cell_statistics <- vector("list", length(blank_rows))
      explanations <- character(length(blank_rows))
      for (k in seq_along(blank_rows)) {
        local_row <- blank_rows[k]
        peers <- peers_by_row[[local_row]]
        conditions <- list()
        group_text <- ""
        block_statistics <- NULL
        tree_conditions <- if (is.null(fit)) list() else
          ScidrAnomalyMissingnessPath(fit, node_numbers[leaves[local_row]], design$Lookup, data, setup, rows, rows[local_row])
        tree_n <- tree_observed <- 0L
        if (!is.null(fit)) {
          member <- leaves == leaves[local_row]
          member[peers] <- FALSE
          tree_n <- sum(member)
          tree_observed <- sum(member & !blank)
        }
        # Two things can explain a blank, and they are not equally informative.
        # A rule on another column's *values* is a skip rule the record appears to
        # break, so it wins -- but only when it actually explains this blank. A
        # tree that merely splits on other columns' missing/observed indicators has
        # rediscovered the co-missing group one split at a time, and the group says
        # the same thing while naming every member, which is the pattern worth
        # printing. A spurious value split (a sequential date separating the
        # skippers, say) explains nothing and must not displace a real group.
        value_rule <- length(tree_conditions) > 0L &&
          any(vapply(tree_conditions, function(condition) !condition$comparison %in% c("is missing", "is observed"), logical(1)))
        tree_explains <- tree_n >= args$min_reference && tree_observed / tree_n >= cutoff
        if (!is.null(fit) && (is.na(block) || (value_rule && tree_explains))) {
          n <- tree_n
          n_observed <- tree_observed
          conditions <- tree_conditions
          group_text <- ScidrAnomalyMissingnessGroupText(conditions, setup)
        } else {
          # Structural: this variable belongs to a group that goes missing
          # together, so the comparison group is the records that carry the
          # group at all, and the question is whether they finish it.
          cover <- coverage[[block]]
          started <- cover$Started
          started[peers] <- FALSE
          n <- sum(started)
          n_observed <- sum(started & cover$Complete)
          gap <- cover$Size - cover$Observed[local_row]
          # Missing most of the group is a skipped form, not a gap in a form.
          if (gap > args$max_block_gap * cover$Size) n <- 0L
          block_statistics <- list(block = block, block_size = cover$Size,
            block_members = blocks[[block]], observed_in_block = cover$Observed[local_row], block_gap = gap)
        }
        supported <- n >= args$min_reference
        score <- if (supported) n_observed / n else NA_real_
        flag <- if (supported) score >= cutoff else NA
        row <- rows[local_row]
        reference_n[k] <- if (supported) n else NA_real_
        reference_observed[k] <- if (supported) n_observed else NA_real_
        scores[k] <- score
        flags[k] <- flag
        cell_conditions[[k]] <- conditions
        cell_statistics[[k]] <- c(list(n_obs = if (supported) n else NA_real_,
          n_observed = if (supported) n_observed else NA_real_,
          n_missing = if (supported) n - n_observed else NA_real_, observed_fraction = score,
          comparison_key = if (nzchar(group_text)) group_text else if (!is.null(block_statistics))
            paste0("co-missing group (", paste(block_statistics$block_members, collapse = ", "), ")") else "unexplained",
          comparison_text = group_text), block_statistics)
        # A tree leaf always carries a path; an empty one means the split could not
        # be parsed back, and claiming "every other record" would overstate the group.
        if (is.null(block_statistics) && !nzchar(group_text)) group_text <- ""
        explanations[k] <- if (!supported) (if (!is.null(block_statistics) && block_statistics$block_gap > args$max_block_gap * block_statistics$block_size)
            "Blank, together with most of the group of variables it belongs to; this reads as a form that was not done, not a gap within one." else
            "Blank; too few comparable records to support a comparison.") else
          if (is.null(block_statistics)) ScidrAnomalyMissingnessReason(n_observed, n, group_text) else
          ScidrAnomalyMissingnessBlockReason(block_statistics$observed_in_block, block_statistics$block_size,
            n_observed, n, block_statistics$block_members, setup$labels)
        if (result$Status[row] != "evaluated") result$Status[row] <- if (supported) "evaluated" else "unavailable"
        if (supported) {
          result$Score[row] <- if (is.na(result$Score[row])) score else max(result$Score[row], score)
          result$Flag[row] <- isTRUE(result$Flag[row]) || flag
          result$Threshold[row] <- cutoff
          result$Explanation[row] <- "One or more fields are blank where comparable records carry a value; inspect variable-level evidence."
        }
      }
      cells[[length(cells) + 1L]] <- tibble::tibble(.Row = rows[blank_rows], Method = "missingness", Stratum = stratum,
        Variable = variable, Label = unname(setup$labels[variable]), Observed = NA_character_, Expected = NA_real_,
        Residual = NA_real_, Score = scores, Flag = flags, Explanation = explanations, Conditions = cell_conditions,
        Status = ifelse(is.finite(scores), "evaluated", "unavailable"), GroupStatistics = cell_statistics)
    }
    unsupported <- sum(result$Status[rows] == "unavailable")
    if (unsupported) setup$Add("missingness", stratum, "insufficient_missingness_support", paste(unsupported, "records lack sufficient reference support."), rows[result$Status[rows] == "unavailable"])
  }
  elapsed <- unname(proc.time()[[3]] - timer_started)
  message("Missingness screening completed in ", round(elapsed, 2), " seconds.")
  list(Results = result, Cells = dplyr::bind_rows(cells), ElapsedSeconds = elapsed,
    Blocks = dplyr::bind_rows(block_records))
}
