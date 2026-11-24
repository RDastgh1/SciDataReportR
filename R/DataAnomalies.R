# Internal infrastructure shared by the anomaly detectors.

# A column name can survive an encoding round trip and come back mangled: UTF-8
# read back as MacRoman turns "IFN-b" (the Greek beta) into "IFN-Œ≤", which is
# what a CSV opened and resaved in Excel on a Mac produces. The mangled name
# matches no column, so the column is skipped and its findings are silently lost.
# The round trip is reversible, so it is undone here and the column the name was
# meant to reach is reported.
ScidrRepairEncoding <- function(x) {
  repairs <- lapply(c("macintosh", "latin1", "CP1252"), function(encoding) {
    bytes <- suppressWarnings(iconv(enc2utf8(x), "UTF-8", encoding, toRaw = TRUE))
    vapply(bytes, function(raw) {
      if (is.null(raw)) return(NA_character_)
      text <- rawToChar(raw)
      Encoding(text) <- "UTF-8"
      if (validUTF8(text)) enc2utf8(text) else NA_character_
    }, character(1), USE.NAMES = FALSE)
  })
  lapply(seq_along(x), function(i) {
    candidates <- unique(vapply(repairs, `[`, character(1), i))
    candidates[!is.na(candidates) & candidates != x[i]]
  })
}

# The mangling can sit on either side: a stale variables vector carrying mangled
# names, or a data frame read in with mangled headers. Both are checked, and a
# name is only suggested when exactly one column answers to it, since repairing
# a truncated name could otherwise point at the wrong column.
ScidrNameSuggestions <- function(absent, present) {
  if (!length(absent) || !length(present)) return(stats::setNames(rep(NA_character_, length(absent)), absent))
  from_absent <- ScidrRepairEncoding(absent)
  from_present <- ScidrRepairEncoding(present)
  suggestions <- vapply(seq_along(absent), function(i) {
    hits <- unique(c(present[present %in% from_absent[[i]]],
      present[vapply(from_present, function(candidates) absent[i] %in% candidates, logical(1))]))
    if (length(hits) == 1L) hits else NA_character_
  }, character(1))
  stats::setNames(suggestions, absent)
}

ScidrNameReport <- function(absent, present) {
  suggestions <- ScidrNameSuggestions(absent, present)
  paste(ifelse(is.na(suggestions), absent,
    paste0(absent, " (did you mean ", suggestions, "?)")), collapse = ", ")
}

ScidrAnomalySetup <- function(data, variables, id_var, stratify_by, methods,
    allowed, thresholds, method_args, Relabel, seed, time_var = NULL,
    covariates = NULL, allow_dates = FALSE) {
  CheckNames <- function(x, argument, optional = FALSE) {
    if (optional && is.null(x)) return(invisible(NULL))
    if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
      stop(argument, " must be a non-empty, unique character vector.", call. = FALSE)
  }
  if (!is.data.frame(data) || anyDuplicated(names(data))) stop("data must have unique column names.")
  CheckNames(variables, "variables")
  for (arg in c("id_var", "stratify_by", "time_var", "covariates"))
    CheckNames(get(arg), arg, TRUE)
  if (length(id_var) > 1L || length(time_var) > 1L) stop("id_var and time_var must name single columns.")
  CheckNames(methods, "methods")
  if (any(!methods %in% allowed)) stop("Unknown methods: ", paste(setdiff(methods, allowed), collapse = ", "))
  # Context columns are structural: a missing id_var or stratum changes what every
  # finding means, so those still stop. Screened variables and covariates are a
  # selection, so a name the frame does not carry is dropped with a warning -- a
  # stale column in a long variables vector should not cost the whole run.
  context_vars <- unique(c(id_var, time_var, stratify_by))
  if (any(!context_vars %in% names(data)))
    stop("Columns not found: ", paste(setdiff(context_vars, names(data)), collapse = ", "))
  absent_columns <- setdiff(unique(c(variables, covariates)), names(data))
  if (length(absent_columns)) {
    warning("Columns not found and skipped: ", ScidrNameReport(absent_columns, names(data)), call. = FALSE)
    variables <- setdiff(variables, absent_columns)
    covariates <- if (is.null(covariates)) NULL else setdiff(covariates, absent_columns)
    if (!length(covariates)) covariates <- NULL
    if (!length(variables)) stop("No requested variables are present in data.", call. = FALSE)
  }
  all_vars <- unique(c(variables, id_var, time_var, stratify_by, covariates))
  requested_variables <- variables
  if (allow_dates) variables <- setdiff(variables, id_var)
  if (!length(variables)) stop("No model variables remain after excluding id_var.")
  if (length(intersect(variables, c(id_var, time_var, stratify_by)))) stop("Model variables cannot include ID, time, or strata columns.")
  reserved_context <- intersect(c(id_var, time_var, stratify_by),
    c("Method", "Stratum", "Status", "Score", "Flag", "Threshold", "Explanation", "Rank"))
  if (length(reserved_context)) stop("Identifier, time, or strata names conflict with result columns: ", paste(reserved_context, collapse = ", "), ". Rename these columns before detection.")
  if (".Row" %in% names(data)) stop(".Row is reserved for original row references.")
  supported <- vapply(data[unique(c(variables, covariates))], function(x)
    is.null(dim(x)) && ((allow_dates && (ScidrAnomalyIsDate(x) || inherits(x, "difftime"))) || (is.numeric(x) && !inherits(x, c("Date", "POSIXt", "difftime"))) || is.factor(x) || is.character(x) || is.logical(x)), logical(1))
  if (any(!supported)) stop("Unsupported selected types: ", paste(names(supported)[!supported], collapse = ", "))
  if (!is.null(thresholds) && (!is.numeric(thresholds) || is.null(names(thresholds)) ||
      anyNA(thresholds) || any(!is.finite(thresholds)) || anyDuplicated(names(thresholds)) ||
      any(!names(thresholds) %in% setdiff(methods, c("conditional", "dates"))))) stop("thresholds must be finite named numeric cutoffs for requested score-based methods.")
  if (!is.list(method_args) || (length(method_args) && (is.null(names(method_args)) ||
      anyDuplicated(names(method_args)) || any(!names(method_args) %in% methods)))) stop("method_args must be a named list keyed by requested methods.")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) stop("seed must be one non-negative integer.")
  labels <- ScidrDisplayLabels(data, all_vars, Relabel)
  df_Review <- tibble::as_tibble(data[all_vars])
  df_Review$.Row <- seq_len(nrow(data))
  df_Groups <- if (length(stratify_by)) data[stratify_by] else data.frame(.All = rep("All", nrow(data)))
  valid <- stats::complete.cases(df_Groups)
  groups <- rep(NA_integer_, nrow(data))
  if (any(valid)) groups[valid] <- as.integer(dplyr::group_indices(dplyr::group_by(
    tibble::as_tibble(df_Groups[valid, , drop = FALSE]), dplyr::across(dplyr::everything()))))
  diagnostics <- list()
  Add <- function(method = NA_character_, stratum = NA_integer_, issue, detail, rows = integer()) {
    diagnostics[[length(diagnostics) + 1L]] <<- tibble::tibble(Method = method,
      Stratum = stratum, Issue = issue, Detail = detail, Rows = list(as.integer(rows)))
  }
  if (length(absent_columns)) Add(issue = "absent_columns",
    detail = paste0(ScidrNameReport(absent_columns, names(data)), ": not present in data and skipped."))
  if (allow_dates && length(intersect(requested_variables, id_var))) Add(issue = "id_excluded_from_model", detail = id_var)
  for (variable in all_vars) {
    label <- attr(data[[variable]], "label", exact = TRUE)
    if (is.null(label) || !length(label) || is.na(label[1]) || !nzchar(label[1]))
      Add(issue = "missing_label", detail = variable)
  }
  if (any(!valid)) Add(issue = "missing_stratum", detail = "Rows cannot be assigned to a stratum.", rows = which(!valid))
  if (length(id_var)) {
    missing_id <- is.na(data[[id_var]])
    if (any(missing_id)) Add(issue = "missing_id", detail = id_var, rows = which(missing_id))
    duplicates <- !missing_id & (duplicated(data[[id_var]]) | duplicated(data[[id_var]], fromLast = TRUE))
    if (any(duplicates)) Add(issue = "repeated_id", detail = "Repeated participant identifiers.", rows = which(duplicates))
  }
  list(data = data, requested_variables = requested_variables, variables = variables, id_var = id_var, time_var = time_var,
    stratify_by = stratify_by, covariates = covariates, methods = methods,
    thresholds = thresholds, method_args = method_args, seed = seed, labels = labels, Relabel = Relabel,
    ReviewData = df_Review, groups = groups, Add = Add,
    diagnostics = function() dplyr::bind_rows(c(list(tibble::tibble(Method = character(), Stratum = integer(), Issue = character(), Detail = character(), Rows = list())), diagnostics)))
}

ScidrAnomalySeed <- function(seed) {
  exists_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (exists_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (exists_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }
}

ScidrAnomalyArgs <- function(args, allowed, defaults = list()) {
  if (is.null(args)) args <- list()
  if (!is.list(args) || (length(args) && (is.null(names(args)) || anyDuplicated(names(args)) || any(!names(args) %in% allowed))))
    stop("Unknown, unnamed, or wrapper-controlled method_args.", call. = FALSE)
  utils::modifyList(defaults, args)
}

ScidrAnomalyEmptyCells <- function() tibble::tibble(.Row = integer(), Method = character(),
  Stratum = integer(), Variable = character(), Label = character(), Observed = character(),
  Expected = double(), Residual = double(), Score = double(), Flag = logical(),
  Explanation = character(), Conditions = list())

ScidrAnomalyFinish <- function(setup, results, cells, intervals, call, type, settings, date_results = NULL,
    flag_percentile = NULL) {
  df_Results <- dplyr::bind_rows(results)
  df_Cells <- dplyr::bind_rows(c(list(ScidrAnomalyEmptyCells()), cells))
  df_Intervals <- dplyr::bind_rows(c(list(tibble::tibble(PreviousRow = integer(), .Row = integer(),
    Method = character(), Stratum = integer(), Variable = character(), Label = character(),
    PreviousValue = character(), Observed = character(), Elapsed = double(), Score = double(), Flag = logical())), intervals))
  df_Results <- df_Results %>% dplyr::group_by(.data$Method, .data$Stratum) %>%
    dplyr::mutate(Rank = rank(-.data$Score, ties.method = "min", na.last = "keep")) %>% dplyr::ungroup()
  # Score-only methods leave every flag NA; supply an explicit review budget instead of nothing.
  percentile_methods <- character()
  if (!is.null(flag_percentile) && nrow(df_Results)) {
    keys <- unique(df_Results[c("Method", "Stratum")])
    for (i in seq_len(nrow(keys))) {
      rows <- which(df_Results$Method == keys$Method[i] &
        (df_Results$Stratum == keys$Stratum[i] | (is.na(df_Results$Stratum) & is.na(keys$Stratum[i]))))
      scores <- df_Results$Score[rows]
      if (keys$Method[i] == "missingness") next
      if (!all(is.na(df_Results$Flag[rows])) || !any(is.finite(scores))) next
      cutoff <- unname(stats::quantile(scores[is.finite(scores)], flag_percentile, type = 7, names = FALSE))
      df_Results$Threshold[rows] <- cutoff
      df_Results$Flag[rows] <- ifelse(is.finite(scores), scores >= cutoff, NA)
      percentile_methods <- c(percentile_methods, keys$Method[i])
    }
  }
  percentile_methods <- unique(percentile_methods)
  context <- setup$ReviewData[unique(c(".Row", setup$id_var, setup$stratify_by, setup$time_var))]
  df_Results <- dplyr::left_join(df_Results, context, by = ".Row")
  df_Dates <- dplyr::bind_rows(list(ScidrAnomalyEmptyDates(), date_results))
  if (nrow(df_Dates)) df_Dates <- dplyr::left_join(df_Dates, context, by = ".Row")
  summary <- df_Results %>% dplyr::group_by(.data$Method, .data$Stratum) %>% dplyr::summarise(
    Records = dplyr::n(), Evaluated = sum(.data$Status == "evaluated"),
    Unavailable = sum(.data$Status == "unavailable"), NotApplicable = sum(.data$Status == "not applicable"),
    Flagged = sum(.data$Flag %in% TRUE), .groups = "drop")
  participant <- NULL
  if (length(setup$id_var)) participant <- df_Results %>%
    dplyr::filter(!is.na(.data[[setup$id_var]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(setup$id_var))) %>% dplyr::summarise(
      EvaluatedRecords = dplyr::n_distinct(.data$.Row[.data$Status == "evaluated"]),
      FlaggedRecords = dplyr::n_distinct(.data$.Row[.data$Flag %in% TRUE]),
      MethodsFlagged = dplyr::n_distinct(.data$Method[.data$Flag %in% TRUE]), .groups = "drop")
  packages <- c("SciDataReportR", "hms", "outliertree", "isotree", "dbscan", "cluster", "robustbase", "cellWise", "nlme")
  versions <- vapply(packages, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1))
  structure(list(Results = df_Results, CellResults = df_Cells,
    DateResults = df_Dates, IntervalResults = df_Intervals, ParticipantSummary = participant, Summary = summary,
    Diagnostics = setup$diagnostics(), ReviewData = setup$ReviewData,
    Metadata = list(Call = call, RequestedVariables = setup$requested_variables,
      Variables = setup$variables, DateConversions = setup$date_conversions, Labels = setup$labels, Relabel = setup$Relabel,
      IDVariable = setup$id_var, TimeVariable = setup$time_var, StrataVariables = setup$stratify_by,
      Covariates = setup$covariates, Seed = setup$seed, Methods = setup$methods,
      Thresholds = setup$thresholds, Settings = settings, PackageVersions = versions,
      ScoreDefinitions = c(dates = "Maximum excess beyond empirical fences divided by max(IQR, 1 day); native flags", conditional = "Native flags; no general record score",
        missingness = "Observed fraction in the learned comparison group for a blank field, leaving the participant out; not error probability",
        isolation = "Native isolation anomaly score", local = "Local outlier factor",
        distance = "Robust Mahalanobis distance", cellwise = "Absolute standardized DDC residual",
        change = "Absolute MAD-standardized change rate", trajectory = "Absolute in-sample Pearson residual",
        transition = "1 - leave-one-participant-out P(current | previous)"),
      FlagPercentile = flag_percentile, FlagPercentileMethods = percentile_methods,
      ThresholdProvenance = paste0("User cutoffs where supplied; native conditional/DDC/date flags otherwise",
        if (length(percentile_methods)) paste0("; default top ", signif(100 * (1 - flag_percentile), 3),
          "% percentile cutoffs within method and stratum for ", paste(percentile_methods, collapse = ", "))
        else "; remaining flags unavailable"))),
    class = c(type, "DataAnomalyResult", "list"))
}

ScidrAnomalyRecords <- function(n, method, groups, thresholds) {
  cutoff <- if (method %in% names(thresholds)) unname(thresholds[[method]]) else NA_real_
  tibble::tibble(.Row = seq_len(n), Method = method, Stratum = groups,
    Status = "unavailable", Score = NA_real_, Flag = NA, Threshold = cutoff,
    Explanation = NA_character_)
}

#' Detect unusual combinations of recorded values
#'
#' Screen explicitly selected variables without changing the input data.
#' Flags are review candidates, not proof of entry errors. Attached variable
#' labels (for example from [RevalueData()]) are used for display.
#' @param data A data frame with uniquely named columns.
#' @param variables Unique character vector of model input columns. NULL (default)
#'   selects all supported columns except supplied ID and strata columns, reporting
#'   unsupported columns in Diagnostics. Names that \code{data} does not carry are
#'   skipped with a warning and recorded in an \code{absent_columns} diagnostic, so
#'   a stale entry in a long selection does not cost the run; \code{id_var},
#'   \code{time_var} and \code{stratify_by} must exist, since a missing one would
#'   change what every finding means. RevalueData() is not required; absent
#'   codebook labels fall back to column names. Numeric codes
#'   remain numeric; explicitly recode categorical variables before detection.
#' @param id_var Optional participant ID column name; retained for review and
#'   automatically excluded from model variables even if selected.
#' @param stratify_by Optional character vector defining independent fits and ranks.
#' @param methods Subset of `conditional`, `isolation`, `local`, `distance`,
#'   `cellwise`, `dates`. Dates screening needs no optional backend package.
#' @param thresholds Named numeric cutoffs on native scores. Scores at least the
#'   cutoff are flagged. Without cutoffs score-only flags are `NA`. Conditional
#'   and date flags are native and cannot be overridden. Cellwise cutoffs use absolute residuals.
#' @param method_args Named lists of backend settings keyed by method. Only explicit
#'   backend formals are accepted, excluding wrapper-controlled data, seeds,
#'   output formats, weights, and input transformations. For cellwise supply DDC
#'   settings directly (not nested inside `DDCpars`). Dates accepts
#'   `min_observations = 20`, `iqr_multiplier = 3`, `majority_fraction = 0.95`.
#' @param date_pairs NULL examines all unique pairs of selected date columns.
#'   Alternatively supply a two-column character matrix. Pair direction is always
#'   second minus first in selected-variable order, even if a supplied pair is
#'   reversed. Self-pairs and duplicate unordered pairs error. A zero-row, two-column
#'   matrix disables pair checks. More than 100 automatic pairs warns.
#' @param codebook Optional dataframe with exact raw names in Variable and optional
#'   Label and Category columns. Nonblank explicit labels override attached labels;
#'   categories default to Uncategorized. Explicit Type = clock/duration plus
#'   TimeUnit/TimeFormat control temporal modeling copies; other Type values do
#'   not trigger recoding or conversion. Original inputs remain unchanged. Identical duplicate
#'   entries collapse; conflicting duplicates error. Missing/unmatched metadata and
#'   unused entries are diagnosed. Resolved metadata is retained in
#'   Metadata$VariableMetadata, outside detection and explanation fingerprints.
#' @param time_types Optional named clock/duration overrides, resolved before
#'   codebook declarations then attached scidr_type attributes. Plain difftime
#'   defaults to duration; unclassified hms requires a declaration. Names must
#'   identify selected model columns. No temporal meaning is inferred from names.
#' @param time_units Named input units for untyped numeric temporal columns:
#'   secs, mins, hours, days or weeks (usual aliases accepted). Typed difftime/hms
#'   objects use their own physical unit attributes.
#' @param time_formats Named HMS, HM or MS parsing formats for text. Clock values
#'   support HMS/HM; durations support all three. See [RevalueData()].
#' @details Durations use local seconds. Clock fields use paired sine/cosine
#'   features over 24 hours; model-to-source maps and exclusions are retained in
#'   Metadata$FeatureMaps and declarations in Metadata$TimeConversions. Invalid
#'   clocks are diagnosed and their field values unavailable, never wrapped.
#'   Text parsing failures mask only affected modeling values; missing parsing
#'   configuration makes the column unavailable. Native missing-data policies
#'   remain in effect. Originals and classes remain in ReviewData.
#'   Sensitivity and optional SHAP explain original columns through the exact
#'   encoded prediction wrapper, keeping each clock pair together. Clock cellwise
#'   evidence is collapsed by source field; component predictions are technical
#'   details, not inferred replacement times. Conditional clock boundaries are
#'   displayed as clock regions, not linear clock limits or clinical ranges.
#'   LOF uses circular distances for clock fields in mixed data and jointly scaled
#'   paired coordinates in numeric data. Dates/time stamps remain separate: clock
#'   fields do not enter date_pairs or imply an elapsed interval. Circular
#'   longitudinal trajectories and stripping dates from timestamps are not supported.
#'   An hms column left undeclared is modelled as plain seconds and recorded in an
#'   `ambiguous_temporal_seconds` diagnostic: correct for a duration, and for a
#'   clock it loses only wrap-around handling, so 23:50 and 00:10 read as far
#'   apart. Declare `Type` to model a clock circularly.
#' @param explain_isolation Automatically screen columns for selected isolation
#'   records? Default TRUE; FALSE skips this cost. Requires optional digest.
#' @param explanation_args Named list. Shared settings: method = "sensitivity"
#'   (or "shap"), top_n = 10 per stratum, rows = NULL (unique original row
#'   references override top_n), background_n = 50. Fast sensitivity settings:
#'   replacements = 3 reference rows, batch_size = 64 features (automatically
#'   reduced to bound wide input batches). SHAP-only settings: nsim = 64,
#'   repeats = 3; requires optional fastshap >= 0.1.1. All counts are positive
#'   integers; nsim and background_n must be >= 2. Every actual model input is
#'   screened. No model is refit or saved; computation is serial and reproducible.
#' @section Conditional pattern:
#' `conditional` uses `outliertree` to find values that are unusual within a
#' subgroup learned from the other columns. A numeric column with fewer than
#' `method_args$conditional$min_distinct` distinct values (default 10) is
#' compared as a category rather than as a number. Such columns -- counts that
#' are usually zero, rarely-endorsed items -- have a near-zero spread, so every
#' departure from the dominant value otherwise reads as an extreme conditional
#' outlier and the method flags much of the sample. As a category, a rare value
#' is judged by how often it occurs instead. Date and clock columns are never
#' converted, since they are few-valued by nature in small studies and the
#' explanations depend on their formatting. A
#' `conditional_discrete_as_categorical` diagnostic names the converted columns;
#' set `min_distinct = 0` to disable the conversion entirely.
#'
#' @section Unexpected missingness:
#' The default `"missingness"` method asks whether a blank is unexpected for
#' that record, not whether the column is mostly complete. For each selected
#' field with both observed and missing values it fits one small `rpart`
#' classification tree per analysis group, predicting whether the field is
#' blank from the other selected fields' values and from their own
#' missing/observed indicators. A blank is flagged only when its leaf holds at
#' least `min_reference` records of which at least `min_observed_fraction`
#' carry a value, so a field that is deliberately skipped for part of the
#' sample is not reported for the records that skip it. The leaf's decision
#' path is retained as the finding's comparison conditions and is what the
#' explanation names. Settings in `method_args$missingness` are
#' `min_reference = 20`, `min_observed_fraction = 0.95`, `cp = 0.01`,
#' `maxdepth = 4`, `max_levels = 20`, and `max_predictors = 200`; categorical
#' fields with more than
#' `max_levels` values contribute their indicator but no value splits. On wide
#' frames, where every other field contributes both a value column and an
#' indicator, fitting against all of them dominates the run time, so when the
#' candidates exceed `max_predictors` only the strongest by marginal
#' association with the blank are offered to each tree and a
#' `missingness_predictors_screened` diagnostic is recorded. The screen ranks
#' candidates only -- the tree still searches the survivors in full -- but a
#' predictor that matters solely through an interaction, with no marginal
#' signal of its own, can be screened out; raise `max_predictors` to widen the
#' search at a proportional cost in time. When no value rule
#' is found, the blank is judged against the group of variables it goes missing
#' with. Those groups are discovered, not declared: columns are clustered on how
#' similarly they go blank (complete-linkage on Jaccard distance over the
#' missing indicator), which recovers the forms a wide merged frame is built
#' from without knowing what any column measures. A blank is reported when the
#' record carries most of its group -- the gap is at most `max_block_gap` of the
#' group -- while at least `min_observed_fraction` of the records carrying any of
#' that group carry all of it. A record missing most of a group has skipped the
#' form and is not reported. Group settings are `min_block_events = 5` (blanks a
#' column needs before it can join a group), `min_block_size = 3`,
#' `block_similarity = 0.8` and `max_block_gap = 0.5`; the discovered groups are
#' returned in `Metadata$Settings$missingness$Blocks` and shown on the explorer
#' page so the grouping can be checked before its findings are believed. A blank
#' explained by neither a value rule nor a group has no comparison group at all
#' and is **not** reported: a `missingness_unexplained` diagnostic is recorded,
#' and the column appears in the explorer's completeness table instead. This is
#' deliberate -- "this column is mostly full and you are blank" is a fact about
#' the column, not an anomaly about the record. The reviewed row is excluded; with an ID,
#' all records with that ID are excluded. Missing IDs use row exclusion.
#' The leaf's observed fraction is a descriptive score,
#' not a validated prediction or error probability. `thresholds["missingness"]`
#' overrides the observed-fraction cutoff, retaining minimum reference support;
#' `flag_percentile` does not apply. Only `is.na()` values are missing: sentinel
#' codes are not automatically recoded. Skip logic and selective testing can
#' explain flags. No values are filled. All-missing columns cannot be assessed.
#' Comparisons are cached across fields, with cost increasing with the number
#' of selected variables and comparison groups. Exclude `"missingness"` from
#' `methods` to disable this screening.
#' @param Relabel Use attached labels; record missing labels and fall back to names.
#' @param flag_percentile Default cutoff for score-only methods, as a fraction
#'   strictly between 0 and 1; default 0.99 flags the top 1 percent of scores
#'   within each method and stratum. Applied only where a method produced scores
#'   but no flags at all, so native conditional, date and cellwise flags and any
#'   explicit `thresholds` entry always take precedence. Small strata can
#'   legitimately yield no flags. NULL leaves score-only methods unflagged, in
#'   which case ranks are the only ordering. A percentile cutoff is a review
#'   budget chosen for convenience, not an estimated error rate.
#' @param seed Non-negative integer seed; the caller's RNG state is restored.
#' @details FeatureSensitivity retains original values, signed ScoreDecrease
#'   (original score minus mean score after replacing one column with reference
#'   values), replacement-score mean/range/count and positive ranks. Positive
#'   effects suggest columns to inspect, not errors; effects are not additive.
#'   The same up-to-three reference rows are used across columns and selected
#'   records within each stratum. Native missing-value handling is retained.
#'   Replacement ranges are not confidence intervals or error probabilities.
#'   FeatureContributions is exclusively for explicit method = "shap" runs:
#'   adjusted additive contributions with simulation SD/range and positive ranks.
#'   ExplanationSummary records Method, score, status (Explained, Not requested
#'   or Unavailable), and SHAP baseline/reconstruction fields (NA for sensitivity).
#'   ExplanationMetadata records settings, model inputs, transformations, sampled
#'   background/reference rows, seeds, elapsed screening time and a separate
#'   fingerprint. Elapsed time is excluded from that fingerprint.
#'   Background observations may contain errors. Replacing correlated or derived
#'   values can create unrealistic combinations; neither method establishes
#'   independent effects, clinical implausibility, interactions or corrections.
#'   SHAP cost grows with records, features, nsim and repeats and can be substantial
#'   for wide datasets. Explanation failures preserve successful detection scores.
#' @return A `DataAnomalyResult` list with `Results` (one record per method),
#'   `CellResults`, `IntervalResults`, `ParticipantSummary`, `Summary`,
#'   `Diagnostics`, `ReviewData`, `DateResults`, `Metadata`, `FeatureSensitivity`,
#'   `FeatureContributions`, `ExplanationSummary`, and `ExplanationMetadata`. Larger scores are more unusual;
#'   ranks are within method and stratum. Unscored rows are not treated as normal.
#' @details Dates and POSIXct timestamps are modeled locally as elapsed days from
#'   per-column reference dates recorded in `Metadata$DateConversions`. Original
#'   classes/time zones remain in `ReviewData`. Date-like strings are not inferred.
#'   The dates method screens columns and signed within-row intervals per stratum
#'   with outer Tukey fences (type-7 quartiles). Negative intervals are not
#'   inherently incorrect. With zero IQR, departures from a modal value are flagged
#'   only if its observed fraction meets majority_fraction; otherwise unavailable.
#'   Scores measure excess beyond limits divided by max(IQR, 1 day), aggregated
#'   by maximum per record. These are review priorities, not error probabilities.
#'   DateResults contains every check, original readable dates, sample sizes,
#'   observed ranges, fence values, rules, status, flags and optional ID/context.
#'   Column ValueDays/ranges/fences use the recorded reference; interval values use
#'   signed days. LowerDate/UpperDate show readable column fences.
#'   Zero dates yields not applicable. One date yields only column checks.
#'   Chronology is learned without names, fixed calendar rules, or date parsing.
#'
#'   Conditional and isolation methods use native missing-data handling.
#'   Local and distance methods use complete rows. DDC uses native filtering.
#'   Numeric-only methods exclude categorical columns with diagnostics. No
#'   imputation or automatic correction is performed. Mixed-data LOF uses a
#'   quadratic-size Gower distance matrix; select manageable strata.
#'   Backend failures are recorded rather than replaced with another method.
#'   `.Row` is reserved. ID, time, and strata columns cannot use result field
#'   names (`Method`, `Stratum`, `Status`, `Score`, `Flag`, `Threshold`,
#'   `Explanation`, `Rank`); rename conflicting context columns first.
#' @examples
#' set.seed(24)
#' df_Example <- data.frame(ID = 1:80, x = rnorm(80))
#' df_Example$y <- df_Example$x + rnorm(80, sd = 0.1)
#' df_Example$y[1] <- 10
#' if (requireNamespace("robustbase", quietly = TRUE)) {
#'   anomalies <- DetectDataAnomalies(df_Example, c("x", "y"),
#'     id_var = "ID", methods = "distance")
#'   anomalies$Results
#' }
#' @family data anomalies
#' @seealso [CreateDataAnomalyReport()], [ExploreDataAnomalies()],
#'   [CreateAnomalyReviewQueue()], [DetectLongitudinalAnomalies()];
#'   `vignette("data-anomalies")` for the workflow.
#' @export
DetectDataAnomalies <- function(data, variables = NULL, id_var = NULL,
    stratify_by = NULL, methods = c("conditional", "isolation", "dates", "missingness"), thresholds = NULL,
    method_args = list(), Relabel = TRUE, seed = 2024, date_pairs = NULL,
    explain_isolation = TRUE, explanation_args = list(), codebook = NULL,
    time_types = NULL, time_units = NULL, time_formats = NULL,
    flag_percentile = 0.99) {
  codebook <- ScidrAnomalyCodebook(codebook)
  flag_percentile <- ScidrAnomalyPercentile(flag_percentile)
  explanation_settings <- ScidrAnomalyExplanationArgs(explanation_args, explain_isolation, methods, nrow(data))
  automatic_variables <- is.null(variables)
  excluded_variables <- character()
  if (automatic_variables) {
    if (!is.data.frame(data)) stop("data must be a dataframe.")
    eligible <- vapply(data, function(x) is.null(dim(x)) &&
      (ScidrAnomalyIsDate(x) || inherits(x, "difftime") || (is.numeric(x) && !inherits(x, "difftime")) || is.factor(x) || is.character(x) || is.logical(x)), logical(1))
    excluded_variables <- names(data)[!eligible]
    variables <- setdiff(names(data)[eligible], c(id_var, stratify_by))
    if (!length(variables)) stop("No eligible model variables remain.")
  }
  missingness_args <- if ("missingness" %in% methods) ScidrAnomalyMissingnessArgs(method_args$missingness) else NULL
  temporal <- ScidrTimeDeclarations(data[intersect(setdiff(variables, id_var), names(data))], codebook, time_types, time_units, time_formats, detection = TRUE)
  setup <- ScidrAnomalySetup(data, variables, id_var, stratify_by, methods,
    c("conditional", "isolation", "local", "distance", "cellwise", "dates", "missingness"), thresholds, method_args, Relabel, seed, allow_dates = TRUE)
  if ("missingness" %in% names(thresholds) && (thresholds[["missingness"]] < 0 || thresholds[["missingness"]] > 1)) stop("The missingness cutoff must be between 0 and 1.")
  if (length(excluded_variables)) setup$Add(issue = "unsupported_auto_excluded",
    detail = paste(excluded_variables, collapse = ", "))
  variables <- setup$variables
  vars_Dates <- variables[vapply(data[variables], ScidrAnomalyIsDate, logical(1))]
  if (!is.null(date_pairs) && !"dates" %in% methods) stop("date_pairs requires the dates method.")
  pairs <- if ("dates" %in% methods) ScidrAnomalyDatePairs(date_pairs, vars_Dates) else NULL
  df_ModelData <- data
  setup$date_conversions <- list()
  for (variable in vars_Dates) {
    days <- ScidrAnomalyDateDays(data[[variable]])
    reference <- if (any(is.finite(days))) floor(min(days[is.finite(days)])) else 0
    conversion <- list(ReferenceDays = reference, ReferenceDate = as.character(as.Date(reference, origin = "1970-01-01")),
      InputClass = class(data[[variable]]), TimeZone = attr(data[[variable]], "tzone", exact = TRUE), Units = "elapsed days")
    setup$date_conversions[[variable]] <- conversion
    df_ModelData[[variable]] <- days - reference
  }
  temporal$ModelUnit <- rep("secs", nrow(temporal))
  temporal$PeriodSeconds <- ifelse(temporal$TimeType == "clock", 86400, NA_real_)
  temporal$Status <- rep("Prepared", nrow(temporal))
  temporal$Reason <- rep(NA_character_, nrow(temporal))
  temporal$UnavailableRows <- rep(list(integer()), nrow(temporal))
  for (i in seq_len(nrow(temporal))) {
    variable <- temporal$Variable[i]
    seconds <- tryCatch(ScidrTimeSeconds(data[[variable]], temporal[i, ], strict = FALSE), error = function(e) {
      temporal$Status[i] <<- "Unavailable"
      temporal$Reason[i] <<- conditionMessage(e)
      setup$Add(issue = "temporal_parse_failure", detail = paste(variable, conditionMessage(e)), rows = which(!is.na(data[[variable]])))
      rep(NA_real_, nrow(data))
    })
    parse_failures <- attr(seconds, "parse_failures", exact = TRUE)
    if (length(parse_failures)) {
      setup$Add(issue = "temporal_parse_failure", detail = paste(variable, "has unparseable values; other rows retained."), rows = parse_failures)
      temporal$Status[i] <- "Partially available"
      temporal$Reason[i] <- "Some text values could not be parsed"
    }
    invalid <- which(!is.na(seconds) & (!is.finite(seconds) | (temporal$TimeType[i] == "clock" & (seconds < 0 | seconds >= 86400))))
    if (length(invalid)) {
      setup$Add(issue = "invalid_temporal_value", detail = paste(variable, "has values outside its declared temporal domain; not wrapped or corrected."), rows = invalid)
      seconds[invalid] <- NA_real_
      temporal$Status[i] <- "Partially available"
      temporal$Reason[i] <- "Values outside the declared temporal domain"
    }
    if (!any(is.finite(seconds))) {temporal$Status[i] <- "Unavailable"; if (is.na(temporal$Reason[i])) temporal$Reason[i] <- "No observed usable temporal values"}
    temporal$UnavailableRows[[i]] <- sort(unique(c(parse_failures, invalid, if (temporal$Status[i] == "Unavailable") which(!is.na(data[[variable]])) else integer())))
    df_ModelData[[variable]] <- as.numeric(seconds)
  }
  # An hms column left undeclared is ambiguous, not unusable: ScidrTimeDeclarations()
  # has already warned that it will be compared as plain seconds, and this is where
  # that happens. Without the conversion the column reaches outliertree or isotree
  # as hms, which reject the whole frame on an unsupported column type and take the
  # entire method down -- so an undeclared time column silently cost every finding
  # from those methods rather than only its own wrap-around handling.
  ambiguous <- intersect(attr(temporal, "scidr_ambiguous", exact = TRUE), names(df_ModelData))
  for (variable in ambiguous) df_ModelData[[variable]] <- as.numeric(data[[variable]])
  if (length(ambiguous)) setup$Add(issue = "ambiguous_temporal_seconds",
    detail = paste0(paste(ambiguous, collapse = ", "),
      ": compared as plain seconds. An hms column can hold a clock time or a duration; declare Type = clock/duration ",
      "in the codebook or supply time_types so a clock is modelled with wrap-around handling."))
  setup$time_metadata <- temporal
  feature_maps <- list()
  restore <- ScidrAnomalySeed(seed)
  on.exit(restore(), add = TRUE)
  packages <- c(conditional = "outliertree", isolation = "isotree", local = "dbscan", distance = "robustbase", cellwise = "cellWise")
  explanations <- list()
  results <- cells <- settings <- list()
  date_results <- NULL
  for (method in methods) {
    if (method == "missingness") {
      fit <- ScidrAnomalyMissingness(setup, data, df_ModelData, missingness_args)
      results[[method]] <- fit$Results
      cells[[length(cells) + 1L]] <- fit$Cells
      settings[[method]] <- list(Arguments = missingness_args, Blocks = fit$Blocks,
        Definition = "Descriptive observed fraction in the comparison group learned from the other columns, excluding the evaluated participant; not a prediction or error probability. The comparison group is either a rule learned from another column's values or a group of variables that go missing together; a blank explained by neither is not reported.")
      next
    }
    if (method == "dates") {
      args <- ScidrAnomalyArgs(method_args[[method]], c("min_observations", "iqr_multiplier", "majority_fraction"),
        list(min_observations = 20L, iqr_multiplier = 3, majority_fraction = 0.95))
      for (argument in names(args)) if (!is.numeric(args[[argument]]) || length(args[[argument]]) != 1L ||
          !is.finite(args[[argument]])) stop("Date settings must be finite numeric scalars.")
      if (args$min_observations < 2 || args$min_observations != floor(args$min_observations) ||
          args$iqr_multiplier <= 0 || args$majority_fraction <= 0.5 || args$majority_fraction > 1) stop("Invalid date settings: require integer min_observations >= 2, positive iqr_multiplier, and 0.5 < majority_fraction <= 1.")
      date_fit <- ScidrAnomalyDateScreen(setup, data, df_ModelData, vars_Dates, pairs, args)
      results[[method]] <- date_fit$Results
      date_results <- date_fit$DateResults
      settings[[method]] <- list(Arguments = args, Pairs = pairs,
        Direction = "Second minus First, ordered by selected columns", Transformations = setup$date_conversions)
      next
    }
    package <- packages[[method]]
    if (!requireNamespace(package, quietly = TRUE)) stop("Install package '", package, "' to use method '", method, "'.")
    backend <- switch(method, conditional = outliertree::outlier.tree, isolation = isotree::isolation.forest,
      local = dbscan::lof, distance = robustbase::covMcd, cellwise = cellWise::DDC)
    allowed <- setdiff(names(formals(backend)), c("data", "df", "x", "X", "cols_ignore", "save_outliers", "...", "seed", "control", "sample_weights", "column_weights", "categ_cols", "output_score", "output_dist", "output_imputations", "build_imputer", "recode_categ", "standardize_data", "DDCpars"))
    if (method == "conditional") allowed <- c(allowed, "min_distinct")
    if (method == "cellwise") allowed <- c("fracNA", "numDiscrete", "precScale", "tolProb", "corrlim", "combinRule", "includeSelf", "fastDDC", "nbngbrs", "transFun", "standType", "cleanNAfirst", "silent", "nLocScale", "returnBigXimp", "corrType")
    defaults <- switch(method, isolation = list(nthreads = 1L), conditional = list(nthreads = 1L, outliers_print = 0L, min_distinct = 10L), local = list(minPts = 5L), cellwise = list(silent = TRUE), list())
    args <- ScidrAnomalyArgs(method_args[[method]], allowed, defaults)
    settings[[method]] <- list(Arguments = args, BackendDefaults = formals(backend), Transformations = list(), DateConversions = setup$date_conversions)
    df_Result <- ScidrAnomalyRecords(nrow(data), method, setup$groups, thresholds)
    for (stratum in sort(unique(stats::na.omit(setup$groups)))) {
      rows <- which(setup$groups == stratum)
      df_X <- df_ModelData[rows, variables, drop = FALSE]
      numeric_only <- method %in% c("distance", "cellwise")
      if (numeric_only) {
        excluded <- names(df_X)[!vapply(df_X, is.numeric, logical(1))]
        if (length(excluded)) setup$Add(method, stratum, "excluded_variables", paste(excluded, collapse = ", "))
        df_X <- df_X[vapply(df_X, is.numeric, logical(1))]
      }
      usable <- vapply(df_X, function(x) length(unique(x[!is.na(x)])) > 1L, logical(1))
      if (any(!usable)) setup$Add(method, stratum, "constant_or_missing", paste(names(df_X)[!usable], collapse = ", "))
      df_X <- df_X[usable]
      if (!ncol(df_X)) {
        df_Result$Status[rows] <- if (numeric_only && !any(vapply(df_ModelData[variables], is.numeric, logical(1)))) "not applicable" else "unavailable"
        setup$Add(method, stratum, "insufficient_dimensions", "No usable variables.", rows)
        next
      }
      finite_rows <- Reduce(`&`, lapply(df_X, function(x) if (is.numeric(x)) is.na(x) | is.finite(x) else rep(TRUE, length(x))))
      keep <- if (method %in% c("local", "distance")) stats::complete.cases(df_X) else rowSums(!is.na(df_X)) > 0
      keep <- keep & finite_rows
      omitted <- rows[!keep]
      if (length(omitted)) setup$Add(method, stratum, "excluded_rows", "Missing or nonfinite inputs.", omitted)
      rows <- rows[keep]
      df_X <- df_X[keep, , drop = FALSE]
      df_X[] <- lapply(df_X, function(x) if (is.character(x) || is.logical(x)) factor(x) else x)
      df_SourceX <- df_X
      feature_map <- ScidrTimeFeatureMap(names(df_SourceX), temporal)
      full_feature_map <- feature_map
      excluded_components <- character()
      df_X <- ScidrTimeEncode(df_SourceX, feature_map)
      if (any(temporal$TimeType == "clock")) {
        clock_features <- unlist(feature_map[intersect(names(feature_map), temporal$Variable[temporal$TimeType == "clock"])], use.names = FALSE)
        flat <- vapply(df_X[clock_features], function(x) {
          x <- x[is.finite(x)]; length(x) < 2L || diff(range(x)) <= 1e-12
        }, logical(1))
        if (any(flat)) {
          excluded_components <- names(flat)[flat]
          setup$Add(method, stratum, "constant_clock_components", paste(names(flat)[flat], collapse = ", "))
          df_X <- df_X[setdiff(names(df_X), names(flat)[flat])]
          feature_map <- lapply(feature_map, intersect, y = names(df_X))
          feature_map <- feature_map[lengths(feature_map) > 0L]
          df_SourceX <- df_SourceX[names(feature_map)]
        }
      }
      feature_maps[[paste(method, stratum, sep = ":")]] <- list(FeatureMap = feature_map, ModelFeatures = names(df_X), SourceVariables = names(df_SourceX), FullFeatureMap = full_feature_map, ExcludedComponents = excluded_components)
      setup$time_adapter <- if (nrow(temporal)) feature_maps[[paste(method, stratum, sep = ":")]] else NULL
      tryCatch(withCallingHandlers({
        if (nrow(df_X) < 3L || (method %in% c("distance", "cellwise") && ncol(df_X) < 2L)) stop("Insufficient rows or dimensions.")
        if (method == "local" && nrow(df_X) <= args$minPts) stop("LOF requires more rows than minPts.")
        score <- rep(NA_real_, length(rows))
        flag <- rep(NA, length(rows))
        if (method == "isolation") {
          fit <- do.call(backend, c(list(data = df_X, seed = seed), args))
          score <- as.numeric(predict(fit, df_X, type = "score"))
          if (explain_isolation) explanations[[as.character(stratum)]] <- ScidrAnomalyExplainIsolation(
            fit, df_SourceX, rows, score, stratum, setup, explanation_settings)
        } else if (method == "local") {
          if (all(vapply(df_X, is.numeric, logical(1)))) {
            sds <- vapply(df_X, stats::sd, numeric(1))
            if (any(!is.finite(sds) | sds == 0)) stop("Constant columns after complete-case filtering.")
            for (variable in intersect(names(feature_map), temporal$Variable[temporal$TimeType == "clock"])) {
              group <- feature_map[[variable]]
              sds[group] <- sqrt(sum(sds[group]^2))
            }
            settings[[method]]$Transformations[[as.character(stratum)]] <- list(center = colMeans(df_X), scale = sds)
            input <- scale(df_X, scale = sds)
          } else {
            if (!requireNamespace("cluster", quietly = TRUE)) stop("Install package 'cluster' for mixed-data LOF.")
            input <- if (any(names(df_SourceX) %in% temporal$Variable[temporal$TimeType == "clock"])) ScidrTimeGower(df_SourceX, temporal) else cluster::daisy(df_X, metric = "gower")
            settings[[method]]$Transformations[[as.character(stratum)]] <- if (any(names(df_SourceX) %in% temporal$Variable[temporal$TimeType == "clock"])) "Gower distance with circular clock distance min(abs(delta), 86400-abs(delta))/43200; one weight per source field" else "Gower distance"
          }
          score <- as.numeric(do.call(backend, c(list(x = input), args)))
        } else if (method == "distance") {
          fit <- do.call(backend, c(list(x = df_X), args))
          score <- sqrt(stats::mahalanobis(df_X, fit$center, fit$cov))
        } else if (method == "conditional") {
          # A mostly-constant numeric column (a count that is usually zero, a
          # rarely-endorsed item) has a near-zero spread, so every departure from
          # the dominant value reads as an extreme conditional outlier and the
          # method flags most of the sample. Below min_distinct values such a
          # column is compared as a category, where a rare value is judged by how
          # often it occurs rather than by its distance in standard deviations.
          backend_args <- args
          min_distinct <- backend_args$min_distinct
          backend_args$min_distinct <- NULL
          df_Fit <- df_X
          # Dates and clock features stay numeric: they are few-valued by nature in
          # small studies, and discretising them would drop the date formatting the
          # explanation depends on.
          temporal_features <- unique(c(vars_Dates, temporal$Variable,
            unlist(feature_map[intersect(names(feature_map), c(vars_Dates, temporal$Variable))], use.names = FALSE)))
          discrete <- vapply(names(df_Fit), function(name) {
            x <- df_Fit[[name]]
            is.numeric(x) && !name %in% temporal_features &&
              length(unique(stats::na.omit(x))) < min_distinct
          }, logical(1))
          if (any(discrete)) {
            for (name in names(df_Fit)[discrete]) df_Fit[[name]] <- factor(as.character(df_Fit[[name]]))
            setup$Add(method, stratum, "conditional_discrete_as_categorical",
              paste(paste(names(df_Fit)[discrete], collapse = ", "),
                "compared as categories; fewer than", min_distinct, "distinct numeric values."))
          }
          fit <- do.call(backend, c(list(df = df_Fit), backend_args))
          findings <- predict(fit, df_Fit, nthreads = 1L, outliers_print = 0L)
          flag[] <- FALSE
          unformatted <- character()
          for (i in seq_along(rows)) if (length(findings[[i]]$suspicous_value$column)) {
            finding <- findings[[i]]
            variable <- finding$suspicous_value$column
            if (is.numeric(variable)) variable <- names(df_X)[variable]
            if (!length(variable) || !variable %in% names(df_X)) stop("Unrecognized OutlierTree target column.")
            flag[i] <- TRUE
            # The backend can fail to format its own report for a finding: an
            # invalid decimal count, or a field width beyond what sprintf allows.
            # That is a printing fault, not a detection fault, and it must not
            # discard this finding or the ones after it, so the report is rebuilt
            # from the same stored evidence the detail panel uses.
            explanation <- tryCatch(paste(utils::capture.output(print(findings[i])), collapse = "\n"),
              error = function(e) {
                unformatted <<- c(unformatted, paste0(variable, ": ", conditionMessage(e)))
                NA_character_
              })
            mapped <- ScidrTimeFinding(finding, variable, feature_map, temporal, data, rows[i])
            source_variable <- mapped$Variable
            if (!is.na(explanation)) for (display_var in names(feature_map)) explanation <- gsub(paste0("[", display_var, "]"),
              paste0("[", setup$labels[[display_var]], "]"), explanation, fixed = TRUE)
            date_involved <- source_variable %in% vars_Dates || any(vapply(mapped$Conditions,
              function(condition) condition$column %in% vars_Dates, logical(1)))
            if (date_involved && identical(source_variable, variable) && !any(vapply(mapped$Conditions, function(x) !is.null(x$native_component), logical(1))))
              explanation <- ScidrAnomalyDateExplanation(finding, variable, rows[i], setup, data)
            if (is.na(explanation)) explanation <- paste0(
              ScidrAnomalyDateExplanation(finding, variable, rows[i], setup, data),
              "\nThe detection backend could not format its own report for this finding; this description is built from the same stored evidence.")
            cells[[length(cells) + 1L]] <- tibble::tibble(.Row = rows[i], Method = method, Stratum = stratum,
              Variable = source_variable, Label = unname(setup$labels[source_variable]),
              Observed = if (source_variable %in% vars_Dates) ScidrAnomalyFormatDate(df_ModelData[[source_variable]][rows[i]], setup$date_conversions[[source_variable]]) else ScidrTimeObserved(data, source_variable, rows[i], temporal),
              Flag = TRUE, Explanation = explanation, Conditions = list(mapped$Conditions), GroupStatistics = list(mapped$Statistics),
              NativeComponents = list(if (nrow(temporal)) finding else NULL))
            df_Result$Explanation[rows[i]] <- explanation
          }
          if (length(unformatted)) setup$Add(method, stratum, "conditional_explanation_rebuilt",
            paste0(length(unformatted), " finding(s) kept with a rebuilt description after a backend formatting error: ",
              paste(utils::head(unformatted, 5), collapse = "; "),
              if (length(unformatted) > 5L) paste0("; and ", length(unformatted) - 5L, " more.") else "."))
        } else {
          fit <- do.call(backend, list(X = as.matrix(df_X), DDCpars = args))
          residuals <- fit$stdResid
          retained_rows <- unname(fit$rowInAnalysis)
          retained_columns <- unname(fit$colInAnalysis)
          if (is.null(residuals) || !identical(dim(residuals), c(length(retained_rows), length(retained_columns)))) stop("Unrecognized DDC row/column mapping.")
          removed <- setdiff(seq_along(rows), retained_rows)
          if (length(removed)) setup$Add(method, stratum, "ddc_excluded_rows", "DDC native filtering.", rows[removed])
          dropped <- setdiff(seq_len(ncol(df_X)), retained_columns)
          if (length(dropped)) setup$Add(method, stratum, "ddc_excluded_columns", paste(names(df_X)[dropped], collapse = ", "))
          rows <- rows[retained_rows]
          df_X <- df_X[retained_rows, retained_columns, drop = FALSE]
          nonfinite <- !is.na(residuals) & !is.finite(residuals)
          if (any(nonfinite)) {
            setup$Add(method, stratum, "nonfinite_residuals",
              paste("DDC produced nonfinite standardized residuals; affected cells are unavailable:",
                paste(names(df_X)[colSums(nonfinite) > 0], collapse = ", ")),
              rows[rowSums(nonfinite) > 0])
            residuals[nonfinite] <- NA_real_
          }
          native <- matrix(FALSE, nrow(df_X), ncol(df_X))
          native[fit$indcells] <- TRUE
          cell_pieces <- list()
          for (j in seq_len(ncol(df_X))) {
            values <- abs(residuals[, j])
            flags <- if (method %in% names(thresholds)) values >= thresholds[[method]] else native[, j]
            flags[is.na(values) | is.na(df_X[[j]])] <- NA
            model_variable <- names(df_X)[j]
            source_variable <- names(feature_map)[vapply(feature_map, function(x) model_variable %in% x, logical(1))]
            declaration <- temporal[match(source_variable, temporal$Variable), , drop = FALSE]
            kind <- if (nrow(declaration)) declaration$TimeType else NA_character_
            expected <- if (source_variable %in% vars_Dates) ScidrAnomalyFormatDate(fit$Xest[, j], setup$date_conversions[[source_variable]]) else
              if (!is.na(kind) && kind == "duration") ScidrTimeFormat(fit$Xest[, j], "duration") else as.character(fit$Xest[, j])
            cell_pieces[[j]] <- tibble::tibble(.Row = rows, Method = method, Stratum = stratum,
              Variable = source_variable, Label = unname(setup$labels[source_variable]),
              Observed = if (source_variable %in% vars_Dates) ScidrAnomalyFormatDate(df_ModelData[[source_variable]][rows], setup$date_conversions[[source_variable]]) else ScidrTimeObserved(data, source_variable, rows, temporal), ExpectedDisplay = expected,
              Expected = fit$Xest[, j], Residual = unname(residuals[, j]), Score = unname(values), Flag = unname(flags),
              TemporalType = kind, ModelFeature = model_variable, NativeComponents = rep(list(NULL), length(rows)))
          }
          df_CurrentCells <- ScidrTimeCollapseCells(dplyr::bind_rows(cell_pieces), temporal)
          if (!nrow(temporal)) df_CurrentCells <- df_CurrentCells[setdiff(names(df_CurrentCells), c("TemporalType", "ModelFeature", "NativeComponents"))]
          cells[[length(cells) + 1L]] <- df_CurrentCells
          score <- apply(abs(residuals), 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
          native[is.na(residuals) | is.na(as.matrix(df_X))] <- FALSE
          flag <- apply(native, 1, any)
          score[rowSums(!is.na(df_X)) == 0] <- NA_real_
        }
        if (method %in% names(thresholds)) flag <- score >= thresholds[[method]]
        evaluated <- if (method == "conditional") rep(TRUE, length(rows)) else !is.na(score)
        if (any(!evaluated)) setup$Add(method, stratum, "unscored_backend_rows", "Backend returned no usable score.", rows[!evaluated])
        df_Result$Status[rows[evaluated]] <- "evaluated"
        df_Result$Score[rows] <- score
        df_Result$Flag[rows[evaluated]] <- flag[evaluated]
      }, warning = function(w) {
        setup$Add(method, stratum, "backend_warning", conditionMessage(w), rows)
        invokeRestart("muffleWarning")
      }), error = function(e) setup$Add(method, stratum, "backend_failure", conditionMessage(e), rows))
    }
    results[[method]] <- df_Result
  }
  result <- ScidrAnomalyFinish(setup, results, cells, list(), match.call(), "CrossSectionalAnomalyResult", settings, date_results, flag_percentile)
  if (!nrow(temporal) && "NativeComponents" %in% names(result$CellResults)) result$CellResults$NativeComponents <- NULL
  result$Metadata$TimeConversions <- temporal
  result$Metadata$FeatureMaps <- feature_maps
  result$FeatureContributions <- dplyr::bind_rows(c(list(ScidrAnomalyEmptyContributions()), lapply(explanations, `[[`, "Contributions")))
  result$FeatureSensitivity <- dplyr::bind_rows(c(list(ScidrAnomalyEmptySensitivity()), lapply(explanations, `[[`, "Sensitivity")))
  isolation <- result$Results[result$Results$Method == "isolation", ]
  result$ExplanationSummary <- tibble::tibble(.Row = isolation$.Row, Stratum = isolation$Stratum,
    Score = isolation$Score, Baseline = NA_real_, SumContributions = NA_real_, ReconstructionResidual = NA_real_,
    ExplainedVariables = 0L, Status = ifelse(isolation$Status == "evaluated", "Not requested", "Unavailable"),
    Reason = ifelse(isolation$Status == "evaluated", "Outside explanation selection or explanations disabled", "Native score unavailable"))
  for (explanation in explanations) {
    index <- match(explanation$Summary$.Row, result$ExplanationSummary$.Row)
    result$ExplanationSummary[index, names(explanation$Summary)] <- explanation$Summary
  }
  result$ExplanationSummary$Method <- explanation_settings$method
  result$ExplanationMetadata <- list(Method = explanation_settings$method, Enabled = explain_isolation && "isolation" %in% methods,
    Settings = explanation_settings, Strata = lapply(explanations, `[[`, "Metadata"),
    PackageVersions = c(result$Metadata$PackageVersions[c("SciDataReportR", "isotree")], fastshap = if (requireNamespace("fastshap", quietly = TRUE)) as.character(utils::packageVersion("fastshap")) else NA_character_),
    Definitions = if (explanation_settings$method == "sensitivity") "Score decrease = original isolation score minus mean score after replacing one column with reference values. Positive suggests a column to inspect, not an error. Effects are not additive; replacement ranges are not confidence intervals. Replacements may create unrealistic combinations. Background may contain errors." else "Approximate adjusted SHAP on the exact native score. Positive raises score relative to background; negative lowers it. Background may contain unusual/error records. Simulation variability is not error probability. Correlated/derived columns share attribution; no independence or interaction claim.")
  if (explain_isolation && "isolation" %in% methods) result$ExplanationMetadata$Fingerprint <- ScidrAnomalyExplanationFingerprint(result)
  ScidrAnomalyAttachCodebook(result, codebook)
}
