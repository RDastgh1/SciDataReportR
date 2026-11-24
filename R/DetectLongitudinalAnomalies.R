#' Detect unusual longitudinal changes and trajectories
#'
#' Screen long-format repeated observations using empirical change rates,
#' leave-one-participant-out categorical transitions, or in-sample mixed-model
#' residuals. These diagnostics do not establish that either endpoint is wrong.
#' @inheritParams DetectDataAnomalies
#' @param id_var Required participant identifier column name.
#' @param thresholds Optional named cutoffs for change, trajectory, or transition
#'   scores. Values at least the cutoff are flagged; otherwise flags are `NA`.
#' @param time_var Required numeric, Date, or POSIXct time column. Dates and
#'   timestamps are converted to elapsed days; numeric times keep their units.
#' @param methods Subset of `change`, `transition`, and `trajectory`.
#' @param codebook Optional labels/categories for presentation. Temporal modeling
#'   declarations apply to DetectDataAnomalies(), not this longitudinal detector.
#'   Circular longitudinal outcomes and clock-only ordering are not supported.
#' @param covariates Optional character vector of fixed-effect covariates for
#'   trajectory models. Time is included automatically.
#' @param method_args Named lists keyed by method. Change and transition have no
#'   tunable settings. Trajectory accepts `control` for [nlme::lmeControl()].
#' @return A `DataAnomalyResult` with the same components as
#'   [DetectDataAnomalies()]. `IntervalResults` links both original endpoint
#'   rows. `CellResults` contains trajectory fitted values and residuals.
#' @details Change scores use median and MAD of change per elapsed time, requiring
#'   ten intervals from five participants. Zero MAD leaves scores unavailable.
#'   Transitions condition on the preceding category and exclude all transitions
#'   from the scored participant, requiring ten reference transitions from five
#'   other participants. Missing observations are bridged to the next observed
#'   value; elapsed time is retained. Transition probabilities are unsmoothed.
#'   Trajectories use linear elapsed time and a random intercept, requiring three
#'   observed complete visits per eligible participant and ten eligible
#'   participants. Scores are absolute in-sample Pearson residuals, not calibrated
#'   prediction errors. Sparse histories and failed fits remain unavailable.
#'   No score-based method flags observations unless a threshold is supplied.
#' @examples
#' df_Visits <- expand.grid(ID = seq_len(12), Time = 0:3)
#' df_Visits$Value <- sin(df_Visits$ID) + df_Visits$Time
#' result <- DetectLongitudinalAnomalies(df_Visits, "Value", "ID", "Time",
#'   methods = "change")
#' result$IntervalResults
#' @param flag_percentile Default cutoff for score-only methods, as in
#'   [DetectDataAnomalies()]; default 0.99. NULL leaves change, trajectory and
#'   transition scores unflagged.
#' @family data anomalies
#' @seealso [DetectDataAnomalies()], [CreateDataAnomalyReport()],
#'   [PlotDataAnomalyTrajectories()]
#' @export
DetectLongitudinalAnomalies <- function(data, variables, id_var, time_var,
    stratify_by = NULL, methods = c("change", "transition"), covariates = NULL,
    thresholds = NULL, method_args = list(), Relabel = TRUE, seed = 2024, codebook = NULL,
    flag_percentile = 0.99) {
  codebook <- ScidrAnomalyCodebook(codebook)
  flag_percentile <- ScidrAnomalyPercentile(flag_percentile)
  if (missing(id_var) || is.null(id_var) || missing(time_var) || is.null(time_var)) stop("id_var and time_var are required.")
  setup <- ScidrAnomalySetup(data, variables, id_var, stratify_by, methods,
    c("change", "transition", "trajectory"), thresholds, method_args, Relabel, seed, time_var, covariates)
  if (!is.numeric(data[[time_var]]) && !inherits(data[[time_var]], c("Date", "POSIXct"))) stop("time_var must be numeric, Date, or POSIXct.")
  time <- if (inherits(data[[time_var]], "POSIXct")) as.numeric(data[[time_var]]) / 86400 else as.numeric(data[[time_var]])
  valid <- !is.na(data[[id_var]]) & is.finite(time) & !is.na(setup$groups)
  key <- data[c(id_var, time_var, stratify_by)]
  if (any(duplicated(key[valid, , drop = FALSE]))) stop("Duplicate participant-time keys within strata.")
  if (any(!valid)) setup$Add(issue = "missing_longitudinal_key", detail = "Missing ID, time, or strata.", rows = which(!valid))
  restore <- ScidrAnomalySeed(seed)
  on.exit(restore(), add = TRUE)
  results <- cells <- intervals <- settings <- list()
  for (method in methods) {
    args <- ScidrAnomalyArgs(method_args[[method]], if (method == "trajectory") "control" else character())
    settings[[method]] <- list(Arguments = args, TimeUnits = if (inherits(data[[time_var]], c("Date", "POSIXct"))) "days" else "input units",
      MissingValues = "Consecutive observed values; no imputation")
    if (method == "trajectory" && !requireNamespace("nlme", quietly = TRUE)) stop("Install package 'nlme' to use trajectory detection.")
    df_Result <- ScidrAnomalyRecords(nrow(data), method, setup$groups, thresholds)
    numeric_vars <- vapply(data[variables], is.numeric, logical(1))
    vars_Eligible <- variables[if (method == "transition") !numeric_vars else numeric_vars]
    if (!length(vars_Eligible)) df_Result$Status[valid] <- "not applicable"
    for (variable in setdiff(variables, vars_Eligible)) setup$Add(method, issue = "excluded_variable", detail = variable)
    for (stratum in sort(unique(setup$groups[valid]))) for (variable in vars_Eligible) {
      rows <- which(valid & setup$groups == stratum & !is.na(data[[variable]]))
      if (is.numeric(data[[variable]])) rows <- rows[is.finite(data[[variable]][rows])]
      if (method == "trajectory" && length(covariates)) rows <- rows[stats::complete.cases(data[rows, covariates, drop = FALSE])]
      excluded <- setdiff(which(valid & setup$groups == stratum), rows)
      if (length(excluded)) setup$Add(method, stratum, "excluded_rows", paste("Missing/nonfinite values or incomplete covariates:", variable), excluded)
      rows <- rows[order(data[[id_var]][rows], time[rows])]
      tryCatch(withCallingHandlers({
        if (method == "trajectory") {
          counts <- table(as.character(data[[id_var]][rows]))
          eligible <- as.character(data[[id_var]][rows]) %in% names(counts)[counts >= 3L]
          if (any(!eligible)) setup$Add(method, stratum, "sparse_history", paste("Fewer than three complete visits:", variable), rows[!eligible])
          rows <- rows[eligible]
          if (length(unique(data[[id_var]][rows])) < 10L) stop("Trajectory requires ten participants with at least three complete visits each.")
          df_Model <- data.frame(.Value = data[[variable]][rows], .Time = time[rows] - min(time[rows]), .ID = factor(data[[id_var]][rows]))
          if (length(covariates)) {
            df_Cov <- data[rows, covariates, drop = FALSE]
            names(df_Cov) <- paste0(".Cov", seq_along(covariates))
            df_Model <- cbind(df_Model, df_Cov)
          }
          formula <- stats::reformulate(c(".Time", if (length(covariates)) paste0(".Cov", seq_along(covariates))), response = ".Value")
          fit <- do.call(nlme::lme, c(list(fixed = formula, random = ~ 1 | .ID, data = df_Model, na.action = stats::na.fail), args))
          score <- abs(as.numeric(stats::residuals(fit, type = "pearson")))
          cells[[length(cells) + 1L]] <- tibble::tibble(.Row = rows, Method = method, Stratum = stratum,
            Variable = variable, Label = unname(setup$labels[variable]), Observed = as.character(data[[variable]][rows]),
            Expected = as.numeric(stats::fitted(fit)), Residual = as.numeric(stats::residuals(fit)),
            Score = score, Flag = if (method %in% names(thresholds)) score >= thresholds[[method]] else NA,
            Explanation = "In-sample random-intercept trajectory residual")
          settings[[method]][[paste(stratum, variable, sep = ":")]] <- list(Formula = formula, TimeOrigin = min(time[rows]))
          next
        }
        pieces <- lapply(split(rows, as.character(data[[id_var]][rows])), function(r) {
          if (length(r) < 2L) return(NULL)
          tibble::tibble(PreviousRow = head(r, -1L), .Row = tail(r, -1L))
        })
        df_Pairs <- dplyr::bind_rows(pieces)
        if (!nrow(df_Pairs)) stop("No consecutive observed pairs.")
        previous <- df_Pairs$PreviousRow
        current <- df_Pairs$.Row
        df_Pairs$Participant <- as.character(data[[id_var]][current])
        df_Pairs$Method <- method
        df_Pairs$Stratum <- stratum
        df_Pairs$Variable <- variable
        df_Pairs$Label <- unname(setup$labels[variable])
        df_Pairs$PreviousValue <- as.character(data[[variable]][previous])
        df_Pairs$Observed <- as.character(data[[variable]][current])
        df_Pairs$Elapsed <- time[current] - time[previous]
        df_Pairs$Score <- NA_real_
        if (method == "change") {
          df_Pairs$Change <- data[[variable]][current] - data[[variable]][previous]
          df_Pairs$Rate <- df_Pairs$Change / df_Pairs$Elapsed
          if (nrow(df_Pairs) >= 10L && length(unique(df_Pairs$Participant)) >= 5L) {
            scale <- stats::mad(df_Pairs$Rate)
            if (is.finite(scale) && scale > 0) df_Pairs$Score <- abs(df_Pairs$Rate - stats::median(df_Pairs$Rate)) / scale
            else setup$Add(method, stratum, "zero_mad", paste("Unscored change rates:", variable), current)
          } else setup$Add(method, stratum, "insufficient_support", paste("Change needs ten intervals from five participants:", variable), current)
        } else {
          df_Pairs$TransitionCount <- df_Pairs$ConditioningCount <- df_Pairs$ReferenceParticipants <- NA_integer_
          for (i in seq_len(nrow(df_Pairs))) {
            reference <- which(df_Pairs$Participant != df_Pairs$Participant[i] & df_Pairs$PreviousValue == df_Pairs$PreviousValue[i])
            denominator <- length(reference)
            numerator <- sum(df_Pairs$Observed[reference] == df_Pairs$Observed[i])
            participants <- length(unique(df_Pairs$Participant[reference]))
            df_Pairs$TransitionCount[i] <- numerator
            df_Pairs$ConditioningCount[i] <- denominator
            df_Pairs$ReferenceParticipants[i] <- participants
            if (denominator >= 10L && participants >= 5L) df_Pairs$Score[i] <- 1 - numerator / denominator
          }
          if (anyNA(df_Pairs$Score)) setup$Add(method, stratum, "insufficient_support", paste("Transition reference support:", variable), current[is.na(df_Pairs$Score)])
        }
        df_Pairs$Flag <- if (method %in% names(thresholds)) df_Pairs$Score >= thresholds[[method]] else NA
        df_Pairs$Status <- ifelse(is.na(df_Pairs$Score), "unavailable", "evaluated")
        intervals[[length(intervals) + 1L]] <- df_Pairs
      }, warning = function(w) {
        setup$Add(method, stratum, "backend_warning", conditionMessage(w), rows)
        invokeRestart("muffleWarning")
      }), error = function(e) setup$Add(method, stratum, "backend_failure", paste(variable, conditionMessage(e)), rows))
    }
    # Aggregate maximum available evidence without attributing an interval to one endpoint.
    evidence <- if (method == "trajectory") dplyr::bind_rows(cells) else dplyr::bind_rows(intervals)
    if (nrow(evidence)) {
      evidence <- evidence[evidence$Method == method & !is.na(evidence$Score), , drop = FALSE]
      if (nrow(evidence) && method != "trajectory") {
        df_Previous <- evidence
        df_Previous$.Row <- df_Previous$PreviousRow
        evidence <- dplyr::bind_rows(evidence, df_Previous)
      }
      for (row in unique(evidence$.Row)) {
        scores <- evidence$Score[evidence$.Row == row]
        df_Result$Status[row] <- "evaluated"
        df_Result$Score[row] <- max(scores)
        df_Result$Flag[row] <- if (method %in% names(thresholds)) max(scores) >= thresholds[[method]] else NA
        df_Result$Explanation[row] <- if (method == "trajectory") "Maximum in-sample trajectory residual across selected variables" else "Maximum score among intervals involving this row; either endpoint may warrant review"
      }
    }
    results[[method]] <- df_Result
  }
  result <- ScidrAnomalyFinish(setup, results, cells, intervals, match.call(), "LongitudinalAnomalyResult", settings, NULL, flag_percentile)
  ScidrAnomalyAttachCodebook(result, codebook)
}
