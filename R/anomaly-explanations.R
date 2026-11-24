ScidrAnomalyExplanationArgs <- function(args, enabled, methods, n) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) stop("explain_isolation must be TRUE or FALSE.")
  defaults <- list(method = "sensitivity", top_n = 10L, rows = NULL, background_n = 50L,
    replacements = 3L, batch_size = 64L, nsim = 64L, repeats = 3L)
  args <- ScidrAnomalyArgs(args, names(defaults), defaults)
  if (!is.character(args$method) || length(args$method) != 1L || is.na(args$method) || !args$method %in% c("sensitivity", "shap")) stop("explanation_args$method must be sensitivity or shap.")
  for (name in c("top_n", "background_n", "replacements", "batch_size", "nsim", "repeats")) {
    x <- args[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < if (name %in% c("background_n", "nsim")) 2 else 1) stop("Invalid explanation_args$", name, ".")
  }
  if (!is.null(args$rows) && (!is.numeric(args$rows) || !length(args$rows) || anyNA(args$rows) || anyDuplicated(args$rows) || any(!is.finite(args$rows) | args$rows != floor(args$rows) | args$rows < 1 | args$rows > n))) stop("explanation_args$rows must be unique original row references.")
  if (!"isolation" %in% methods && !identical(args, defaults)) stop("explanation_args requires the isolation method.")
  if (enabled && "isolation" %in% methods && args$method == "shap" && (!requireNamespace("fastshap", quietly = TRUE) || utils::packageVersion("fastshap") < "0.1.1")) stop("Install fastshap >= 0.1.1 for isolation explanations, or set explain_isolation = FALSE.")
  args
}

ScidrAnomalyEmptyContributions <- function() tibble::tibble(.Row = integer(), Stratum = integer(),
  Variable = character(), Label = character(), Observed = character(), Contribution = double(),
  SimulationSD = double(), SimulationMin = double(), SimulationMax = double(), PositiveRank = double())

ScidrAnomalyExplainIsolation <- function(fit, df_X, rows, scores, stratum, setup, args) {
  if (args$method == "sensitivity") return(ScidrAnomalyIsolationSensitivity(fit, df_X, rows, scores, stratum, setup, args))
  # Keep background/newdata containers identical for older fastshap versions.
  df_X <- as.data.frame(df_X, optional = TRUE)
  restore <- ScidrAnomalySeed(setup$seed)
  on.exit(restore(), add = TRUE)
  Seed <- function(row, repetition) as.integer((as.double(setup$seed) + 1000003 * as.double(stratum) + 1009 * row + repetition) %% (.Machine$integer.max - 1))
  eligible <- which(is.finite(scores))
  selected <- if (is.null(args$rows)) head(eligible[order(-scores[eligible], rows[eligible])], args$top_n) else eligible[rows[eligible] %in% args$rows]
  df_Summary <- tibble::tibble(.Row = rows, Stratum = as.integer(stratum), Score = scores,
    Baseline = NA_real_, SumContributions = NA_real_, ReconstructionResidual = NA_real_,
    ExplainedVariables = 0L, Status = ifelse(is.finite(scores), "Not requested", "Unavailable"),
    Reason = ifelse(is.finite(scores), "Outside explanation selection", "Native score unavailable"))
  metadata <- list(SelectedRows = rows[selected], Settings = args, SourceVariables = names(df_X), ModelingVariables = if (is.null(setup$time_adapter)) names(df_X) else setup$time_adapter$ModelFeatures,
    InputClasses = lapply(df_X, class), FactorLevels = lapply(df_X[vapply(df_X, is.factor, logical(1))], levels),
    DateConversions = setup$date_conversions, TimeConversions = setup$time_metadata, ModelAdapter = setup$time_adapter, PredictionThreads = 1L, BackgroundRows = integer(), Seeds = list())
  contributions <- list(ScidrAnomalyEmptyContributions())
  if (!length(selected)) return(list(Contributions = contributions[[1]], Summary = df_Summary, Metadata = metadata))
  message("Explaining isolation: ", length(selected), " records x ", ncol(df_X), " features x ", args$nsim,
    " simulations x ", args$repeats, " repeats (stratum ", stratum, ").")
  df_Summary$Status[selected] <- "Unavailable"
  df_Summary$Reason[selected] <- "Explanation did not complete"
  background_failure <- tryCatch({
    set.seed(Seed(0, 0))
    background <- sample.int(nrow(df_X), min(args$background_n, nrow(df_X)))
    df_Background <- df_X[background, , drop = FALSE]
    metadata$BackgroundRows <- rows[background]
    metadata$BackgroundSeed <- Seed(0, 0)
    Predict <- function(object, newdata) as.numeric(stats::predict(object, ScidrTimePredictData(newdata, setup$time_adapter), type = "score", nthreads = 1L))
    baseline <- mean(Predict(fit, df_Background))
    if (!is.finite(baseline)) stop("Background baseline is nonfinite.")
    df_Summary$Baseline[selected] <- baseline
    for (i in selected) {
      message("Explaining source row ", rows[i], " (", match(i, selected), "/", length(selected), ").")
      seeds <- vapply(seq_len(args$repeats), function(repetition) Seed(rows[i], repetition), integer(1))
      metadata$Seeds[[as.character(rows[i])]] <- seeds
      row_failure <- tryCatch({
        values <- vapply(seq_len(args$repeats), function(repetition) {
          set.seed(seeds[repetition])
          result <- withCallingHandlers(fastshap::explain(fit, X = df_Background, newdata = df_X[i, , drop = FALSE],
            pred_wrapper = Predict, nsim = args$nsim, adjust = TRUE, baseline = baseline, parallel = FALSE), warning = function(w) {
              setup$Add("isolation", stratum, "explanation_warning", conditionMessage(w), rows[i])
              invokeRestart("muffleWarning")
            })
          as.numeric(result[1, names(df_X)])
        }, numeric(ncol(df_X)))
        values <- matrix(values, nrow = ncol(df_X), ncol = args$repeats)
        if (any(!is.finite(values))) stop("Nonfinite feature contributions.")
        means <- rowMeans(values)
        residual <- scores[i] - baseline - sum(means)
        if (abs(residual) > 1e-6 * max(1, abs(scores[i]))) stop("Contributions do not reconstruct the isolation score within tolerance.")
        contributions[[length(contributions) + 1L]] <- tibble::tibble(.Row = rows[i], Stratum = as.integer(stratum),
          Variable = names(df_X), Label = unname(setup$labels[names(df_X)]),
          Observed = unname(vapply(names(df_X), function(variable) ScidrTimeObserved(setup$data, variable, rows[i], setup$time_metadata), character(1))),
          Contribution = means, SimulationSD = if (args$repeats > 1) apply(values, 1, stats::sd) else rep(NA_real_, nrow(values)),
          SimulationMin = apply(values, 1, min), SimulationMax = apply(values, 1, max),
          PositiveRank = ifelse(means > 0, rank(-replace(means, means <= 0, NA_real_), ties.method = "min", na.last = "keep"), NA_real_))
        df_Summary$Baseline[i] <- baseline
        df_Summary$SumContributions[i] <- sum(means)
        df_Summary$ReconstructionResidual[i] <- residual
        df_Summary$ExplainedVariables[i] <- ncol(df_X)
        df_Summary$Status[i] <- "Explained"
        df_Summary$Reason[i] <- ""
        NULL
      }, error = function(e) {
        detail <- if (ncol(df_X) == 1L) paste0("fastshap could not explain this single-feature fit (version ", utils::packageVersion("fastshap"), "). ", conditionMessage(e), ". Native score retained; update fastshap or review the score directly.") else conditionMessage(e)
        detail
      })
      if (!is.null(row_failure)) {
        df_Summary$Reason[i] <- row_failure
        setup$Add("isolation", stratum, "explanation_failure", row_failure, rows[i])
      }
    }
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(background_failure)) {
    pending <- selected[df_Summary$Status[selected] != "Explained"]
    df_Summary$Reason[pending] <- background_failure
    setup$Add("isolation", stratum, "explanation_failure", background_failure, rows[pending])
  }
  message("Isolation explanations finished: ", sum(df_Summary$Status == "Explained"), "/", length(selected), " records.")
  list(Contributions = dplyr::bind_rows(contributions), Summary = df_Summary, Metadata = metadata)
}

ScidrAnomalyExplanationFingerprint <- function(object) {
  metadata <- object$ExplanationMetadata
  metadata$Fingerprint <- NULL
  if (identical(metadata$Method, "sensitivity")) metadata$Strata <- lapply(metadata$Strata, function(x) { x$ElapsedSeconds <- NULL; x })
  evidence <- list(object$FeatureContributions, object$ExplanationSummary, metadata)
  # Legacy explanation fingerprints used exactly these three elements.
  if (!is.null(object$FeatureSensitivity)) evidence$FeatureSensitivity <- object$FeatureSensitivity
  ScidrReviewHash(evidence)
}

ScidrReviewContributions <- function(object) {
  pieces <- lapply(names(object$Sources), function(analysis) {
    source <- object$Sources[[analysis]]
    df_Sensitivity <- source$FeatureSensitivity
    if (!is.null(df_Sensitivity)) names(df_Sensitivity)[names(df_Sensitivity) == "PositiveRank"] <- "PositiveSensitivityRank"
    df_Values <- dplyr::bind_rows(source$FeatureContributions, df_Sensitivity)
    if (is.null(df_Values) || !nrow(df_Values)) return(NULL)
    if (!is.null(source$ExplanationMetadata$Fingerprint) && !identical(source$ExplanationMetadata$Fingerprint,
        ScidrAnomalyExplanationFingerprint(source))) stop("Explanation fingerprint changed: ", analysis)
    if (!is.null(source$ExplanationMetadata$Method)) df_Values$ExplanationMethod <- source$ExplanationMetadata$Method
    catalog <- ScidrReviewCatalog(object, analysis)
    df_Values <- ScidrAnomalyDisplayFields(df_Values, catalog)
    df_Values$CodebookCategory <- catalog$Category[match(df_Values$Variable, catalog$Variable)]
    df_Values$Analysis <- analysis
    df_Values$IssueID <- vapply(df_Values$.Row, function(row) paste(unique(object$Evidence$IssueID[
      object$Evidence$Analysis == analysis & object$Evidence$SourceRow == row & object$Evidence$Method == "isolation"]), collapse = "; "), character(1))
    df_Values$AdditionalEvidence <- vapply(seq_len(nrow(df_Values)), function(i) {
      e <- object$Evidence
      matches <- e$Analysis == analysis & e$Method != "isolation" &
        ((e$SourceRow == df_Values$.Row[i] & !is.na(e$Variable) & e$Variable == df_Values$Variable[i]) |
          (!is.na(e$OtherRow) & e$OtherRow == df_Values$.Row[i] & !is.na(e$OtherVariable) & e$OtherVariable == df_Values$Variable[i]))
      paste(unique(e$IssueID[which(matches)]), collapse = "; ")
    }, character(1))
    df_Summary <- source$ExplanationSummary
    index <- match(df_Values$.Row, df_Summary$.Row)
    for (column in c("Score", "Baseline", "SumContributions", "ReconstructionResidual", "Status")) df_Values[[column]] <- df_Summary[[column]][index]
    df_Values$ContributionID <- vapply(seq_len(nrow(df_Values)), function(i) ScidrReviewHash(list(analysis, df_Values$.Row[i], df_Values$Variable[i])), character(1))
    df_Values
  })
  empty <- ScidrAnomalyEmptyContributions()
  empty$Analysis <- empty$IssueID <- empty$AdditionalEvidence <- character()
  empty$Score <- empty$Baseline <- empty$SumContributions <- empty$ReconstructionResidual <- double()
  empty$Status <- empty$ContributionID <- character()
  output <- dplyr::bind_rows(c(list(empty), pieces))
  if (nrow(output) && "ExplanationMethod" %in% names(output) && all(!is.na(output$ExplanationMethod) & output$ExplanationMethod == "sensitivity")) {
    output <- output[setdiff(names(output), c("Contribution", "SimulationSD", "SimulationMin", "SimulationMax", "PositiveRank", "Baseline", "SumContributions", "ReconstructionResidual"))]
  }
  output
}
