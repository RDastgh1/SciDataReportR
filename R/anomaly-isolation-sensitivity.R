ScidrAnomalyEmptySensitivity <- function() tibble::tibble(.Row = integer(), Stratum = integer(),
  Variable = character(), Label = character(), Observed = character(), ScoreDecrease = double(),
  ReplacementScoreMean = double(), ReplacementScoreMin = double(), ReplacementScoreMax = double(),
  ReplacementN = integer(), PositiveRank = double())

ScidrAnomalyIsolationSensitivity <- function(fit, df_X, rows, scores, stratum, setup, args) {
  df_X <- as.data.frame(df_X, optional = TRUE)
  restore <- ScidrAnomalySeed(setup$seed)
  on.exit(restore(), add = TRUE)
  started <- proc.time()[["elapsed"]]
  Seed <- function(offset) as.integer((as.double(setup$seed) + 1000003 * as.double(stratum) + offset) %% (.Machine$integer.max - 1))
  eligible <- which(is.finite(scores))
  selected <- if (is.null(args$rows)) head(eligible[order(-scores[eligible], rows[eligible])], args$top_n) else eligible[rows[eligible] %in% args$rows]
  df_Summary <- tibble::tibble(.Row = rows, Stratum = as.integer(stratum), Score = scores,
    Baseline = NA_real_, SumContributions = NA_real_, ReconstructionResidual = NA_real_,
    ExplainedVariables = 0L, Status = ifelse(is.finite(scores), "Not requested", "Unavailable"),
    Reason = ifelse(is.finite(scores), "Outside explanation selection", "Native score unavailable"))
  # Limit hybrid batches to approximately one million input cells as well as columns.
  batch_size <- max(1L, min(args$batch_size, floor(1e6 / (ncol(df_X) * min(args$replacements, args$background_n, length(eligible))))))
  metadata <- list(Method = "sensitivity", SelectedRows = rows[selected], Settings = args,
    SourceVariables = names(df_X), ModelingVariables = if (is.null(setup$time_adapter)) names(df_X) else setup$time_adapter$ModelFeatures, InputClasses = lapply(df_X, class),
    FactorLevels = lapply(df_X[vapply(df_X, is.factor, logical(1))], levels),
    DateConversions = setup$date_conversions, TimeConversions = setup$time_metadata, ModelAdapter = setup$time_adapter, PredictionThreads = 1L,
    BackgroundRows = integer(), ReferenceRows = integer(), BackgroundSeed = Seed(0),
    ReferenceSeed = Seed(1), BatchSize = batch_size)
  pieces <- list(ScidrAnomalyEmptySensitivity())
  if (length(selected)) {
    message("Screening isolation columns: ", length(selected), " records x ", ncol(df_X), " features x up to ",
      args$replacements, " replacements (stratum ", stratum, ").")
    df_Summary$Status[selected] <- "Unavailable"
    df_Summary$Reason[selected] <- "Screening did not complete"
    failure <- tryCatch({
      set.seed(metadata$BackgroundSeed)
      background <- eligible[sample.int(length(eligible), min(args$background_n, length(eligible))) ]
      set.seed(metadata$ReferenceSeed)
      reference <- background[sample.int(length(background), min(args$replacements, length(background)))]
      metadata$BackgroundRows <- rows[background]
      metadata$ReferenceRows <- rows[reference]
      count <- length(reference)
      Predict <- function(newdata) as.numeric(stats::predict(fit, ScidrTimePredictData(newdata, setup$time_adapter), type = "score", nthreads = 1L))
      for (i in selected) {
        message("Screening source row ", rows[i], " (", match(i, selected), "/", length(selected), ").")
        row_failure <- tryCatch({
          means <- lower <- upper <- rep(NA_real_, ncol(df_X))
          for (start in seq.int(1L, ncol(df_X), by = batch_size)) {
            features <- seq.int(start, min(ncol(df_X), start + batch_size - 1L))
            df_Hybrid <- df_X[rep(i, length(features) * count), , drop = FALSE]
            for (j in seq_along(features)) {
              positions <- (j - 1L) * count + seq_len(count)
              df_Hybrid[[features[j]]][positions] <- df_X[[features[j]]][reference]
            }
            predicted <- Predict(df_Hybrid)
            if (length(predicted) != count * length(features) || any(!is.finite(predicted))) stop("Nonfinite or incomplete replacement scores.")
            values <- matrix(predicted, nrow = count)
            means[features] <- colMeans(values)
            lower[features] <- apply(values, 2, min)
            upper[features] <- apply(values, 2, max)
          }
          effect <- scores[i] - means
          pieces[[length(pieces) + 1L]] <- tibble::tibble(.Row = rows[i], Stratum = as.integer(stratum),
            Variable = names(df_X), Label = unname(setup$labels[names(df_X)]),
            Observed = unname(vapply(names(df_X), function(variable) ScidrTimeObserved(setup$data, variable, rows[i], setup$time_metadata), character(1))),
            ScoreDecrease = effect, ReplacementScoreMean = means, ReplacementScoreMin = lower,
            ReplacementScoreMax = upper, ReplacementN = as.integer(count),
            PositiveRank = ifelse(effect > 0, rank(-replace(effect, effect <= 0, NA_real_), ties.method = "min", na.last = "keep"), NA_real_))
          df_Summary$ExplainedVariables[i] <- ncol(df_X)
          df_Summary$Status[i] <- "Explained"
          df_Summary$Reason[i] <- ""
          NULL
        }, error = function(e) conditionMessage(e))
        if (!is.null(row_failure)) {
          df_Summary$Reason[i] <- row_failure
          setup$Add("isolation", stratum, "explanation_failure", row_failure, rows[i])
        }
      }
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(failure)) {
      pending <- selected[df_Summary$Status[selected] != "Explained"]
      df_Summary$Reason[pending] <- failure
      setup$Add("isolation", stratum, "explanation_failure", failure, rows[pending])
    }
  }
  metadata$ElapsedSeconds <- unname(proc.time()[["elapsed"]] - started)
  if (length(selected)) message("Isolation column screening finished: ", sum(df_Summary$Status == "Explained"), "/", length(selected),
    " records in ", round(metadata$ElapsedSeconds, 2), " seconds.")
  list(Contributions = ScidrAnomalyEmptyContributions(), Sensitivity = dplyr::bind_rows(pieces), Summary = df_Summary, Metadata = metadata)
}
