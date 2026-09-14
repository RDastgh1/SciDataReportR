#' Assess missingness in a wide feature matrix
#'
#' @param data A data frame with samples in rows.
#' @param variables Character vector of numeric feature columns.
#' @param missing_threshold Numeric proportion above which a feature is flagged.
#' @param intensity_heuristic Logical. Whether to calculate the exploratory
#'   missingness-versus-sample-intensity correlation.
#' @param intensity_cutoff Minimum absolute correlation for an
#'   intensity-dependent flag.
#'
#' @return A list with a feature-level `summary` tibble and `plot`.
#' @export
AssessFeatureMissingness <- function(
    data,
    variables = NULL,
    missing_threshold = 0.5,
    intensity_heuristic = TRUE,
    intensity_cutoff = 0.3
) {
  variables <- ValidateFeatureMatrix(data, variables)
  if (!is.numeric(missing_threshold) || length(missing_threshold) != 1 ||
      is.na(missing_threshold) || missing_threshold < 0 || missing_threshold > 1) {
    stop("`missing_threshold` must be a single value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(intensity_cutoff) || length(intensity_cutoff) != 1 ||
      is.na(intensity_cutoff) || intensity_cutoff < 0 || intensity_cutoff > 1) {
    stop("`intensity_cutoff` must be a single value between 0 and 1.", call. = FALSE)
  }

  feature_mat <- as.matrix(data[variables])
  storage.mode(feature_mat) <- "double"
  na_indicator <- is.na(feature_mat)
  pct_missing <- colMeans(na_indicator) * 100
  n_samples <- nrow(feature_mat)
  intensity_correlation <- rep(NA_real_, length(variables))
  missingness_type <- rep(NA_character_, length(variables))
  note <- rep(NA_character_, length(variables))

  for (i in seq_along(variables)) {
    n_missing <- sum(na_indicator[, i])
    if (n_missing == 0) {
      missingness_type[i] <- "None"
      next
    }
    if (n_missing == n_samples) {
      note[i] <- "All values missing; mechanism cannot be assessed."
      next
    }
    if (!isTRUE(intensity_heuristic) || ncol(feature_mat) < 2) {
      next
    }

    other_intensity <- apply(feature_mat[, -i, drop = FALSE], 1, stats::median, na.rm = TRUE)
    indicator <- as.numeric(na_indicator[, i])
    if (stats::sd(indicator) == 0 || stats::sd(other_intensity, na.rm = TRUE) == 0 ||
        sum(is.finite(other_intensity)) < 3) {
      note[i] <- "Insufficient variation to compute the intensity heuristic."
      next
    }
    this_correlation <- suppressWarnings(stats::cor(
      indicator, other_intensity, method = "spearman", use = "pairwise.complete.obs"
    ))
    intensity_correlation[i] <- this_correlation
    missingness_type[i] <- if (!is.na(this_correlation) && this_correlation <= -abs(intensity_cutoff)) {
      "Intensity-dependent"
    } else {
      "Not intensity-dependent"
    }
  }

  list(
    summary = tibble::tibble(
      Variable = variables,
      PctMissing = round(unname(pct_missing), 3),
      IntensityCorrelation = round(intensity_correlation, 3),
      MissingnessType = missingness_type,
      Flagged = unname(pct_missing / 100 > missing_threshold),
      Note = note
    ),
    plot = PlotMissingData(data = data, variables = variables)
  )
}
