#' Apply a normative T-score model to new data
#'
#' Applies a previously fitted normative regression model to new data and
#' computes predicted values, z-scores, and T-scores using the same
#' preprocessing settings used during model development.
#'
#' This function is designed to work with the output of
#' `CreateNormativeTScoreModel()`. It uses the saved model and preprocessing
#' settings to score new observations consistently.
#'
#' @param data A data frame containing the test variable, count variable, and
#'   all predictors required by the normative model.
#' @param normative_obj A list returned by `CreateNormativeTScoreModel()`.
#' @param score_prefix A character string prefix used when naming output
#'   columns. Defaults to `"Norm"`.
#'
#' @return A tibble containing the original data plus scored columns, each
#'   named using the `score_prefix` string as a prefix:
#' \describe{
#'   \item{`{score_prefix}Raw`}{The raw input score.}
#'   \item{`{score_prefix}Scaled`}{The transformed analysis-scale score.}
#'   \item{`{score_prefix}Predicted`}{The predicted score from the normative model.}
#'   \item{`{score_prefix}Z`}{The z-score.}
#'   \item{`{score_prefix}T`}{The T-score.}
#' }
#'
#' @examples
#' # A reference group and a clinical group tested on Trail Making A
#' set.seed(206)
#' n_reference <- 220
#' n_clinical <- 80
#'
#' df <- tibble::tibble(
#'   Group = c(rep("Reference", n_reference), rep("Clinical", n_clinical)),
#'   Age = round(c(
#'     stats::rnorm(n_reference, 52, 12),
#'     stats::rnorm(n_clinical, 58, 12)
#'   )),
#'   Education = factor(sample(
#'     c("HighSchool", "College", "Graduate"), n_reference + n_clinical,
#'     replace = TRUE
#'   )),
#'   Sex = factor(sample(c("F", "M"), n_reference + n_clinical, replace = TRUE)),
#'   Visit = sample(1:3, n_reference + n_clinical, replace = TRUE)
#' )
#' df$TrailsA <- round(1000 * exp(stats::rnorm(
#'   nrow(df),
#'   mean = log(28) + 0.011 * (df$Age - 52) - 0.05 * (df$Visit - 1) +
#'     ifelse(df$Group == "Clinical", 0.35, 0),
#'   sd = 0.22
#' )))
#'
#' # Fit the norms on the reference group only
#' norm_obj <- CreateNormativeTScoreModel(
#'   data = df,
#'   test_var = "TrailsA",
#'   count_var = "Visit",
#'   covariates = c("Age", "Education", "Sex"),
#'   reference_var = "Group",
#'   reference_value = "Reference",
#'   include_practice_effect = TRUE,
#'   reverse_score = TRUE,
#'   convert_seconds = TRUE,
#'   log_transform = TRUE,
#'   return_plots = FALSE
#' )
#'
#' # Score everyone through the same model
#' scored_df <- ApplyNormativeTScores(
#'   data = df,
#'   normative_obj = norm_obj
#' )
#'
#' # Before: raw completion times, grouped by clinical status
#' attr(scored_df$NormRaw, "label") <- "Trail Making A completion time (ms)"
#' PlotContinuousDistributions(
#'   scored_df, variables = "NormRaw", Fill = "Group", ncol = 1
#' )
#'
#' # After: demographically adjusted T-scores. The reference group is centered
#' # near 50, and the clinical shift is now directly interpretable in SD units.
#' attr(scored_df$NormT, "label") <- "Demographically adjusted Trail Making A T-score"
#' PlotContinuousDistributions(
#'   scored_df, variables = "NormT", Fill = "Group", ncol = 1
#' )
#'
#' @section What the example shows:
#' The example builds a cohort with a healthy reference group and a clinical
#' group tested on Trail Making A. Times are recorded in milliseconds, get
#' slower with age, and improve a little each time the test is repeated. The
#' norms are fitted on the reference group only, adjusting for demographics and
#' for practice across visits, and then everyone is scored through that same
#' model.
#'
#' The two density plots are the payoff. Raw completion times are right-skewed
#' and in test-specific units, and whether 40 seconds is unusual depends on who
#' the participant is, so the groups are hard to compare by eye. T-scores put
#' everyone on one interpretable scale: the reference group is centered at 50
#' with a standard deviation of 10 by construction, so the clinical group's
#' shift can be read straight off the axis - here about 1.7 SD below
#' expectation for their age, education, sex, and visit number.
#'
#' @param df \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
ApplyNormativeTScores <- function(data,
    normative_obj,
    score_prefix = "Norm",
    df = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(df)) {
    lifecycle::deprecate_warn("19.15.0", "ApplyNormativeTScores(df)", "ApplyNormativeTScores(data)")
    data <- df
  }
  if (!missing(data)) df <- data


  # Validate inputs

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }

  if (!is.list(normative_obj)) {
    stop("`normative_obj` must be a list returned by `CreateNormativeTScoreModel()`.")
  }

  if (is.null(normative_obj$model) || !inherits(normative_obj$model, "lm")) {
    stop("`normative_obj$model` must be an `lm` object.")
  }

  if (is.null(normative_obj$settings) || !is.list(normative_obj$settings)) {
    stop("`normative_obj$settings` must be present and must be a list.")
  }

  if (!is.character(score_prefix) || length(score_prefix) != 1 || nchar(score_prefix) == 0) {
    stop("`score_prefix` must be a single non-empty character string.")
  }

  settings <- normative_obj$settings
  model <- normative_obj$model
  sigma_resid <- summary(model)$sigma

  required_settings <- c(
    "test_var",
    "count_var",
    "covariates",
    "include_practice_effect",
    "reverse_score",
    "convert_seconds",
    "seconds_divisor",
    "log_transform"
  )

  missing_settings <- setdiff(required_settings, names(settings))

  if (length(missing_settings) > 0) {
    stop(
      paste0(
        "`normative_obj$settings` is missing required elements: ",
        paste(missing_settings, collapse = ", ")
      )
    )
  }

  model_terms <- settings$covariates

  if (isTRUE(settings$include_practice_effect)) {
    model_terms <- c(model_terms, settings$count_var)
  }

  required_cols <- unique(c(settings$test_var, model_terms))
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in `df`: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  raw_col <- paste0(score_prefix, "Raw")
  scaled_col <- paste0(score_prefix, "Scaled")
  predicted_col <- paste0(score_prefix, "Predicted")
  z_col <- paste0(score_prefix, "Z")
  t_col <- paste0(score_prefix, "T")

  # Prepare data

  df_scored <- df %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      !!raw_col := .data[[settings$test_var]]
    )

  if (isTRUE(settings$convert_seconds)) {
    df_scored <- df_scored %>%
      dplyr::mutate(
        !!scaled_col := .data[[raw_col]] / settings$seconds_divisor
      )
  } else {
    df_scored <- df_scored %>%
      dplyr::mutate(
        !!scaled_col := .data[[raw_col]]
      )
  }

  if (isTRUE(settings$log_transform)) {
    non_positive_n <- sum(df_scored[[scaled_col]] <= 0, na.rm = TRUE)

    if (non_positive_n > 0) {
      warning(
        paste0(
          "Found ", non_positive_n,
          " values <= 0 in the analysis variable. ",
          "These values will become non-finite after log10 transformation ",
          "and will receive missing scored outputs."
        )
      )
    }

    df_scored <- df_scored %>%
      dplyr::mutate(
        !!scaled_col := log10(.data[[scaled_col]])
      )
  }

  if (isTRUE(settings$reverse_score)) {
    df_scored <- df_scored %>%
      dplyr::mutate(
        !!scaled_col := -1 * .data[[scaled_col]]
      )
  }

  valid_prediction_rows <- stats::complete.cases(
    df_scored[, model_terms, drop = FALSE]
  )

  predicted <- rep(NA_real_, nrow(df_scored))
  predicted[valid_prediction_rows] <- stats::predict(
    model,
    newdata = df_scored[valid_prediction_rows, , drop = FALSE]
  )

  # Return result

  df_scored %>%
    dplyr::mutate(
      !!predicted_col := predicted,
      !!z_col := dplyr::if_else(
        is.finite(.data[[scaled_col]]) & is.finite(.data[[predicted_col]]),
        (.data[[scaled_col]] - .data[[predicted_col]]) / sigma_resid,
        NA_real_
      ),
      !!t_col := dplyr::if_else(
        is.finite(.data[[z_col]]),
        (.data[[z_col]] * 10) + 50,
        NA_real_
      )
    )
}
