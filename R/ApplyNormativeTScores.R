#' Apply a normative T-score model to new data
#'
#' Applies a previously fitted normative regression model to new data and
#' computes predicted values, z-scores, and T-scores using the same
#' preprocessing settings used during model development.
#'
#' This function is designed to work with the output of
#' `CreateNormativeTScores()`. It uses the saved model and preprocessing
#' settings to score new observations consistently.
#'
#' @param df A data frame containing the test variable, count variable, and
#'   all predictors required by the normative model.
#' @param normative_obj A list returned by `CreateNormativeTScores()`.
#' @param score_prefix A character string prefix used when naming output
#'   columns. Defaults to `"Norm"`.
#'
#' @return A tibble containing the original data plus scored columns:
#' \describe{
#'   \item{<prefix>Raw}{The raw input score.}
#'   \item{<prefix>Scaled}{The transformed analysis-scale score.}
#'   \item{<prefix>Predicted}{The predicted score from the normative model.}
#'   \item{<prefix>Z}{The z-score.}
#'   \item{<prefix>T}{The T-score.}
#' }
#'
#' @examples
#' df <- tibble::tibble(
#'   Group = c(
#'     rep("Reference", 8),
#'     rep("Clinical", 2)
#'   ),
#'   Age = c(30, 34, 38, 42, 46, 50, 54, 58, 40, 52),
#'   Education = factor(c(
#'     "College", "College", "Graduate", "Graduate",
#'     "College", "Graduate", "College", "Graduate",
#'     "College", "Graduate"
#'   )),
#'   Sex = factor(c(
#'     "F", "M", "F", "M", "F", "M", "F", "M", "F", "M"
#'   )),
#'   Visit = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2),
#'   TrailsA = c(35, 38, 40, 43, 36, 39, 41, 44, 47, 49) * 1000
#' )
#'
#' norm_obj <- CreateNormativeTScores(
#'   df = df,
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
#' scored_df <- ApplyNormativeTScores(
#'   df = df,
#'   normative_obj = norm_obj
#' )
#'
#' @export
ApplyNormativeTScores <- function(df,
                                  normative_obj,
                                  score_prefix = "Norm") {

  # Validate inputs

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }

  if (!is.list(normative_obj)) {
    stop("`normative_obj` must be a list returned by `CreateNormativeTScores()`.")
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
