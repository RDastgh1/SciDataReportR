#' Create normative T-scores from a regression model
#'
#' Fits a normative regression model in a user-defined reference subgroup and
#' uses the model residual standard deviation to convert observed scores into
#' z-scores and T-scores. This is useful for creating demographically adjusted
#' cognitive norms with optional practice effect adjustment and optional
#' preprocessing such as unit conversion, log transformation, and reverse
#' scoring.
#'
#' The reference group is defined by `reference_var == reference_value`.
#' When `include_practice_effect = FALSE`, the normative model is fit only on
#' rows where `count_var == baseline_count_value`. When
#' `include_practice_effect = TRUE`, `count_var` is added to the model and all
#' available visits in the reference group are used.
#'
#' @param df A data frame containing the test variable, count variable,
#'   reference group variable, and covariates.
#' @param test_var A character string naming the raw test score variable.
#' @param count_var A character string naming the visit count or practice count
#'   variable.
#' @param covariates A character vector of covariate column names to include in
#'   the normative model.
#' @param reference_var A character string naming the variable used to define
#'   the normative reference group.
#' @param reference_value The value of `reference_var` that defines the
#'   normative reference group.
#' @param include_practice_effect Logical. If `TRUE`, `count_var` is included as
#'   a predictor and the model is fit using all available visits in the
#'   reference group. If `FALSE`, the model is fit only on rows where
#'   `count_var == baseline_count_value`.
#' @param baseline_count_value The value of `count_var` used to define the
#'   baseline visit when `include_practice_effect = FALSE`. Defaults to `1`.
#' @param reverse_score Logical. If `TRUE`, the analysis-scale score is
#'   multiplied by `-1` so that higher values reflect better performance.
#' @param convert_seconds Logical. If `TRUE`, the raw score is divided by
#'   `seconds_divisor` before further processing.
#' @param seconds_divisor Numeric divisor used when `convert_seconds = TRUE`.
#'   Defaults to `1000`.
#' @param log_transform Logical. If `TRUE`, applies `log10()` to the analysis
#'   score after optional unit conversion and before optional reverse scoring.
#' @param codebook Optional data frame with columns `Variable` and `Label`.
#'   If supplied, plot labels use variable labels when available.
#' @param return_plots Logical. If `TRUE`, returns a list of diagnostic plots.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{data}{A tibble containing the original data plus `NormRaw`,
#'   `NormScaled`, `NormPredicted`, `NormZ`, and `NormT`.}
#'   \item{model}{The fitted `lm` object.}
#'   \item{model_summary}{A tibble of coefficient estimates.}
#'   \item{model_fit}{A one-row tibble containing model fit statistics.}
#'   \item{training_data}{The rows used to fit the normative model.}
#'   \item{plots}{A named list of ggplot objects when `return_plots = TRUE`.}
#'   \item{settings}{A list of preprocessing and modeling settings used.}
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
#' out <- CreateNormativeTScores(
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
#' out$data
#' out$model
#' @export
CreateNormativeTScores <- function(df,
                                   test_var,
                                   count_var,
                                   covariates,
                                   reference_var,
                                   reference_value,
                                   include_practice_effect = FALSE,
                                   baseline_count_value = 1,
                                   reverse_score = FALSE,
                                   convert_seconds = FALSE,
                                   seconds_divisor = 1000,
                                   log_transform = TRUE,
                                   codebook = NULL,
                                   return_plots = TRUE) {

  # Validate inputs

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }

  if (!is.character(test_var) || length(test_var) != 1) {
    stop("`test_var` must be a single character string.")
  }

  if (!is.character(count_var) || length(count_var) != 1) {
    stop("`count_var` must be a single character string.")
  }

  if (!is.character(reference_var) || length(reference_var) != 1) {
    stop("`reference_var` must be a single character string.")
  }

  if (!is.character(covariates) || length(covariates) < 1) {
    stop("`covariates` must be a character vector with at least one variable.")
  }

  required_cols <- unique(c(test_var, count_var, reference_var, covariates))
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in `df`: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (!is.logical(include_practice_effect) || length(include_practice_effect) != 1) {
    stop("`include_practice_effect` must be TRUE or FALSE.")
  }

  if (!is.logical(reverse_score) || length(reverse_score) != 1) {
    stop("`reverse_score` must be TRUE or FALSE.")
  }

  if (!is.logical(convert_seconds) || length(convert_seconds) != 1) {
    stop("`convert_seconds` must be TRUE or FALSE.")
  }

  if (!is.logical(log_transform) || length(log_transform) != 1) {
    stop("`log_transform` must be TRUE or FALSE.")
  }

  if (!is.logical(return_plots) || length(return_plots) != 1) {
    stop("`return_plots` must be TRUE or FALSE.")
  }

  if (!is.numeric(seconds_divisor) || length(seconds_divisor) != 1 || seconds_divisor <= 0) {
    stop("`seconds_divisor` must be a single positive number.")
  }

  if (!is.null(codebook)) {
    if (!is.data.frame(codebook) || !all(c("Variable", "Label") %in% names(codebook))) {
      stop("`codebook` must be a data frame with columns `Variable` and `Label`.")
    }
  }

  # Prepare data

  df_prepped <- df %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      NormRaw = .data[[test_var]],
      NormCount = .data[[count_var]]
    )

  if (convert_seconds) {
    df_prepped <- df_prepped %>%
      dplyr::mutate(NormScaled = NormRaw / seconds_divisor)
  } else {
    df_prepped <- df_prepped %>%
      dplyr::mutate(NormScaled = NormRaw)
  }

  if (log_transform) {
    non_positive_n <- sum(df_prepped$NormScaled <= 0, na.rm = TRUE)

    if (non_positive_n > 0) {
      warning(
        paste0(
          "Found ", non_positive_n,
          " values <= 0 in the analysis variable. ",
          "These values will become non-finite after log10 transformation ",
          "and will be excluded from model fitting and scoring where applicable."
        )
      )
    }

    df_prepped <- df_prepped %>%
      dplyr::mutate(NormScaled = log10(NormScaled))
  }

  if (reverse_score) {
    df_prepped <- df_prepped %>%
      dplyr::mutate(NormScaled = -1 * NormScaled)
  }

  model_terms <- covariates

  if (include_practice_effect) {
    model_terms <- c(model_terms, count_var)
  }

  model_formula <- stats::reformulate(
    termlabels = model_terms,
    response = "NormScaled"
  )

  training_data <- df_prepped %>%
    dplyr::filter(.data[[reference_var]] == reference_value)

  if (!include_practice_effect) {
    training_data <- training_data %>%
      dplyr::filter(.data[[count_var]] == baseline_count_value)
  }

  training_data <- training_data %>%
    dplyr::filter(
      stats::complete.cases(
        dplyr::select(., dplyr::all_of(c("NormScaled", model_terms)))
      )
    ) %>%
    dplyr::filter(is.finite(.data$NormScaled))

  if (nrow(training_data) == 0) {
    stop("No complete cases were available for model fitting in the reference sample.")
  }

  if (nrow(training_data) <= length(model_terms) + 1) {
    stop(
      paste0(
        "Not enough complete cases to fit the model. ",
        "The reference sample must have more rows than predictors."
      )
    )
  }

  # Apply labels

  get_label <- function(var_name, codebook) {
    if (is.null(codebook)) {
      return(var_name)
    }

    label_value <- codebook %>%
      dplyr::filter(.data$Variable == var_name) %>%
      dplyr::pull(.data$Label)

    if (length(label_value) == 0 || is.na(label_value[1]) || label_value[1] == "") {
      return(var_name)
    }

    label_value[1]
  }

  test_label <- get_label(test_var, codebook)
  reference_label <- get_label(reference_var, codebook)
  count_label <- get_label(count_var, codebook)

  # Build outputs

  model <- stats::lm(
    formula = model_formula,
    data = training_data
  )

  sigma_resid <- summary(model)$sigma

  predicted <- rep(NA_real_, nrow(df_prepped))
  valid_prediction_rows <- stats::complete.cases(
    df_prepped[, unique(model_terms), drop = FALSE]
  )

  predicted[valid_prediction_rows] <- stats::predict(
    model,
    newdata = df_prepped[valid_prediction_rows, , drop = FALSE]
  )

  df_scored <- df_prepped %>%
    dplyr::mutate(
      NormPredicted = predicted,
      NormZ = dplyr::if_else(
        is.finite(NormScaled) & is.finite(NormPredicted),
        (NormScaled - NormPredicted) / sigma_resid,
        NA_real_
      ),
      NormT = dplyr::if_else(
        is.finite(NormZ),
        (NormZ * 10) + 50,
        NA_real_
      )
    )

  model_summary <- broom::tidy(model) %>%
    dplyr::as_tibble()

  model_fit <- broom::glance(model) %>%
    dplyr::transmute(
      n_train = stats::nobs(model),
      sigma = .data$sigma,
      r_squared = .data$r.squared,
      adj_r_squared = .data$adj.r.squared,
      aic = .data$AIC,
      bic = .data$BIC
    )

  plots <- NULL

  if (isTRUE(return_plots)) {
    plot_data_v1 <- df_scored %>%
      dplyr::filter(.data[[count_var]] == baseline_count_value)

    p_raw <- ggplot2::ggplot(
      plot_data_v1,
      ggplot2::aes(x = NormRaw, fill = as.factor(.data[[reference_var]]))
    ) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = test_label,
        y = "Count",
        title = "Raw score distribution",
        fill = reference_label
      )

    p_scaled <- ggplot2::ggplot(
      plot_data_v1,
      ggplot2::aes(x = NormScaled, fill = as.factor(.data[[reference_var]]))
    ) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Analysis scale",
        y = "Count",
        title = "Scaled score distribution",
        fill = reference_label
      )

    p_tscore <- ggplot2::ggplot(
      plot_data_v1,
      ggplot2::aes(x = NormT, fill = as.factor(.data[[reference_var]]))
    ) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity") +
      ggplot2::geom_vline(xintercept = 50, linetype = "dashed") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "T-score",
        y = "Count",
        title = "T-score distribution",
        fill = reference_label
      )

    p_reference <- ggplot2::ggplot(
      plot_data_v1,
      ggplot2::aes(
        x = as.factor(.data[[reference_var]]),
        y = NormT,
        fill = as.factor(.data[[reference_var]])
      )
    ) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = reference_label,
        y = "T-score",
        title = "T-scores by reference group"
      )

    p_practice <- ggplot2::ggplot(
      df_scored,
      ggplot2::aes(
        x = .data[[count_var]],
        y = NormT,
        color = as.factor(.data[[reference_var]])
      )
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", se = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = count_label,
        y = "T-score",
        title = "T-scores by visit count",
        color = reference_label
      )

    covariate_plots <- purrr::map(
      covariates,
      function(var_name) {
        x_vec <- df_scored[[var_name]]
        var_label <- get_label(var_name, codebook)

        if (is.numeric(x_vec)) {
          ggplot2::ggplot(
            plot_data_v1,
            ggplot2::aes(
              x = .data[[var_name]],
              y = NormT,
              color = as.factor(.data[[reference_var]])
            )
          ) +
            ggplot2::geom_point(alpha = 0.7) +
            ggplot2::geom_smooth(method = "lm", se = FALSE) +
            ggplot2::theme_bw() +
            ggplot2::labs(
              x = var_label,
              y = "T-score",
              title = paste("T-scores by", var_label),
              color = reference_label
            )
        } else {
          ggplot2::ggplot(
            plot_data_v1,
            ggplot2::aes(
              x = as.factor(.data[[var_name]]),
              y = NormT,
              fill = as.factor(.data[[reference_var]])
            )
          ) +
            ggplot2::geom_boxplot(alpha = 0.7) +
            ggplot2::theme_bw() +
            ggplot2::labs(
              x = var_label,
              y = "T-score",
              title = paste("T-scores by", var_label),
              fill = reference_label
            )
        }
      }
    )

    names(covariate_plots) <- covariates

    plots <- c(
      list(
        raw = p_raw,
        scaled = p_scaled,
        tscore = p_tscore,
        reference = p_reference,
        practice = p_practice
      ),
      covariate_plots
    )
  }

  # Return result

  list(
    data = df_scored,
    model = model,
    model_summary = model_summary,
    model_fit = model_fit,
    training_data = training_data,
    plots = plots,
    settings = list(
      test_var = test_var,
      count_var = count_var,
      covariates = covariates,
      reference_var = reference_var,
      reference_value = reference_value,
      include_practice_effect = include_practice_effect,
      baseline_count_value = baseline_count_value,
      reverse_score = reverse_score,
      convert_seconds = convert_seconds,
      seconds_divisor = seconds_divisor,
      log_transform = log_transform,
      model_formula = model_formula
    )
  )
}
