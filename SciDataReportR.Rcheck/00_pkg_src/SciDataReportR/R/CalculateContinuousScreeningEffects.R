#' Calculate standardized continuous screening effects
#'
#' Fits one covariate-adjusted linear model per predictor using independently
#' selected complete cases. This shared engine keeps continuous screening
#' statistics aligned across profile and volcano visualizations.
#'
#' @param data A data frame.
#' @param predictor_vars Character vector of numeric predictor names.
#' @param outcome_var Character string naming a numeric outcome.
#' @param covariates Optional character vector of covariate names.
#' @param label_lookup Named character vector of display labels.
#'
#' @return A tibble containing one row per requested predictor, including
#'   diagnostics for models that could not be estimated.
#' @noRd
CalculateContinuousScreeningEffects <- function(data,
    predictor_vars,
    outcome_var,
    covariates = NULL,
    label_lookup = stats::setNames(predictor_vars, predictor_vars)) {
  purrr::map_dfr(predictor_vars, function(this_var) {
    model_vars <- unique(c(this_var, outcome_var, covariates))
    model_data <- data %>%
      dplyr::select(dplyr::all_of(model_vars)) %>%
      dplyr::filter(stats::complete.cases(.))
    final_n <- nrow(model_data)

    EmptyResult <- function(note) {
      tibble::tibble(
        Variable = this_var,
        Label = unname(label_lookup[[this_var]]),
        Beta = NA_real_,
        SE = NA_real_,
        CILow = NA_real_,
        CIHigh = NA_real_,
        PValue = NA_real_,
        N = final_n,
        R = NA_real_,
        AdjustedR = NA_real_,
        Note = note
      )
    }

    if (final_n < 3) {
      return(EmptyResult("Too few complete observations"))
    }

    predictor_sd <- stats::sd(model_data[[this_var]])
    outcome_sd <- stats::sd(model_data[[outcome_var]])
    if (!is.finite(predictor_sd) || predictor_sd == 0 ||
        !is.finite(outcome_sd) || outcome_sd == 0) {
      return(EmptyResult("Zero variance in predictor or outcome"))
    }

    model_data$.scidr_y_scaled <- as.numeric(scale(model_data[[outcome_var]]))
    model_data$.scidr_x_scaled <- as.numeric(scale(model_data[[this_var]]))
    model_formula <- stats::reformulate(
      c(".scidr_x_scaled", covariates),
      response = ".scidr_y_scaled"
    )

    model_error <- NULL
    model_fit <- tryCatch(
      stats::lm(model_formula, data = model_data),
      error = function(e) {
        model_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(model_fit)) {
      note <- if (!is.null(model_error) && grepl("contrasts", model_error, fixed = TRUE)) {
        "Covariate has fewer than two usable levels"
      } else {
        "Model failed"
      }
      return(EmptyResult(note))
    }

    coefficient_table <- summary(model_fit)$coefficients
    if (!".scidr_x_scaled" %in% rownames(coefficient_table) ||
        !is.finite(coefficient_table[".scidr_x_scaled", "Estimate"])) {
      return(EmptyResult("Predictor coefficient not estimable"))
    }

    df_resid <- stats::df.residual(model_fit)
    if (!is.finite(df_resid) || df_resid <= 0) {
      return(EmptyResult("Insufficient residual degrees of freedom"))
    }

    beta <- unname(coefficient_table[".scidr_x_scaled", "Estimate"])
    standard_error <- unname(coefficient_table[".scidr_x_scaled", "Std. Error"])
    p_value <- unname(coefficient_table[".scidr_x_scaled", "Pr(>|t|)"])
    critical_value <- stats::qt(0.975, df = df_resid)
    pearson_r <- tryCatch(
      stats::cor(model_data[[this_var]], model_data[[outcome_var]]),
      error = function(e) NA_real_
    )

    adjusted_r <- NA_real_
    if (length(covariates) > 0) {
      t_value <- unname(coefficient_table[".scidr_x_scaled", "t value"])
      if (is.finite(t_value)) {
        adjusted_r <- unname(t_value / sqrt(t_value^2 + df_resid))
      }
    }

    tibble::tibble(
      Variable = this_var,
      Label = unname(label_lookup[[this_var]]),
      Beta = beta,
      SE = standard_error,
      CILow = beta - critical_value * standard_error,
      CIHigh = beta + critical_value * standard_error,
      PValue = p_value,
      N = final_n,
      R = pearson_r,
      AdjustedR = adjusted_r,
      Note = NA_character_
    )
  })
}
