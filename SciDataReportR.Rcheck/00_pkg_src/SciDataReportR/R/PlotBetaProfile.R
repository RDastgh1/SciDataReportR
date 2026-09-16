#' Plot standardized beta profiles across continuous predictors
#'
#' Fits a separate linear model for every predictor and displays the resulting
#' standardized regression coefficients with 95% confidence intervals. This is
#' the continuous-outcome companion to [PlotZScore()]: use it to compare the
#' direction and magnitude of associations across a panel of predictors.
#'
#' Each model has the form `scale(outcome) ~ scale(predictor) + covariates`.
#' Covariates retain their original representation and are not standardized.
#' The plotted beta is therefore the standard-deviation change in the outcome
#' associated with a one-standard-deviation increase in the predictor,
#' conditional on the supplied covariates. Complete cases are selected
#' independently for every predictor model.
#'
#' @param data A data frame containing the outcome, predictors, and covariates.
#' @param predictor_vars Nonempty character vector of numeric predictor names.
#' @param outcome_var Character string naming a numeric continuous outcome.
#' @param covariates Optional character vector of covariate names. Covariates
#'   are included without automatic standardization.
#' @param VariableCategories Optional categories for the predictors. Supply a
#'   vector corresponding to `predictor_vars`, a named vector keyed by predictor
#'   name, or a data frame with `Variable` and `Category` columns.
#' @param Sort Variable ordering. One of `"original"`, `"pvalue"`, `"fdr"`,
#'   `"effect"`, `"within_category_pvalue"`, or
#'   `"within_category_effect"`.
#' @param AdjustMethod Multiple-testing method passed to [stats::p.adjust()].
#' @param Alpha Numeric significance threshold recorded in the returned
#'   metadata. Significance does not control plot color.
#' @param Relabel Logical. If `TRUE`, display labels are resolved from the
#'   supplied codebook, then variable label attributes, with variable names as
#'   the fallback.
#' @param codebook Optional data frame containing `Variable` and `Label`.
#' @param RemoveXAxisLabels Logical. If `TRUE`, x-axis labels are hidden.
#' @param InteractiveLabels Logical. If `TRUE`, the point layer contains a
#'   `text` aesthetic for `plotly::ggplotly(..., tooltip = "text")`.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`Plot`}{A ggplot object showing standardized betas and 95% CIs.}
#'     \item{`ResultsTable`}{A tibble with one successfully analyzed predictor
#'       per row and columns `Variable`, `Label`, `Category`, `Beta`, `SE`,
#'       `CILow`, `CIHigh`, `PValue`, `FDR`, `N`, `R`, `AdjustedR`, and
#'       `Tooltip`.}
#'     \item{`Metadata`}{A named list describing the outcome, requested and
#'       analyzed predictors, covariates, adjustment method, alpha threshold,
#'       sorting mode, confidence level, and number of fitted models.}
#'   }
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' predictors <- c("Adiponectin", "C_Reactive_Protein", "Ferritin", "tau")
#'
#' # Unadjusted
#' unadjusted <- PlotBetaProfile(
#'   data = Labelled,
#'   predictor_vars = predictors,
#'   outcome_var = "AXL"
#' )
#' unadjusted$Plot
#'
#' # Adjusted
#' adjusted <- PlotBetaProfile(
#'   data = Labelled,
#'   predictor_vars = predictors,
#'   outcome_var = "AXL",
#'   covariates = c("age", "sex")
#' )
#' adjusted$ResultsTable
#'
#' # Categories and within-category sorting
#' categories <- c("Metabolic", "Inflammation", "Inflammation", "Neurology")
#' categorized <- PlotBetaProfile(
#'   data = Labelled,
#'   predictor_vars = predictors,
#'   outcome_var = "AXL",
#'   covariates = c("age", "sex"),
#'   VariableCategories = categories,
#'   Sort = "within_category_pvalue"
#' )
#' categorized$Plot
#' @export
PlotBetaProfile <- function(data,
    predictor_vars,
    outcome_var,
    covariates = NULL,
    VariableCategories = NULL,
    Sort = c(
      "original",
      "pvalue",
      "fdr",
      "effect",
      "within_category_pvalue",
      "within_category_effect"
    ),
    AdjustMethod = "fdr",
    Alpha = 0.05,
    Relabel = TRUE,
    codebook = NULL,
    RemoveXAxisLabels = TRUE,
    InteractiveLabels = TRUE) {
  # Validate inputs

  Sort <- match.arg(Sort)
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (!is.character(predictor_vars) || length(predictor_vars) == 0 ||
      anyNA(predictor_vars) || any(!nzchar(predictor_vars))) {
    stop("predictor_vars must be a nonempty character vector of variable names.", call. = FALSE)
  }
  if (!is.character(outcome_var) || length(outcome_var) != 1 || is.na(outcome_var) ||
      !nzchar(outcome_var)) {
    stop("outcome_var must be one nonempty character string.", call. = FALSE)
  }
  if (!outcome_var %in% names(data)) {
    stop("outcome_var was not found in data: ", outcome_var, ".", call. = FALSE)
  }
  if (!is.numeric(data[[outcome_var]])) {
    stop("outcome_var must be numeric and continuous.", call. = FALSE)
  }
  if (!is.null(covariates) &&
      (!is.character(covariates) || anyNA(covariates) || any(!nzchar(covariates)))) {
    stop("covariates must be NULL or a character vector of variable names.", call. = FALSE)
  }
  missing_covariates <- setdiff(covariates, names(data))
  if (length(missing_covariates) > 0) {
    stop(
      "The following covariates were not found in data: ",
      paste(missing_covariates, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  covariates <- unique(covariates)
  if (!is.character(AdjustMethod) || length(AdjustMethod) != 1 ||
      !AdjustMethod %in% stats::p.adjust.methods) {
    stop(
      "AdjustMethod must be one of: ",
      paste(stats::p.adjust.methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.numeric(Alpha) || length(Alpha) != 1 || is.na(Alpha) ||
      Alpha <= 0 || Alpha >= 1) {
    stop("Alpha must be one number greater than 0 and less than 1.", call. = FALSE)
  }
  logical_args <- list(
    Relabel = Relabel,
    RemoveXAxisLabels = RemoveXAxisLabels,
    InteractiveLabels = InteractiveLabels
  )
  invalid_logicals <- names(logical_args)[!vapply(
    logical_args,
    function(x) is.logical(x) && length(x) == 1 && !is.na(x),
    logical(1)
  )]
  if (length(invalid_logicals) > 0) {
    stop(paste(invalid_logicals, collapse = ", "), " must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(codebook)) {
    if (!is.data.frame(codebook)) {
      stop("codebook must be a data frame when supplied.", call. = FALSE)
    }
    missing_codebook_columns <- setdiff(c("Variable", "Label"), names(codebook))
    if (length(missing_codebook_columns) > 0) {
      stop("codebook must contain columns Variable and Label.", call. = FALSE)
    }
  }

  predictors_requested <- predictor_vars
  category_lookup <- stats::setNames(rep(NA_character_, length(predictor_vars)), predictor_vars)
  category_order <- character(0)

  if (!is.null(VariableCategories)) {
    if (is.data.frame(VariableCategories)) {
      missing_category_columns <- setdiff(c("Variable", "Category"), names(VariableCategories))
      if (length(missing_category_columns) > 0) {
        stop("VariableCategories data frames must contain Variable and Category.", call. = FALSE)
      }
      df_categories <- VariableCategories %>%
        dplyr::transmute(
          Variable = as.character(.data$Variable),
          Category = as.character(.data$Category)
        )
    } else if (!is.null(names(VariableCategories))) {
      df_categories <- tibble::tibble(
        Variable = names(VariableCategories),
        Category = as.character(VariableCategories)
      )
    } else {
      if (length(VariableCategories) != length(predictor_vars)) {
        stop("An unnamed VariableCategories vector must correspond to predictor_vars.", call. = FALSE)
      }
      df_categories <- tibble::tibble(
        Variable = predictor_vars,
        Category = as.character(VariableCategories)
      )
    }
    if (anyNA(df_categories$Variable) || any(!nzchar(df_categories$Variable))) {
      stop("VariableCategories contains a missing or empty variable name.", call. = FALSE)
    }
    conflicting_categories <- df_categories %>%
      dplyr::group_by(.data$Variable) %>%
      dplyr::summarise(
        NCategories = dplyr::n_distinct(.data$Category, na.rm = FALSE),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$NCategories > 1) %>%
      dplyr::pull(.data$Variable)
    if (length(conflicting_categories) > 0) {
      stop(
        "VariableCategories contains conflicting mappings for: ",
        paste(conflicting_categories, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    df_categories <- df_categories %>%
      dplyr::distinct(.data$Variable, .keep_all = TRUE)
    category_lookup <- stats::setNames(df_categories$Category, df_categories$Variable)
    category_order <- unique(df_categories$Category[!is.na(df_categories$Category)])
  }

  duplicate_predictors <- unique(predictor_vars[duplicated(predictor_vars)])
  if (length(duplicate_predictors) > 0) {
    warning(
      "Duplicated predictors were analyzed once and later occurrences were removed: ",
      paste(duplicate_predictors, collapse = ", "),
      ".",
      call. = FALSE
    )
    predictor_vars <- unique(predictor_vars)
  }
  missing_predictors <- setdiff(predictor_vars, names(data))
  if (length(missing_predictors) > 0) {
    warning(
      "Predictors not found in data were removed: ",
      paste(missing_predictors, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  predictor_vars <- predictor_vars[predictor_vars %in% names(data)]
  nonnumeric_predictors <- predictor_vars[!vapply(data[predictor_vars], is.numeric, logical(1))]
  if (length(nonnumeric_predictors) > 0) {
    warning(
      "Nonnumeric predictors were removed: ",
      paste(nonnumeric_predictors, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  predictor_vars <- setdiff(predictor_vars, nonnumeric_predictors)
  outcome_predictors <- intersect(predictor_vars, outcome_var)
  if (length(outcome_predictors) > 0) {
    warning("The outcome cannot also be a predictor and was removed.", call. = FALSE)
  }
  predictor_vars <- setdiff(predictor_vars, outcome_predictors)
  covariate_predictors <- intersect(predictor_vars, covariates)
  if (length(covariate_predictors) > 0) {
    warning(
      "Predictors also supplied as covariates were removed: ",
      paste(covariate_predictors, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  predictor_vars <- setdiff(predictor_vars, covariate_predictors)
  if (length(predictor_vars) == 0) {
    stop("No valid predictors remain to analyze.", call. = FALSE)
  }

  # Apply labels

  label_lookup <- stats::setNames(predictor_vars, predictor_vars)
  if (Relabel) {
    attribute_labels <- ScidrDisplayLabels(data, predictor_vars, Relabel = TRUE)
    label_lookup[names(attribute_labels)] <- attribute_labels
    if (!is.null(codebook)) {
      df_codebook_labels <- codebook %>%
        dplyr::transmute(
          Variable = as.character(.data$Variable),
          Label = as.character(.data$Label)
        ) %>%
        dplyr::filter(
          .data$Variable %in% predictor_vars,
          !is.na(.data$Label),
          nzchar(.data$Label)
        ) %>%
        dplyr::distinct(.data$Variable, .keep_all = TRUE)
      label_lookup[df_codebook_labels$Variable] <- df_codebook_labels$Label
    }
  }

  # Prepare models

  results <- CalculateContinuousScreeningEffects(
    data = data,
    predictor_vars = predictor_vars,
    outcome_var = outcome_var,
    covariates = covariates,
    label_lookup = label_lookup
  )
  invalid_results <- results %>%
    dplyr::filter(!is.na(.data$Note))
  if (nrow(invalid_results) > 0) {
    invalid_messages <- paste0(
      invalid_results$Variable,
      " (",
      invalid_results$Note,
      ")"
    )
    warning(
      "Predictors that could not be modeled were removed: ",
      paste(invalid_messages, collapse = "; "),
      ".",
      call. = FALSE
    )
  }
  results <- results %>%
    dplyr::filter(is.na(.data$Note))
  if (nrow(results) == 0) {
    stop("No predictors could be modeled successfully.", call. = FALSE)
  }

  results <- results %>%
    dplyr::mutate(
      Category = unname(category_lookup[.data$Variable]),
      FDR = stats::p.adjust(.data$PValue, method = AdjustMethod)
    )

  category_rank <- match(results$Category, category_order)
  if (length(category_order) == 0) {
    category_rank <- rep(1L, nrow(results))
  } else {
    category_rank[is.na(category_rank)] <- length(category_order) + 1L
  }
  results$.OriginalOrder <- match(results$Variable, predictor_vars)
  results$.CategoryOrder <- category_rank

  if (Sort == "original") {
    results <- results %>% dplyr::arrange(.data$.OriginalOrder)
  } else if (Sort == "pvalue") {
    results <- results %>% dplyr::arrange(.data$PValue, .data$.OriginalOrder)
  } else if (Sort == "fdr") {
    results <- results %>% dplyr::arrange(.data$FDR, .data$.OriginalOrder)
  } else if (Sort == "effect") {
    results <- results %>% dplyr::arrange(dplyr::desc(abs(.data$Beta)), .data$.OriginalOrder)
  } else if (Sort == "within_category_pvalue") {
    results <- results %>%
      dplyr::arrange(.data$.CategoryOrder, .data$PValue, .data$.OriginalOrder)
  } else {
    results <- results %>%
      dplyr::arrange(.data$.CategoryOrder, dplyr::desc(abs(.data$Beta)), .data$.OriginalOrder)
  }

  covariate_text <- if (length(covariates) == 0) "None" else paste(covariates, collapse = ", ")
  results <- results %>%
    dplyr::mutate(
      Tooltip = paste0(
        "<b>", .data$Label, "</b>",
        "<br>Variable: ", .data$Variable,
        "<br>Standardized beta: ", signif(.data$Beta, 3),
        "<br>95% CI: ", signif(.data$CILow, 3), " to ", signif(.data$CIHigh, 3),
        "<br>p-value: ", format.pval(.data$PValue, digits = 3, eps = 0.001),
        "<br>FDR: ", format.pval(.data$FDR, digits = 3, eps = 0.001),
        "<br>R: ", signif(.data$R, 3),
        ifelse(is.na(.data$AdjustedR), "", paste0("<br>Adjusted R: ", signif(.data$AdjustedR, 3))),
        "<br>N: ", .data$N,
        ifelse(is.na(.data$Category), "", paste0("<br>Category: ", .data$Category)),
        "<br>Covariates: ", covariate_text
      )
    )

  # Build plot

  plot_data <- results %>%
    dplyr::mutate(
      .AxisVariable = factor(.data$Variable, levels = .data$Variable),
      .PlotCategory = factor(.data$Category, levels = category_order)
    )
  axis_labels <- stats::setNames(results$Label, results$Variable)
  has_categories <- any(!is.na(results$Category))

  if (has_categories) {
    category_colors <- stats::setNames(
      .SciDataColorValues(length(category_order)),
      category_order
    )
    profile_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$.AxisVariable,
        y = .data$Beta,
        color = .data$.PlotCategory
      )
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$CILow, ymax = .data$CIHigh),
        width = 0.18,
        alpha = 0.5,
        linewidth = 0.55
      ) +
      ggplot2::scale_color_manual(
        values = category_colors,
        breaks = category_order,
        na.value = "grey70",
        na.translate = FALSE,
        name = "Category"
      )
  } else {
    profile_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$.AxisVariable, y = .data$Beta)
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$CILow, ymax = .data$CIHigh),
        width = 0.18,
        alpha = 0.5,
        linewidth = 0.55,
        color = .SciDataColorValues(1)[1]
      )
  }

  point_layer <- if (InteractiveLabels) {
    ggplot2::layer(
      geom = "point",
      stat = "identity",
      position = "identity",
      mapping = ggplot2::aes(text = .data$Tooltip),
      params = if (has_categories) {
        list(size = 2.4)
      } else {
        list(size = 2.4, color = .SciDataColorValues(1)[1])
      },
      inherit.aes = TRUE,
      check.aes = FALSE,
      show.legend = NA
    )
  } else if (has_categories) {
    ggplot2::geom_point(size = 2.4)
  } else {
    ggplot2::geom_point(size = 2.4, color = .SciDataColorValues(1)[1])
  }

  profile_plot <- profile_plot +
    point_layer +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey65",
      linewidth = 0.4
    ) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::labs(x = NULL, y = "Standardized Beta") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (has_categories) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (RemoveXAxisLabels) {
    profile_plot <- profile_plot +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else {
    profile_plot <- profile_plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
      )
  }

  # Return result

  results <- results %>%
    dplyr::select(dplyr::all_of(c(
      "Variable", "Label", "Category", "Beta", "SE", "CILow", "CIHigh",
      "PValue", "FDR", "N", "R", "AdjustedR", "Tooltip"
    )))

  list(
    Plot = profile_plot,
    ResultsTable = results,
    Metadata = list(
      Outcome = outcome_var,
      PredictorsRequested = predictors_requested,
      PredictorsAnalyzed = results$Variable,
      Covariates = covariates,
      AdjustMethod = AdjustMethod,
      Alpha = Alpha,
      Sort = Sort,
      ConfidenceLevel = 0.95,
      NModels = nrow(results)
    )
  )
}
