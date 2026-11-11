#' Plot Single Interaction Effect
#'
#' Creates a scatter plot with regression lines showing the interaction between a predictor
#' and outcome variable, moderated by either a continuous or categorical variable.
#'
#' @param data A data frame containing the variables to be analyzed
#' @param interVar Character string specifying the interaction variable (moderator)
#' @param outcome_var Character string specifying the outcome variable
#' @param predictor_var Character string specifying the predictor variable
#' @param covariates Character vector of covariate names to include in the model
#' @param n_lines For continuous moderators, number of lines to plot (default: 3 for low/med/high)
#' @param alpha Transparency level for points (default: 0.6)
#' @param point_size Size of points (default: 2)
#'
#' @return A ggplot object showing the interaction effect
#'
#' @details
#' For categorical moderators, separate regression lines are plotted for each category.
#' For continuous moderators, regression lines are plotted at the mean and +/- 1 SD.
#'
#' Model-based partial residuals preserve the fitted interaction slopes.
#' Ribbons show 95% confidence intervals for the displayed model component.
#'
#' Variable labels are used if available in the data frame.
#'
#' @param adjust Display model-based covariate-adjusted partial residuals.
#' @param standardize Use pooled model-sample SD units for both axes.
#' @param covariate_interactions Allow moderator-specific covariate effects.
#' @param TreatOrdinalAs Ordinal treatment passed to PlotInteractionEffectsMatrix.
#' @param source_formula Optional stored matrix formula; validated against the
#'   supplied model specification and used for adjusted displays.
#' @param interaction_results Optional source interaction row, retaining its
#'   original multiple-testing q-value in the annotation.
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param outcomeVar \strong{Deprecated} (since 19.15.0). Use \code{outcome_var} instead.
#' @param predictorVar \strong{Deprecated} (since 19.15.0). Use \code{predictor_var} instead.
#' @param covars \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels for readable axis titles and legend
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Categorical moderator: one regression line per Diagnosis group
#' PlotInteractionEffectsContinuous(
#'   Labelled,
#'   interVar = "Diagnosis",
#'   outcome_var = "AXL",
#'   predictor_var = "age"
#' )
#'
#' # Continuous moderator: lines at mean and +/- 1 SD of the moderator
#' PlotInteractionEffectsContinuous(
#'   Labelled,
#'   interVar = "Adiponectin",
#'   outcome_var = "AXL",
#'   predictor_var = "age"
#' )
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank
#' @importFrom dplyr mutate case_when
#' @importFrom stats lm coef sd quantile
#' @importFrom paletteer paletteer_d paletteer_c

#' @inheritParams inference-results
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
PlotInteractionEffectsContinuous <- function(data,
    interVar = NULL,
    outcome_var = NULL,
    predictor_var = NULL,
    covariates = NULL,
    n_lines = 3,
    alpha = 0.6,
    point_size = 2,
    Data = lifecycle::deprecated(),
    outcomeVar = lifecycle::deprecated(),
    predictorVar = lifecycle::deprecated(),
    covars = lifecycle::deprecated(),
    Pthresh = 0.05,
    PVar = "PValue",
    adjust = TRUE, standardize = FALSE, covariate_interactions = FALSE,
    TreatOrdinalAs = "Exclude", interaction_results = NULL, source_formula = NULL) {
  .scidr_inference <- ScidrStartInference("PlotInteractionEffectsContinuous", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsContinuous(Data)", "PlotInteractionEffectsContinuous(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(outcomeVar)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsContinuous(outcomeVar)", "PlotInteractionEffectsContinuous(outcome_var)")
    outcome_var <- outcomeVar
  }
  outcomeVar <- outcome_var
  if (lifecycle::is_present(predictorVar)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsContinuous(predictorVar)", "PlotInteractionEffectsContinuous(predictor_var)")
    predictor_var <- predictorVar
  }
  predictorVar <- predictor_var
  if (lifecycle::is_present(covars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsContinuous(covars)", "PlotInteractionEffectsContinuous(covariates)")
    covariates <- covars
  }
  covars <- covariates


  # Input validation
  if (is.null(Data) || !is.data.frame(Data)) {
    stop("Data must be a non-null data frame")
  }

  if (is.null(interVar) || is.null(outcomeVar) || is.null(predictorVar)) {
    stop("interVar, outcomeVar, and predictorVar must all be specified")
  }

  if (!all(c(interVar, outcomeVar, predictorVar) %in% names(Data))) {
    stop("All specified variables must exist in the data frame")
  }

  if (!is.null(covars) && !all(covars %in% names(Data))) {
    stop("All covariates must exist in the data frame")
  }

  if (!is.logical(adjust) || length(adjust) != 1L || is.na(adjust) ||
      !is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) stop("adjust and standardize must be TRUE or FALSE.")
  if (length(n_lines) != 1L || !is.finite(n_lines) || n_lines < 2 || n_lines != as.integer(n_lines)) stop("n_lines must be an integer >= 2.")
  # Reuse the matrix inference engine, including pooled scaling, exact term
  # identification, joint tests, and complete-case selection.
  inference <- PlotInteractionEffectsMatrix(Data, interVar, outcomeVar, predictorVar,
    covariates = if (adjust) covars else NULL,
    covariate_interactions = covariate_interactions, TreatOrdinalAs = TreatOrdinalAs)
  if (any(inference$Results$Status == "failed")) {
    p <- ggplot2::ggplot() + ggplot2::theme_void()
    p$Results <- inference$Results
    p$Metadata <- inference$Metadata
    return(p)
  }
  prep <- ConvertOrdinalToNumeric(Data, unique(c(outcomeVar, predictorVar, interVar, covars)),
    TreatOrdinalAs = TreatOrdinalAs, ReturnMetadata = TRUE)
  df_Input <- prep$data
  formula_str <- inference$Metadata$Models$Formula[1]
  if (adjust && !is.null(source_formula)) {
    inferred_terms <- stats::terms(stats::as.formula(formula_str))
    source_terms <- stats::terms(stats::as.formula(source_formula))
    if (!setequal(attr(inferred_terms, "term.labels"), attr(source_terms, "term.labels")) ||
        !identical(all.vars(stats::formula(inferred_terms))[[1]], all.vars(stats::formula(source_terms))[[1]])) {
      stop("Source formula differs from the supported matrix specification; regenerate InteractionObject.")
    }
    formula_str <- source_formula
  }
  model <- ScidrObserveInference(stats::lm(stats::as.formula(formula_str), data = df_Input))
  df_Model <- stats::model.frame(model)
  continuous <- is.numeric(df_Model[[interVar]])
  labels <- ScidrDisplayLabels(Data, c(outcomeVar, predictorVar, interVar))
  x_sd <- stats::sd(df_Model[[predictorVar]])
  y_sd <- stats::sd(df_Model[[outcomeVar]])
  x_mean <- mean(df_Model[[predictorVar]])
  y_mean <- mean(df_Model[[outcomeVar]])
  # Model-matrix derivatives avoid dependence on factor contrast coding.
  Design <- function(df_New) stats::model.matrix(stats::delete.response(stats::terms(model)),
    df_New, contrasts.arg = model$contrasts, xlev = model$xlevels)
  Derivative <- function(df_New) {
    df_Zero <- df_New
    df_One <- df_New
    df_Zero[[predictorVar]] <- 0
    df_One[[predictorVar]] <- 1
    Design(df_One) - Design(df_Zero)
  }
  if (continuous) {
    mod_mean <- mean(df_Model[[interVar]])
    mod_sd <- stats::sd(df_Model[[interVar]])
    mod_values <- seq(mod_mean - mod_sd, mod_mean + mod_sd, length.out = n_lines)
    group_labels <- if (n_lines == 2L) c("Mean - 1 SD", "Mean + 1 SD") else paste0("Moderator = ", signif(mod_values, 3))
  } else {
    mod_values <- levels(factor(df_Model[[interVar]]))
    group_labels <- mod_values
  }
  df_Points <- data.frame(predictor = df_Model[[predictorVar]], outcome = df_Model[[outcomeVar]],
    moderator = df_Model[[interVar]])
  if (adjust) {
    nuisance <- setdiff(covars, c(interVar, outcomeVar, predictorVar))
    rhs <- c(ScidrQuoteFormulaNames(interVar), ScidrQuoteFormulaNames(nuisance))
    if (covariate_interactions && length(nuisance)) rhs <- c(rhs, paste(ScidrQuoteFormulaNames(nuisance), "*", ScidrQuoteFormulaNames(interVar)))
    nuisance_formula <- stats::as.formula(paste(ScidrQuoteFormulaNames(predictorVar), "~", paste(rhs, collapse = " + ")))
    predictor_residual <- stats::residuals(stats::lm(nuisance_formula, data = df_Model))
    slopes <- drop(Derivative(df_Model) %*% stats::coef(model))
    df_Points$predictor <- predictor_residual
    # Partial residuals retain the original model's predictor contribution,
    # after removing nuisance effects and the removed predictor component.
    df_Points$outcome <- stats::residuals(model) + slopes * predictor_residual
  }
  if (standardize) {
    df_Points$predictor <- (df_Points$predictor - if (adjust) 0 else x_mean) / x_sd
    df_Points$outcome <- (df_Points$outcome - if (adjust) 0 else y_mean) / y_sd
  }
  grid_x <- seq(min(df_Points$predictor), max(df_Points$predictor), length.out = 80)
  curves <- lapply(seq_along(mod_values), function(i) {
    df_New <- df_Model[rep(1L, length(grid_x)), , drop = FALSE]
    df_New[[interVar]] <- if (continuous) rep(mod_values[i], length(grid_x)) else factor(rep(mod_values[i], length(grid_x)), levels = model$xlevels[[interVar]])
    if (adjust) {
      weights <- Derivative(df_New) * (grid_x * if (standardize) x_sd else 1)
    } else {
      df_New[[predictorVar]] <- if (standardize) grid_x * x_sd + x_mean else grid_x
      weights <- Design(df_New)
    }
    fit <- drop(weights %*% stats::coef(model))
    se <- sqrt(pmax(0, rowSums((weights %*% stats::vcov(model)) * weights)))
    if (standardize) {
      fit <- (fit - if (adjust) 0 else y_mean) / y_sd
      se <- se / y_sd
    }
    critical <- stats::qt(0.975, stats::df.residual(model))
    data.frame(predictor = grid_x, fit = fit, lower = fit - critical * se,
      upper = fit + critical * se, Group = group_labels[i])
  })
  df_Curves <- dplyr::bind_rows(curves)
  palette <- stats::setNames(.SciDataColorValues(length(group_labels)), group_labels)
  p <- ggplot2::ggplot(df_Points, ggplot2::aes(x = predictor, y = outcome))
  if (continuous) {
    # A continuous point colour scale is separate from categorical line colours.
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = moderator), alpha = alpha, size = point_size) +
      ggplot2::scale_colour_gradient(low = palette[1], high = palette[length(palette)], name = labels[[interVar]])
    line_colors <- scales::col_numeric(c(palette[1], palette[length(palette)]),
      domain = range(df_Points$moderator))(mod_values)
    for (i in seq_along(group_labels)) {
      df_Curve <- df_Curves[df_Curves$Group == group_labels[i], ]
      p <- p + ggplot2::geom_ribbon(data = df_Curve, ggplot2::aes(x = predictor, ymin = lower, ymax = upper),
        inherit.aes = FALSE, fill = line_colors[i], alpha = 0.15) +
        ggplot2::geom_line(data = df_Curve, ggplot2::aes(x = predictor, y = fit), inherit.aes = FALSE, colour = line_colors[i])
    }
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = moderator), alpha = alpha, size = point_size) +
      ggplot2::geom_ribbon(data = df_Curves, ggplot2::aes(x = predictor, ymin = lower, ymax = upper, fill = Group, group = Group), inherit.aes = FALSE, alpha = 0.15) +
      ggplot2::geom_line(data = df_Curves, ggplot2::aes(x = predictor, y = fit, colour = Group, group = Group), inherit.aes = FALSE) +
      ggplot2::scale_colour_manual(values = palette, name = labels[[interVar]]) +
      ggplot2::scale_fill_manual(values = palette, name = labels[[interVar]])
  }
  df_Slopes <- inference$Results[inference$Results$Term != "Interaction", ]
  columns <- if (standardize) c("EstimateStd", "ConfLowStd", "ConfHighStd") else c("Estimate", "ConfLow", "ConfHigh")
  units <- if (standardize) "SD/SD" else "outcome units per predictor unit"
  slope_notes <- paste0(df_Slopes$Level, ": slope (", units, ") = ", signif(df_Slopes[[columns[1]]], 3),
    " [", signif(df_Slopes[[columns[2]]], 3), ", ", signif(df_Slopes[[columns[3]]], 3), "]; n = ",
    if (continuous) df_Slopes$N else df_Slopes$NGroup)
  df_Interaction <- if (is.null(interaction_results)) inference$Results[inference$Results$Term == "Interaction", ] else interaction_results
  direction <- if (all(c("ComparisonLevel", "ReferenceLevel") %in% names(df_Interaction)) && !is.na(df_Interaction$ComparisonLevel[1])) paste0(" (", df_Interaction$ComparisonLevel[1], " - ", df_Interaction$ReferenceLevel[1], ")") else ""
  note <- paste(c(slope_notes, paste0("Interaction", direction, ", ", units, " = ", signif(df_Interaction[[columns[1]]], 3),
    " [", signif(df_Interaction[[columns[2]]], 3), ", ", signif(df_Interaction[[columns[3]]], 3), "]; p = ",
    format.pval(df_Interaction$PValue, digits = 3), "; q = ", format.pval(df_Interaction$PValueFDR, digits = 3))), collapse = "\n")
  p <- p + ggplot2::labs(caption = stringr::str_wrap(note, width = 85),
      x = paste0(labels[[predictorVar]], if (standardize) " (SD units)" else ""),
      y = paste0(labels[[outcomeVar]], if (standardize) " (SD units)" else "")) +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title.position = "plot",
      plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0))
  p$Results <- inference$Results
  p$Metadata <- inference$Metadata
  p$Metadata$Display <- list(adjust = adjust, standardize = standardize,
    PartialResidualDefinition = "x_resid = x - fitted nuisance model; y_partial = model residual + slope(moderator) * x_resid",
    SD = c(Predictor = x_sd, Outcome = y_sd))
  attr(p, "model") <- model
  attr(p, "interaction_p") <- df_Interaction$PValue
  attr(p, "formula") <- formula_str
  attr(p, "LineData") <- df_Curves
  return(p)
}
