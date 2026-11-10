#' Plot significant interaction pairs
#'
#' Select interaction rows from a matrix result and reuse
#' [PlotInteractionEffectsContinuous()] for each single-pair display.
#' @param data Data used to fit the source matrix.
#' @param InteractionObject Output of [PlotInteractionEffectsMatrix()] with
#'   Results and Metadata. Plot layer data is never used for selection.
#' @param PVar Selection column: "PValue" or "PValueFDR".
#' @param Pthresh Select p-values strictly below this cutoff.
#' @param adjust Display model-based covariate-adjusted partial residuals.
#'   The predictor is residualized on nuisance covariates and moderator terms;
#'   the outcome partial residual is the full-model residual plus the fitted
#'   predictor contribution on that residualized predictor. Independent
#'   covariate-only residualization of both axes would not preserve interaction
#'   slopes in general. FALSE displays raw values and an unadjusted model.
#' @param standardize Scale axes using full fitted-sample SDs, pooled over groups.
#' @param ncol Optional positive integer requesting a patchwork layout.
#' @param max_plots Maximum number of plots, ordered by selected p-value.
#' @return Named list of ggplots, with selected interaction rows in attribute
#'   Results. When ncol is supplied, attribute Layout holds a patchwork.
#'   Empty selections return an empty list with a message.
#' @details Annotations report displayed-model standardized simple slopes,
#'   intervals and sample sizes. Interaction estimate, p and q retain the
#'   original source inference even when adjust is FALSE. Metadata and source
#'   formulas are retained as attributes. Continuous moderator lines are at
#'   the model-sample mean minus/plus one SD; point colours remain continuous.
#' @export
plotSigInteractions <- function(data, InteractionObject, PVar = "PValue",
    Pthresh = 0.05, adjust = TRUE, standardize = TRUE, ncol = NULL,
    max_plots = 25) {
  if (!is.list(InteractionObject)) stop("InteractionObject must contain Results and Metadata from PlotInteractionEffectsMatrix().")
  df_Results <- InteractionObject$Results
  metadata <- InteractionObject$Metadata
  required <- c("Outcome", "Predictor", "Term", "Estimate", "PValue", "PValueFDR", "OutcomeLabel", "PredictorLabel")
  if (!is.data.frame(df_Results) || !all(required %in% names(df_Results)) ||
      !is.list(metadata) || !is.list(metadata$Arguments) ||
      !is.data.frame(metadata$Models) || !all(c("Outcome", "Predictor", "Formula") %in% names(metadata$Models))) {
    stop("InteractionObject must contain Results and Metadata from PlotInteractionEffectsMatrix(); regenerate the matrix with the current package.", call. = FALSE)
  }
  spec <- metadata$Arguments
  if (!all(c("interVar", "covariates", "covariate_interactions", "TreatOrdinalAs") %in% names(spec)) || is.null(spec$interVar) || is.null(spec$covariate_interactions) || is.null(spec$TreatOrdinalAs)) {
    stop("InteractionObject Metadata is missing the moderator or model specification.", call. = FALSE)
  }
  if (!PVar %in% c("PValue", "PValueFDR")) stop("PVar must be PValue or PValueFDR.")
  if (length(Pthresh) != 1L || !is.finite(Pthresh) || Pthresh <= 0 || Pthresh >= 1) stop("Pthresh must lie between 0 and 1.")
  if (length(max_plots) != 1L || !is.finite(max_plots) || max_plots < 1 || max_plots != floor(max_plots)) stop("max_plots must be a positive integer.")
  if (!is.null(ncol) && (length(ncol) != 1L || !is.finite(ncol) || ncol < 1 || ncol != floor(ncol))) stop("ncol must be a positive integer.")
  df_Selected <- df_Results[which(df_Results$Term == "Interaction" &
    is.finite(df_Results[[PVar]]) & df_Results[[PVar]] < Pthresh), , drop = FALSE]
  df_Selected <- df_Selected[order(df_Selected[[PVar]]), , drop = FALSE]
  if (nrow(df_Selected) > max_plots) {
    warning("More than max_plots interactions passed; keeping the ", max_plots, " smallest ", PVar, " values.", call. = FALSE)
    df_Selected <- utils::head(df_Selected, max_plots)
  }
  plots <- stats::setNames(list(), character())
  if (!nrow(df_Selected)) {
    message("No interactions with ", PVar, " < ", Pthresh, ".")
    attr(plots, "Results") <- df_Selected
    attr(plots, "Metadata") <- metadata
    return(plots)
  }
  if (!is.null(ncol) && !requireNamespace("patchwork", quietly = TRUE)) stop("Install 'patchwork' to request a combined layout.")
  for (i in seq_len(nrow(df_Selected))) {
    row <- df_Selected[i, ]
    models <- metadata$Models
    formula <- models$Formula[models$Outcome == row$Outcome & models$Predictor == row$Predictor]
    if (!length(formula) || anyNA(formula)) stop("Source model formula is unavailable for ", row$Outcome, " / ", row$Predictor, ".")
    plot <- PlotInteractionEffectsContinuous(data, interVar = spec$interVar,
      outcome_var = row$Outcome, predictor_var = row$Predictor,
      covariates = spec$covariates, covariate_interactions = spec$covariate_interactions,
      TreatOrdinalAs = spec$TreatOrdinalAs, adjust = adjust, standardize = standardize,
      n_lines = 2, interaction_results = row, source_formula = formula[1])
    plot <- plot + ggplot2::labs(
      x = paste0(row$PredictorLabel, if (standardize) " (SD units)" else ""),
      y = paste0(row$OutcomeLabel, if (standardize) " (SD units)" else ""))
    attr(plot, "SourceFormula") <- formula
    plots[[paste(row$Outcome, row$Predictor, sep = "__")]] <- plot
  }
  attr(plots, "Results") <- df_Selected
  attr(plots, "Metadata") <- metadata
  if (!is.null(ncol)) attr(plots, "Layout") <- patchwork::wrap_plots(plots, ncol = ncol)
  plots
}
