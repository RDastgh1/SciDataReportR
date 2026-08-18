#' PlotMiningMatrix
#'
#' Generate a matrix of statistical relationships between variables.
#'
#' @param data A data frame.
#' @param outcome_vars Outcome variables.
#' @param predictor_vars Predictor variables. If NULL, uses OutcomeVars.
#' @param covariates Optional covariates (reserved for future use).
#' @param Relabel Use labels instead of names.
#' @param TreatOrdinalAs How ordinal variables are handled: `"Categorical"`,
#' `"Continuous"`, `"Both"`, or `"Exclude"`.
#' @param Parametric Use parametric tests.
#'
#' @return List with tables and plots.
#' @param fdr_scope Either `"matrix"` (default) or `"per_outcome"`, passed to
#'   [ApplyFDRCorrection()]. `"matrix"` corrects across all pairwise p-values
#'   at once (historical behavior, computed on the symmetrized pair table).
#'   `"per_outcome"` corrects separately within each x-axis variable (`XVar`,
#'   ordered by `outcome_vars`).
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param OutcomeVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @param PredictorVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param Covariates \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels for readable axes
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # A mining matrix over 11 mixed-type variables (categorical + continuous)
#' result <- PlotMiningMatrix(
#'   Labelled,
#'   outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
#'   predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
#'                      "Apolipoprotein_B", "C_Reactive_Protein",
#'                      "Cortisol", "Insulin")
#' )
#'
#' # The relationship plot (point shape/size encode significance from raw p)
#' result$Unadjusted$plot
#'
#' # The p-value table carries both unadjusted (p) and FDR-adjusted (p_adj)
#' # p-values so results can be inspected with and without FDR correction
#' result$Unadjusted$PvalTable[, c("XVar", "YVar", "p", "p_adj", "Test")]
#'
#' # Per-outcome FDR correction instead of matrix-wide
#' result_perout <- PlotMiningMatrix(
#'   Labelled,
#'   outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
#'   predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
#'                      "Apolipoprotein_B", "C_Reactive_Protein",
#'                      "Cortisol", "Insulin"),
#'   fdr_scope = "per_outcome"
#' )
#'
#' # An interactive version of the matrix
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   plotly::ggplotly(result$Unadjusted$plot)
#' }
#' @export
PlotMiningMatrix <- function(data,
    outcome_vars,
    predictor_vars = NULL,
    covariates = NULL,
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    Parametric = TRUE,
    Data = lifecycle::deprecated(),
    OutcomeVars = lifecycle::deprecated(),
    PredictorVars = lifecycle::deprecated(),
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    Covariates = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMiningMatrix(Data)", "PlotMiningMatrix(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(OutcomeVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMiningMatrix(OutcomeVars)", "PlotMiningMatrix(outcome_vars)")
    outcome_vars <- OutcomeVars
  }
  if (!missing(outcome_vars)) OutcomeVars <- outcome_vars
  if (lifecycle::is_present(PredictorVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMiningMatrix(PredictorVars)", "PlotMiningMatrix(predictor_vars)")
    predictor_vars <- PredictorVars
  }
  PredictorVars <- predictor_vars
  if (lifecycle::is_present(Covariates)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMiningMatrix(Covariates)", "PlotMiningMatrix(covariates)")
    covariates <- Covariates
  }
  Covariates <- covariates
  fdr_scope <- match.arg(fdr_scope)


  method <- ifelse(Parametric, "pearson", "spearman")

  # Supplied names must exist. Intersecting them away silently shrank the matrix
  # and, for covariates, produced unadjusted tests reported as adjusted.
  OutcomeVars <- ScidrValidateVariables(Data, OutcomeVars, "outcome_vars")
  ScidrValidateVariables(Data, Covariates, "covariates")

  if (is.null(PredictorVars)) {
    PredictorVars <- OutcomeVars
  } else {
    PredictorVars <- ScidrValidateVariables(Data, PredictorVars, "predictor_vars")
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) {
    empty <- data.frame()
    p0 <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = p0)
    ))
  }

  ordinal <- ConvertOrdinalToNumeric(
    Data, unique(c(OutcomeVars, PredictorVars)), TreatOrdinalAs = TreatOrdinalAs,
    Relabel = Relabel, ReturnMetadata = TRUE
  )
  Data <- ordinal$data
  OutcomeVars <- unique(unlist(ordinal$variable_map[OutcomeVars], use.names = FALSE))
  PredictorVars <- unique(unlist(ordinal$variable_map[PredictorVars], use.names = FALSE))

  coalesce_p <- function(df) {
    out <- rep(NA_real_, nrow(df))
    for (nm in c("p", "P", "pval", "p.value", "p_value")) {
      if (nm %in% names(df)) {
        out <- dplyr::coalesce(out, suppressWarnings(as.numeric(df[[nm]])))
      }
    }
    out
  }

  coalesce_effect <- function(df) {
    out <- rep(NA_real_, nrow(df))
    for (nm in c("r", "R", "estimate", "cor", "correlation", "ges", "cramers_v")) {
      if (nm %in% names(df)) {
        out <- dplyr::coalesce(out, suppressWarnings(as.numeric(df[[nm]])))
      }
    }
    out
  }

  labels <- ScidrDisplayLabels(Data, names(Data), Relabel)

  # SAFE LOOKUP FUNCTION (never returns NA)
  safe_lookup <- function(vars, labels) {
    out <- labels[vars]
    missing_idx <- is.na(out) | out == ""
    out[missing_idx] <- vars[missing_idx]
    return(out)
  }

  num_vars <- SciDataReportR::getNumVars(Data)
  cat_vars <- SciDataReportR::getCatVars(Data)

  num_Outcomes   <- intersect(OutcomeVars, num_vars)
  cat_Outcomes   <- intersect(OutcomeVars, cat_vars)
  num_Predictors <- intersect(PredictorVars, num_vars)
  cat_Predictors <- intersect(PredictorVars, cat_vars)

  if (!is.null(Covariates)) {
    num_Outcomes   <- setdiff(num_Outcomes, Covariates)
    num_Predictors <- setdiff(num_Predictors, Covariates)
    cat_Outcomes   <- setdiff(cat_Outcomes, Covariates)
    cat_Predictors <- setdiff(cat_Predictors, Covariates)
  }

  results_list <- list()

  # Correlations
  if (length(num_Outcomes) > 0 && length(num_Predictors) > 0) {

    cor_res <- tryCatch(
      SciDataReportR::PlotCorrelationsHeatmap(
        Data,
        predictor_vars = num_Predictors,
        outcome_vars = num_Outcomes,
        covariates = Covariates,
        method = method,
        Relabel = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(cor_res) && !is.null(cor_res$Unadjusted$plot$data)) {

      df_cor <- cor_res$Unadjusted$plot$data %>%
        dplyr::mutate(
          p = coalesce_p(.),
          EffectSize = coalesce_effect(.),
          Test = "Correlation"
        ) %>%
        dplyr::select(XVar, YVar, p, EffectSize, Test) %>%
        dplyr::filter(XVar != YVar)

      results_list[[length(results_list) + 1]] <- df_cor
    }
  }

  # ANOVA
  if (length(cat_Predictors) > 0 && length(num_Outcomes) > 0) {

    anova_res <- tryCatch(
      SciDataReportR::PlotAnovaRelationshipsMatrix(
        Data,
        CatVars = cat_Predictors,
        ContVars = num_Outcomes,
        covariates = Covariates,
        Relabel = FALSE,
        Parametric = Parametric
      ),
      error = function(e) NULL
    )

    if (!is.null(anova_res) && !is.null(anova_res$Unadjusted$PvalTable)) {

      df_anova <- anova_res$Unadjusted$PvalTable %>%
        dplyr::mutate(
          XVar = ContinuousVariable,
          YVar = CategoricalVariable,
          p = coalesce_p(.),
          EffectSize = coalesce_effect(.),
          Test = "ANOVA"
        ) %>%
        dplyr::select(XVar, YVar, p, EffectSize, Test) %>%
        dplyr::filter(XVar != YVar)

      results_list[[length(results_list) + 1]] <- df_anova
    }
  }

  # ChiSq
  if (length(cat_Predictors) > 0 && length(cat_Outcomes) > 0) {

    chi_res <- tryCatch(
      SciDataReportR::PlotChiSqCovar(
        Data,
        predictor_vars = cat_Predictors,
        outcome_vars = cat_Outcomes,
        covariates = Covariates,
        Relabel = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(chi_res) && !is.null(chi_res$p$data)) {

      df_chi <- chi_res$p$data %>%
        dplyr::mutate(
          p = coalesce_p(.),
          EffectSize = coalesce_effect(.),
          Test = "ChiSq"
        ) %>%
        dplyr::select(XVar, YVar, p, EffectSize, Test) %>%
        dplyr::filter(XVar != YVar)

      results_list[[length(results_list) + 1]] <- df_chi
    }
  }

  if (length(results_list) == 0) {
    return(list(
      Unadjusted = list(PvalTable = data.frame(), plot = ggplot2::ggplot() + ggplot2::theme_void())
    ))
  }

  # Symmetry fix
  results <- dplyr::bind_rows(results_list) %>%
    dplyr::filter(
      !is.na(p),
      !is.na(EffectSize),
      is.finite(EffectSize)
    ) %>%
    dplyr::mutate(
      VarA = pmin(XVar, YVar),
      VarB = pmax(XVar, YVar)
    ) %>%
    dplyr::group_by(VarA, VarB) %>%
    dplyr::summarise(
      p = min(p, na.rm = TRUE),
      EffectSize = max(abs(EffectSize), na.rm = TRUE),
      Test = dplyr::first(Test),
      .groups = "drop"
    )

  results <- dplyr::bind_rows(
    results %>% dplyr::transmute(XVar = VarA, YVar = VarB, p, EffectSize, Test),
    results %>% dplyr::transmute(XVar = VarB, YVar = VarA, p, EffectSize, Test)
  )

  # The pair table is symmetrized (each pair appears in both orientations)
  # before correction; "matrix" scope reproduces the historical behavior on
  # that duplicated table. For "per_outcome", p-values are grouped by XVar,
  # the x-axis of the final plot, which is ordered by outcome_vars.
  results$p_adj <- ApplyFDRCorrection(
    results$p,
    fdr_scope = fdr_scope,
    outcome_ids = results$XVar,
    predictor_ids = results$YVar
  )

  results <- results %>%
    rstatix::add_significance(p.col = "p", output.col = "stars")

  # Safe label application
  if (Relabel) {
    results$XLabel <- safe_lookup(results$XVar, labels)
    results$YLabel <- safe_lookup(results$YVar, labels)
  } else {
    results$XLabel <- results$XVar
    results$YLabel <- results$YVar
  }

  # Safe ordering with no NA levels
  x_order <- if (Relabel) safe_lookup(OutcomeVars, labels) else OutcomeVars
  y_order <- if (Relabel) safe_lookup(PredictorVars, labels) else PredictorVars

  results$XLabel <- factor(results$XLabel, levels = unique(x_order))
  results$YLabel <- factor(results$YLabel, levels = rev(unique(y_order)))

  # FINAL CLEANUP
  results <- results %>%
    dplyr::filter(!is.na(XLabel), !is.na(YLabel))

  results$EffectSizeAbs <- abs(results$EffectSize)

  # rstatix also emits "****"; fold it into "***" so every band has a shape and
  # the legend matches the four documented cutpoints.
  results$stars <- as.character(results$stars)
  results$stars[results$stars == "****"] <- "***"
  results$stars <- factor(results$stars, levels = c("ns", "*", "**", "***"))

  # Sizes are compensated for the differing apparent area of shapes 15-18 so
  # that the marker grows monotonically with significance.
  size_map  <- c("ns" = 2, "*" = 3.5, "**" = 4, "***" = 6)
  shape_map <- c("ns" = 16, "*" = 17, "**" = 15, "***" = 18)

  results$size_val <- size_map[as.character(results$stars)]

  shape_labels <- c(
    "ns" = "ns (p >= 0.05)",
    "*"  = "* (p <= 0.05)",
    "**" = "** (p <= 0.01)",
    "***"= "*** (p <= 0.001)"
  )

  p <- ggplot2::ggplot(results, ggplot2::aes(x = XLabel, y = YLabel)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = EffectSizeAbs, size = size_val, shape = stars)
    ) +
    ggplot2::scale_colour_gradient(
      low = "grey85", high = "#841B37", limits = c(0, NA),
      name = "Effect size", na.value = "grey70") +
    ggplot2::scale_shape_manual(
      values = shape_map, labels = shape_labels, name = "Significance",
      drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank()
    )

  out <- list(
    Unadjusted = list(PvalTable = results, plot = p),
    Metadata = list(
      TreatOrdinalAs = TreatOrdinalAs,
      DisplayLabels = labels[unique(c(OutcomeVars, PredictorVars))]
    )
  )
  # Standardized p-value element alias (old name kept)
  out$p <- out$Unadjusted
  out
}
