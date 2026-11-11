#' PlotMiningMatrix
#'
#' Generate a matrix of statistical relationships between variables.
#'
#' @param data A data frame.
#' @param outcome_vars Outcome variables.
#' @param predictor_vars Predictor variables. If NULL, uses OutcomeVars.
#' @param covariates Optional variables to exclude from the tested axes.
#' @param Relabel Use labels instead of names.
#' @param TreatOrdinalAs How ordinal variables are handled: `"Categorical"`,
#' `"Continuous"`, `"Both"`, or `"Exclude"`.
#' @param Parametric Use parametric tests.
#'
#' @return List with tables and plots.
#' @details Numeric pairs use Pearson correlation when `Parametric = TRUE`,
#' and Spearman otherwise. Categorical/numeric pairs use ANOVA or Kruskal-Wallis
#' in either orientation; categorical pairs use chi-squared tests. Logical
#' columns are converted locally to FALSE/TRUE factors, preserving labels.
#' Missing effect sizes retain their p-values and plot in neutral grey.
#' Colors represent absolute effect sizes on a fixed 0--1 scale with a
#' square-root transformation to emphasize smaller effects. Values outside
#' these limits are squished to the endpoints for color mapping only;
#' legend labels and returned estimates remain in original effect-size units.
#' Significance uses the p-values selected by `PVar`: circles for p > 0.05,
#' triangles for 0.01 < p <= 0.05, squares for 0.001 < p <= 0.01, and
#' diamonds for p <= 0.001. With `IncludeTrendSymbol = TRUE`, plus signs
#' represent 0.05 < p <= 0.10 and circles represent p > 0.10.
#' Marker sizes are fixed across plots: circle 2, plus 2.7, triangle 3.5,
#' square 4, and diamond 6. The legend includes all enabled bands at these sizes.
#' `Unadjusted$Missing` lists requested non-self pairs without a finite p-value.
#' @param fdr_scope One of `"matrix"`, `"per_outcome"`, or `"per_predictor"`.
#' Correction uses the symmetrized table before selecting the displayed grid.
#' `"matrix"` corrects all p-values together (historical behavior).
#' `"per_outcome"` groups by x-axis outcome (`XVar`), correcting over predictors;
#' `"per_predictor"` groups by y-axis predictor (`YVar`), correcting over outcomes.
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
#' @param IncludeTrendSymbol Logical; add a plus sign for `0.05 < p <= 0.10`.
#' Defaults to FALSE; the trend legend entry is omitted when disabled.
#' @export
#' @inheritParams inference-results
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
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
    IncludeTrendSymbol = FALSE,
    Covariates = lifecycle::deprecated(),
    Pthresh = 0.05,
    PVar = "PValue") {
  .scidr_inference <- ScidrStartInference("PlotMiningMatrix", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

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
      Unadjusted = list(PvalTable = empty, plot = p0,
                        Missing = data.frame(XVar = character(), YVar = character()))
    ))
  }

  for (variable in names(Data)[vapply(Data, is.logical, logical(1))]) {
    label <- attr(Data[[variable]], "label", exact = TRUE)
    Data[[variable]] <- factor(Data[[variable]], levels = c(FALSE, TRUE))
    attr(Data[[variable]], "label") <- label
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
    for (nm in c("r", "R", "estimate", "cor", "correlation", "ges", "eta2_H", "cramers_v")) {
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

  WarnBranch <- function(branch, variables, e) {
    warning(branch, " branch failed for ", paste(variables, collapse = ", "),
            ": ", conditionMessage(e), call. = FALSE)
    NULL
  }
  df_Requested <- expand.grid(
    XVar = setdiff(OutcomeVars, Covariates),
    YVar = setdiff(PredictorVars, Covariates), stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  ) %>% dplyr::filter(XVar != YVar)

  MiningPairTable <- function(object, type) {
    if (is.data.frame(object$Results)) {
      return(dplyr::transmute(object$Results, XVar = .data$Outcome, YVar = .data$Predictor,
        p = .data$PValue, EffectSize = .data$Estimate, Test = type))
    }
    # Compatibility for older dependency outputs. Current functions expose Results.
    df_Table <- if (type == "ANOVA") object$Unadjusted$PvalTable else if (type == "Correlation") object$Unadjusted$plot$data else object$p$data
    if (!is.data.frame(df_Table) || !nrow(df_Table)) return(data.frame(
      XVar = character(), YVar = character(), p = numeric(), EffectSize = numeric(), Test = character()))
    if (type == "ANOVA") {
      df_Table$XVar <- as.character(df_Table$ContinuousVariable)
      df_Table$YVar <- as.character(df_Table$CategoricalVariable)
    }
    df_Table %>% dplyr::mutate(p = coalesce_p(.), EffectSize = coalesce_effect(.), Test = type) %>%
      dplyr::select(XVar, YVar, p, EffectSize, Test)
  }

  results_list <- list()

  # Correlations
  if (length(num_Outcomes) > 0 && length(num_Predictors) > 0) {

    cor_res <- tryCatch(
      SciDataReportR::PlotCorrelationsHeatmap(
        Data,
        predictor_vars = num_Predictors,
        outcome_vars = num_Outcomes,
        method = method,
        Relabel = FALSE
      ),
      error = function(e) WarnBranch("correlation", c(num_Outcomes, num_Predictors), e)
    )

    if (!is.null(cor_res)) {

      df_cor <- MiningPairTable(cor_res, "Correlation") %>%
        dplyr::select(XVar, YVar, p, EffectSize, Test) %>%
        dplyr::filter(XVar != YVar)

      results_list[[length(results_list) + 1]] <- df_cor
    }
  }

  # Test mixed pairs in both requested orientations.
  for (direction in c("forward", "reverse")) {
    vars_Categorical <- if (direction == "forward") cat_Predictors else cat_Outcomes
    vars_Continuous <- if (direction == "forward") num_Outcomes else num_Predictors
    if (length(vars_Categorical) == 0 || length(vars_Continuous) == 0) next

    anova_res <- tryCatch(
      SciDataReportR::PlotAnovaRelationshipsMatrix(
        Data,
        CatVars = vars_Categorical,
        ContVars = vars_Continuous,
        Relabel = FALSE,
        Parametric = Parametric
      ),
      error = function(e) WarnBranch("ANOVA", c(vars_Categorical, vars_Continuous), e)
    )

    if (!is.null(anova_res)) {

      df_anova <- MiningPairTable(anova_res, "ANOVA") %>%
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
        Relabel = FALSE
      ),
      error = function(e) WarnBranch("ChiSq", c(cat_Outcomes, cat_Predictors), e)
    )

    if (!is.null(chi_res)) {

      df_chi <- MiningPairTable(chi_res, "ChiSq") %>%
        dplyr::select(XVar, YVar, p, EffectSize, Test) %>%
        dplyr::filter(XVar != YVar)

      results_list[[length(results_list) + 1]] <- df_chi
    }
  }

  # Symmetry fix
  results <- dplyr::bind_rows(c(list(data.frame(
      XVar = character(), YVar = character(), p = numeric(),
      EffectSize = numeric(), Test = character())), results_list)) %>%
    dplyr::filter(
      is.finite(p)
    ) %>%
    dplyr::mutate(
      VarA = pmin(XVar, YVar),
      VarB = pmax(XVar, YVar)
    ) %>%
    dplyr::group_by(VarA, VarB) %>%
    dplyr::summarise(
      p = min(p, na.rm = TRUE),
      EffectSize = if (any(is.finite(EffectSize))) max(abs(EffectSize[is.finite(EffectSize)])) else NA_real_,
      Test = dplyr::first(Test),
      .groups = "drop"
    )

  results <- dplyr::bind_rows(
    results %>% dplyr::transmute(XVar = VarA, YVar = VarB, p, EffectSize, Test),
    results %>% dplyr::transmute(XVar = VarB, YVar = VarA, p, EffectSize, Test)
  )

  df_Missing <- dplyr::anti_join(df_Requested, results, by = c("XVar", "YVar"))
  if (nrow(df_Missing) > 0) {
    pairs <- paste(df_Missing$XVar, df_Missing$YVar, sep = " vs ")
    warning("Missing requested pairs: ", paste(utils::head(pairs, 10), collapse = ", "),
            if (length(pairs) > 10) paste0(" and ", length(pairs) - 10, " more") else "",
            call. = FALSE)
  }

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

  results$stars <- ScidrPValueStars(
    ScidrInferenceP(results$p, results$p_adj), IncludeTrendSymbol = IncludeTrendSymbol
  )

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
  results$stars <- factor(results$stars, levels = c("ns", "+", "*", "**", "***"))

  # Sizes are compensated for the differing apparent area of shapes 15-18 so
  # that the marker grows monotonically with significance.
  size_map  <- c("ns" = 2, "+" = 2.7, "*" = 3.5, "**" = 4, "***" = 6)
  shape_map <- c("ns" = 16, "+" = 3, "*" = 17, "**" = 15, "***" = 18)

  results$size_val <- size_map[as.character(results$stars)]

  shape_labels <- c(
    "ns" = if (IncludeTrendSymbol) "ns (p > 0.10)" else "ns (p > 0.05)",
    "+" = "+ (0.05 < p <= 0.10)",
    "*"  = "* (0.01 < p <= 0.05)",
    "**" = "** (0.001 < p <= 0.01)",
    "***"= "*** (p <= 0.001)"
  )
  significance_bands <- if (IncludeTrendSymbol) names(shape_map) else setdiff(names(shape_map), "+")

  p <- ggplot2::ggplot(results, ggplot2::aes(x = XLabel, y = YLabel)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = EffectSizeAbs, size = stars, shape = stars),
      stroke = 1, show.legend = TRUE
    ) +
    ggplot2::scale_colour_gradient(
      low = "grey85", high = "#841B37", limits = c(0, 1),
      trans = "sqrt", oob = scales::squish,
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Effect size", na.value = "grey70") +
    ggplot2::scale_shape_manual(
      values = shape_map, labels = shape_labels[significance_bands], name = "Significance",
      limits = significance_bands, breaks = significance_bands, drop = FALSE,
      guide = ggplot2::guide_legend(override.aes = list(
        size = unname(size_map[significance_bands]), colour = "black", stroke = 1))) +
    ggplot2::scale_size_manual(values = size_map, limits = significance_bands,
                              drop = FALSE, guide = "none") +
    ggplot2::labs(subtitle = if (PVar == "PValue")
      "Significance: raw p-values (no multiple-comparison correction)" else ScidrAdjustedSignificanceSubtitle(fdr_scope)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 7, color = "grey40"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank()
    )

  out <- list(
    Unadjusted = list(PvalTable = results, plot = p, Missing = df_Missing),
    Metadata = list(
      TreatOrdinalAs = TreatOrdinalAs,
      DisplayLabels = labels[unique(c(OutcomeVars, PredictorVars))]
    )
  )
  # Standardized p-value element alias (old name kept)
  out$p <- out$Unadjusted
  out
}
