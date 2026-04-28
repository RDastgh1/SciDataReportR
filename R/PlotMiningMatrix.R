#' PlotMiningMatrix
#'
#' Generate a matrix of statistical relationships between variables.
#'
#' @param Data A data frame.
#' @param OutcomeVars Outcome variables.
#' @param PredictorVars Predictor variables. If NULL, uses OutcomeVars.
#' @param Covariates Optional covariates (reserved for future use).
#' @param Relabel Use labels instead of names.
#' @param Parametric Use parametric tests.
#'
#' @return List with tables and plots.
#' @export
PlotMiningMatrix <- function(
    Data,
    OutcomeVars,
    PredictorVars = NULL,
    Covariates = NULL,
    Relabel = TRUE,
    Parametric = TRUE
) {

  method <- ifelse(Parametric, "pearson", "spearman")

  OutcomeVars <- unique(intersect(as.character(OutcomeVars), names(Data)))

  if (is.null(PredictorVars)) {
    PredictorVars <- OutcomeVars
  } else {
    PredictorVars <- unique(intersect(as.character(PredictorVars), names(Data)))
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) {
    empty <- data.frame()
    p0 <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = p0)
    ))
  }

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

  # FORCE LABEL COMPLETENESS
  labels <- sjlabelled::get_label(Data)
  names(labels) <- names(Data)

  missing_idx <- is.na(labels) | labels == ""
  labels[missing_idx] <- names(labels)[missing_idx]

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
        xVars = num_Predictors,
        yVars = num_Outcomes,
        covars = Covariates,
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
        Covariates = Covariates,
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
        xVars = cat_Predictors,
        yVars = cat_Outcomes,
        covars = Covariates,
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

  # 🔥 SYMMETRY FIX
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

  results$p_adj <- stats::p.adjust(results$p, method = "fdr")

  results <- results %>%
    rstatix::add_significance(p.col = "p", output.col = "stars")

  # 🔥 SAFE LABEL APPLICATION
  if (Relabel) {
    results$XLabel <- safe_lookup(results$XVar, labels)
    results$YLabel <- safe_lookup(results$YVar, labels)
  } else {
    results$XLabel <- results$XVar
    results$YLabel <- results$YVar
  }

  # 🔥 SAFE ORDERING (no NA levels ever)
  x_order <- if (Relabel) safe_lookup(OutcomeVars, labels) else OutcomeVars
  y_order <- if (Relabel) safe_lookup(PredictorVars, labels) else PredictorVars

  results$XLabel <- factor(results$XLabel, levels = unique(x_order))
  results$YLabel <- factor(results$YLabel, levels = rev(unique(y_order)))

  # FINAL CLEANUP
  results <- results %>%
    dplyr::filter(!is.na(XLabel), !is.na(YLabel))

  results$EffectSizeAbs <- abs(results$EffectSize)

  size_map  <- c("ns" = 2, "*" = 3, "**" = 4, "***" = 5)
  shape_map <- c("ns" = 16, "*" = 17, "**" = 15, "***" = 18)

  results$size_val <- size_map[as.character(results$stars)]

  shape_labels <- c(
    "ns" = "ns (p ≥ 0.05)",
    "*"  = "* (p < 0.05)",
    "**" = "** (p < 0.01)",
    "***"= "*** (p < 0.001)"
  )

  p <- ggplot2::ggplot(results, ggplot2::aes(x = XLabel, y = YLabel)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = EffectSizeAbs, size = size_val, shape = stars)
    ) +
    viridis::scale_color_viridis(option = "plasma", limits = c(0,1), name = "Effect Size", direction = -1) +
    ggplot2::scale_shape_manual(values = shape_map, labels = shape_labels, name = "Significance") +
    ggplot2::scale_size_continuous(range = c(2, 5), guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank()
    )

  list(
    Unadjusted = list(PvalTable = results, plot = p)
  )
}
