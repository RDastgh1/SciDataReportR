#' Plot ANOVA Relationships Matrix
#'
#' This function plots the relationship between continuous and categorical variables
#' using ANOVA or Kruskal-Wallis tests. It generates a "heatmap" with points
#' colored and shaped based on statistical significance and effect size.
#'
#' @param Data The data frame containing the variables of interest.
#' @param CatVars Character vector of categorical variable names.
#' @param ContVars Character vector of continuous variable names.
#' @param Covariates Optional character vector of covariate names for ANCOVA analysis.
#' @param Relabel Logical indicating whether to relabel variables with their labels (default is TRUE).
#' @param Ordinal Logical, indicating whether ordinal variables should be considered.
#' @param Parametric Logical indicating whether to use parametric (ANOVA) or non-parametric (Kruskal-Wallis) tests (default is TRUE).
#' @return A list containing three ggplot objects: p (scatter plot without multiple comparison correction), p_FDR (scatter plot with FDR correction), and pvaltable (data frame of p-values and significance).
#' @import dplyr
#' @importFrom ggplot2 aes geom_point labs scale_shape_manual scale_color_gradientn guides theme
#' @importFrom sjlabelled get_label
#' @importFrom tidyr pivot_longer
#' @importFrom rstatix anova_test kruskal_test get_summary_stats add_significance adjust_pvalue
#' @importFrom stats var
#' @importFrom utils na.omit
#' @importFrom rstatix add_significance adjust_pvalue anova_test kruskal_test get_summary_stats
#' @export
PlotAnovaRelationshipsMatrix <- function(
    Data,
    CatVars,
    ContVars,
    Covariates = NULL,
    Relabel = TRUE,
    Parametric = TRUE,
    Ordinal = FALSE
) {


  CatVars <- unique(as.character(CatVars))
  ContVars <- unique(as.character(ContVars))
  Covariates <- if (!is.null(Covariates) && length(Covariates) > 0) unique(as.character(Covariates)) else character(0)

  vars_keep <- unique(c(CatVars, ContVars, Covariates))
  missing <- setdiff(vars_keep, names(Data))
  if (length(missing) > 0) {
    stop("PlotAnovaRelationshipsMatrix: Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }


  df0 <- Data |>
    dplyr::select(dplyr::all_of(vars_keep)) |>
    dplyr::mutate(.rowID = dplyr::row_number())


  cov_df <- if (length(Covariates) > 0) {
    df0 |> dplyr::select(.rowID, dplyr::all_of(Covariates))
  } else {
    df0 |> dplyr::select(.rowID)
  }

  # Factors for categorical predictors
  df0 <- df0 |>
    dplyr::mutate(dplyr::across(dplyr::all_of(CatVars), ~ as.factor(.x)))

  # Pivot CatVars long
  mData <- df0 |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(CatVars),
      names_to = "CategoricalVariable",
      values_to = "CategoricalValue"
    ) |>
    # Pivot ContVars long
    tidyr::pivot_longer(
      cols = dplyr::all_of(ContVars),
      names_to = "ContinuousVariable",
      values_to = "ContinuousValue"
    ) |>
    # drop the original covariate columns if they were carried through pivoting
    dplyr::select(-dplyr::any_of(Covariates)) |>
    # reattach covariates deterministically
    dplyr::left_join(cov_df, by = ".rowID") |>
    dplyr::mutate(
      ContinuousVariable = as.factor(ContinuousVariable),
      CategoricalValue = as.factor(CategoricalValue),
      ContinuousValue = suppressWarnings(as.numeric(as.character(ContinuousValue)))
    ) |>
    tidyr::drop_na(ContinuousValue, CategoricalValue)


  mData <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
    dplyr::filter(dplyr::n_distinct(CategoricalValue) > 1) |>
    dplyr::filter(stats::sd(ContinuousValue, na.rm = TRUE) > 0) |>
    dplyr::ungroup()


  rhs <- "CategoricalValue"
  if (length(Covariates) > 0) rhs <- paste(rhs, paste(Covariates, collapse = " + "), sep = " + ")
  fml <- stats::as.formula(paste("ContinuousValue ~", rhs))

  if (Parametric) {
    method <- "Anova"
    stat.test <- mData |>
      dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
      rstatix::anova_test(fml) |>
      rstatix::add_significance() |>
      rstatix::adjust_pvalue(method = "fdr") |>
      rstatix::add_significance()
  } else {
    method <- "Kruskal"
    stat.test <- mData |>
      dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
      rstatix::kruskal_test(fml) |>
      rstatix::add_significance() |>
      rstatix::adjust_pvalue(method = "fdr") |>
      rstatix::add_significance()
  }


  summstats <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable, CategoricalValue) |>
    rstatix::get_summary_stats() |>
    dplyr::filter(variable == "ContinuousValue")

  ngroups <- length(unique(summstats$CategoricalValue))
  if (ngroups == 2) {
    FCStats <- summstats |>
      dplyr::select(CategoricalVariable, CategoricalValue, ContinuousVariable, n, mean, sd) |>
      tidyr::pivot_wider(names_from = CategoricalValue, values_from = n:sd)

    Groups <- levels(summstats$CategoricalValue)
    GroupMeanLabels <- paste0("mean_", Groups)
    FCStats$FoldChange <- FCStats[[GroupMeanLabels[2]]] / FCStats[[GroupMeanLabels[1]]]
    FCStats$Log2FC <- log2(FCStats$FoldChange)
  } else {
    FCStats <- data.frame()
  }


  stat.test$logp <- -log10(stat.test$p)
  stat.test$`p<.05` <- factor(stat.test$p.signif, levels = c("ns","*","**","***","****"))
  stat.test$logp_FDR <- -log10(stat.test$p.adj)
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif, levels = c("ns","*","**","***","****"))

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)

    xlabels <- sjlabelled::get_label(Data[as.character(stat.test$CategoricalVariable)],
                                     def.value = stat.test$CategoricalVariable) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(stat.test$ContinuousVariable)],
                                     def.value = stat.test$ContinuousVariable) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  stat.test$PlotText <- paste(
    "</br>Cat Var Label:", stat.test$XLabel,
    "</br> Cat Var:", stat.test$CategoricalVariable,
    "</br> Cont Var Label:", stat.test$YLabel,
    "</br> Cont Var:", stat.test$ContinuousVariable,
    "</br> P-Value: ", stat.test$p, stat.test$p.signif,
    "</br> FDR-corrected P: ", stat.test$p.adj, stat.test$p.adj.signif,
    "</br> GES Effect size: ", stat.test$ges,
    "</br> npairs: ", stat.test$DFd
  )


  stat.test <- stat.test |>
    dplyr::ungroup() |>
    as.data.frame() |>
    dplyr::filter(Effect == "CategoricalValue")


  p <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = `p<.05`, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = `p<.05`, colour = ges)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::scale_color_gradientn(
      transform = "log1p",
      colours = rev(paletteer::paletteer_c("grDevices::Viridis", n = 20)),
      values = scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)),
      limits = c(0, 0.3),
      breaks = c(0.01, 0.06, 0.14, 0.3),
      labels = c("0.01: Small", "0.06: Medium", "0.14: Large", "0.3-1"),
      oob = scales::squish
    ) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  p_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = p.adj.signif, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = p.adj.signif, colour = ges)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::scale_color_gradientn(
      transform = "log1p",
      colours = rev(paletteer::paletteer_c("grDevices::Viridis", n = 20)),
      values = scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)),
      limits = c(0, 0.3),
      breaks = c(0.01, 0.06, 0.14, 0.3),
      labels = c("0.01: Small", "0.06: Medium", "0.14: Large", "0.3-1"),
      oob = scales::squish
    ) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  M <- list(PvalTable = stat.test |> dplyr::mutate(Test = method), plot = p)
  M_FDR <- list(PvalTable = stat.test |> dplyr::mutate(Test = method), plot = p_FDR)

  list(Unadjusted = M, FDRCorrected = M_FDR, method = method, Relabel = Relabel, Covariates = Covariates)
}
