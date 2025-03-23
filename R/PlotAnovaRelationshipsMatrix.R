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
PlotAnovaRelationshipsMatrix<- function (Data, CatVars, ContVars, Covariates = NULL, Relabel = TRUE,
          Parametric = TRUE, Ordinal = FALSE)
{
  DataSubset <- Data[c(CatVars, ContVars, Covariates)]
  DataSubset[CatVars] <- lapply(DataSubset[CatVars], factor,
                                ordered = FALSE) %>% as.data.frame()
  mData <- pivot_longer(DataSubset, cols = all_of(CatVars),
                        names_to = "CategoricalVariable", values_to = "CategoricalValue")
  mData <- pivot_longer(mData, cols = all_of(ContVars), names_to = "ContinuousVariable",
                        values_to = "ContinuousValue")
  mData$ContinuousVariable <- as.factor(mData$ContinuousVariable)
  mData$CategoricalValue <- as.factor(mData$CategoricalValue)
  mData <- mData %>%
    mutate(ContinuousValue = as.numeric(as.character(ContinuousValue)))

  mData <- na.omit(mData)
  nlevels.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
    mutate(n = length(unique(CategoricalValue)))
  mData <- mData[nlevels.test$n > 1, ]
  nvar.test <- mData %>% group_by(ContinuousVariable) %>% mutate(v = var(ContinuousValue))
  mData <- mData[nvar.test$v > 0, ]

  # Remove value-group combinations with only one unique value
  mData <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable) %>%
    filter(n_distinct(ContinuousValue) > 1) %>%
    ungroup()



  # Remove value-group combinations with no variability
  mData <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable) %>%
    filter(sd(ContinuousValue, na.rm = TRUE) > 0) %>%
    ungroup()

  mData <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable) %>%
    filter(sd(ContinuousValue, na.rm = TRUE) > 0) %>%  # Ensure ContinuousValue varies
    filter(n_distinct(CategoricalValue, na.rm = TRUE) > 1) %>%  # Ensure CategoricalValue varies
    ungroup()

  # check the variance of ContinuousValue within each CategoricalValue in every (ContinuousVariable, CategoricalVariable) group.
  # may need to adjust this so it's not so stringent.
  mData <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable, CategoricalValue) %>%
    filter(sd(ContinuousValue, na.rm = TRUE) > 0) %>%
    ungroup()

  mData <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable) %>%
    filter(n_distinct(CategoricalValue) > 1) %>%
    ungroup()

  if (Parametric) {
    method <- "Anova"
    stat.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
      rstatix::anova_test(as.formula(paste("ContinuousValue ~ CategoricalValue",
                                           if (!is.null(Covariates))
                                             paste("+", Covariates, collapse = "+")))) %>%
      rstatix::add_significance() %>% rstatix::adjust_pvalue(method = "fdr") %>%
      rstatix::add_significance()
  }else {
    method <- "Kruskal"
    stat.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
      rstatix::kruskal_test(as.formula(paste("ContinuousValue ~ CategoricalValue",
                                             if (!is.null(Covariates))
                                               paste("+", Covariates, collapse = "+")))) %>%
      rstatix::add_significance() %>% rstatix::adjust_pvalue(method = "fdr") %>%
      rstatix::add_significance()
  }
  summstats <- mData %>% group_by(ContinuousVariable, CategoricalVariable,
                                  CategoricalValue) %>% rstatix::get_summary_stats() %>%
    filter(variable == "ContinuousValue")
  ngroups <- length(unique(summstats$CategoricalValue))
  if (ngroups == 2) {
    FCStats <- summstats %>% select(CategoricalVariable,
                                    CategoricalValue, ContinuousVariable, n, mean, sd) %>%
      pivot_wider(names_from = CategoricalValue, values_from = n:sd)
    Groups <- levels(summstats$CategoricalValue)
    GroupMeanLabels <- paste0("mean_", Groups)
    FCStats$FoldChange <- FCStats[[GroupMeanLabels[2]]]/FCStats[[GroupMeanLabels[1]]]
    FCStats$Log2FC <- log2(FCStats$FoldChange)
  }else {
    FCStats <- data.frame()
  }
  stat.test$logp <- -log10(stat.test$p)
  stat.test$`p<.05` <- stat.test$p.signif
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns",
                                                            "*", "**", "***", "****"))
  stat.test$logp_FDR <- -log10(stat.test$p.adj)
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns",
                                                            "*", "**", "***", "****"))
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif,
                                   levels = c("ns", "*", "**", "***", "****"))
  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[as.character(stat.test$CategoricalVariable)],
                                     def.value = stat.test$CategoricalVariable) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")
    ylabels <- sjlabelled::get_label(Data[as.character(stat.test$ContinuousVariable)],
                                     def.value = stat.test$ContinuousVariable) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")
    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }
  PlotText <- paste("</br>Cat Var Label:", stat.test$XLabel,
                    "</br> Cat Var:", stat.test$CategoricalVariable, "</br> Cont Var Label:",
                    stat.test$YLabel, "</br> Cont Var:", stat.test$ContinuousVariable,
                    "</br> P-Value: ", stat.test$p, stat.test$p.signif, "</br> FDR-corrected P: ",
                    stat.test$p.adj, stat.test$p.adj.signif, "</br> GES Effect size: ",
                    stat.test$ges, "</br> npairs: ", stat.test$DFd)
  stat.test$PlotText <- PlotText
  if (Relabel) {
    stat.test$YLabel <- factor(stat.test$ContinuousVariable,
                               levels = ContVars, labels = sjlabelled::get_label(Data[as.character(ContVars)],
                                                                                 def.value = ContVars), )
    stat.test$XLabel <- factor(stat.test$CategoricalVariable,
                               levels = CatVars, labels = sjlabelled::get_label(Data[as.character(CatVars)],
                                                                                def.value = CatVars), )
  }else {
    stat.test$YLabel <- factor(stat.test$ContinuousVariable,
                               levels = ContVars)
    stat.test$XLabel <- factor(stat.test$CategoricalVariable,
                               levels = CatVars)
  }
  stat.test <- stat.test %>% ungroup() %>% ungroup() %>% as.data.frame() %>%
    dplyr::filter(Effect == "CategoricalValue")
  p <- ggplot(stat.test, aes(y = YLabel, x = XLabel, shape = `p<.05`,
                             text = PlotText)) + geom_point(aes(size = `p<.05`, colour = ges)) +
    theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
          axis.title.y = element_blank()) + scale_shape_manual(values = c(7,
                                                                          16, 17, 15, 18), drop = FALSE) + scale_color_gradientn(transform = "log1p",
                                                                                                                                 colours = rev(paletteer::paletteer_c("grDevices::Viridis",
                                                                                                                                                                      n = 20)), values = scales::rescale(c(0, 0.01, 0.05,
                                                                                                                                                                                                           0.1, 0.2, 0.5, 1)), limits = c(0, 0.3), breaks = c(0.01,
                                                                                                                                                                                                                                                              0.06, 0.14, 0.3), labels = c("0.01: Small", "0.06: Medium",
                                                                                                                                                                                                                                                                                           "0.14: Large", "0.3-1"), oob = scales::squish) +
    guides(size = "none") + labs(subtitle = "No Multiple Comparison Correction") +
    xlab("") + ylab("")
  p_FDR <- ggplot(stat.test, aes(y = YLabel, x = XLabel, shape = p.adj.signif,
                                 text = PlotText)) + geom_point(aes(size = p.adj.signif,
                                                                    colour = ges)) + theme(axis.text.x = element_text(angle = 90),
                                                                                           axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    scale_color_gradientn(transform = "log1p", colours = rev(paletteer::paletteer_c("grDevices::Viridis",
                                                                                    n = 20)), values = scales::rescale(c(0, 0.01, 0.05,
                                                                                                                         0.1, 0.2, 0.5, 1)), limits = c(0, 0.3), breaks = c(0.01,
                                                                                                                                                                            0.06, 0.14, 0.3), labels = c("0.01: Small", "0.06: Medium",
                                                                                                                                                                                                         "0.14: Large", "0.3-1"), oob = scales::squish) +
    guides(size = "none") + labs(subtitle = "FDR Correction") +
    xlab("") + ylab("")
  M <- list()
  M$PvalTable <- stat.test %>% ungroup() %>% as.data.frame() %>%
    mutate(Test = method)
  M$plot <- p
  M_FDR <- list()
  M_FDR$PvalTable <- stat.test %>% ungroup() %>% as.data.frame() %>%
    mutate(Test = method)
  M_FDR$plot <- p_FDR
  return(list(Unadjusted = M, FDRCorrected = M_FDR, method = method,
              Relabel = Relabel, Covariates = Covariates))
}
