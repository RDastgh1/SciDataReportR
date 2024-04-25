#' Plot ANOVA Relationships Matrix
#'
#' This function plots the relationship between continuous and categorical variables
#' using ANOVA or Kruskal-Wallis tests. It generates scatter plots with points
#' colored and shaped based on statistical significance and effect size.
#'
#' @param Data The data frame containing the variables of interest.
#' @param CatVars Character vector of categorical variable names.
#' @param ContVars Character vector of continuous variable names.
#' @param Covariates Optional character vector of covariate names for ANCOVA analysis.
#' @param Relabel Logical indicating whether to relabel variables with their labels (default is TRUE).
#' @param Parametric Logical indicating whether to use parametric (ANOVA) or non-parametric (Kruskal-Wallis) tests (default is TRUE).
#' @return A list containing three ggplot objects: p (scatter plot without multiple comparison correction), p_FDR (scatter plot with FDR correction), and pvaltable (data frame of p-values and significance).
#' @import dplyr
#' @importFrom ggplot2 aes geom_point labs scale_shape_manual scale_color_gradientn guides theme
#' @importFrom sjlabelled get_label
#' @importFrom tidyr pivot_longer
#' @importFrom rstatix anova_test kruskal_test get_summary_stats add_significance adjust_pvalue
#' @importFrom stats var
#' @importFrom utils na.omit
#' @export
PlotAnovaRelationshipsMatrix <- function(Data, CatVars, ContVars, Covariates = NULL, Relabel = TRUE, Parametric = TRUE) {

  # Subset data
  DataSubset <- Data[c(CatVars, ContVars, Covariates)]

  # Pivot data into long format
  mData <- pivot_longer(DataSubset, cols = all_of(CatVars), names_to = "CategoricalVariable", values_to = "CategoricalValue")
  mData <- pivot_longer(mData, cols = all_of(ContVars), names_to = "ContinuousVariable", values_to = "ContinuousValue")

  # Convert to factors and remove NA values
  mData$ContinuousVariable <- as.factor(mData$ContinuousVariable)
  mData$CategoricalValue <- as.factor(mData$CategoricalValue)
  mData <- na.omit(mData)

  # Filter out variables with only one level
  nlevels.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>% mutate(n = length(unique(CategoricalValue)))
  mData <- mData[nlevels.test$n > 1,]

  # Perform ANOVA or Kruskal-Wallis tests
  if (Parametric) {
    stat.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
      anova_test(as.formula(paste("ContinuousValue ~ CategoricalValue", if (!is.null(Covariates)) paste("+", Covariates, collapse = "+")))) %>%
      add_significance() %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance()
  } else {
    stat.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
      kruskal_test(as.formula(paste("ContinuousValue ~ CategoricalValue", if (!is.null(Covariates)) paste("+", Covariates, collapse = "+")))) %>%
      add_significance() %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance()
  }

  # Calculate summary statistics
  summstats <- mData %>% group_by(ContinuousVariable, CategoricalVariable, CategoricalValue) %>% get_summary_stats() %>% filter(variable == "ContinuousValue")

  # Calculate fold change for two-group comparisons
  ngroups <- length(unique(summstats$CategoricalValue))
  if (ngroups == 2) {
    FCStats <- summstats %>% select(CategoricalVariable, CategoricalValue, ContinuousVariable, n, mean, sd) %>%
      pivot_wider(names_from = CategoricalValue, values_from = n:sd)
    Groups <- levels(summstats$CategoricalValue)
    GroupMeanLabels <- paste0("mean_", Groups)
    FCStats$FoldChange <- FCStats[[GroupMeanLabels[2]]] / FCStats[[GroupMeanLabels[1]]]
    FCStats$Log2FC <- log2(FCStats$FoldChange)
  } else {
    FCStats <- data.frame()
  }

  # Modify significance levels for plotting
  stat.test$logp <- -log10(stat.test$p)
  stat.test$`p<.05` <- stat.test$p.signif
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns", "*", "**", "***", "****"))
  stat.test$logp_FDR <- -log10(stat.test$p.adj)
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns", "*", "**", "***", "****"))
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif, levels = c("ns", "*", "**", "***", "****"))

  # Retrieve variable labels if Relabel is TRUE
  if (Relabel) {
    Data<- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[stat.test$CategoricalVariable], def.value = colnames(Data[stat.test$CategoricalVariable])) %>% as.data.frame() %>% rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")
    ylabels <- sjlabelled::get_label(Data[stat.test$ContinuousVariable], def.value = colnames(Data[stat.test$ContinuousVariable])) %>% as.data.frame() %>% rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")
    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  # Create plot text for hover-over display
  PlotText <- paste("</br>Cat Var Label:", stat.test$XLabel,
                    "</br> Cat Var:", stat.test$CategoricalVariable,
                    "</br> Cont Var Label:", stat.test$YLabel,
                    "</br> Cont Var:", stat.test$ContinuousVariable,
                    "</br> P-Value: ", stat.test$p, stat.test$p.signif,
                    "</br> FDR-corrected P: ", stat.test$p.adj, stat.test$p.adj.signif)
  stat.test$PlotText <- PlotText

  # Convert variable labels to factors if Relabel is TRUE
  if (Relabel) {
    stat.test$ContinuousVariable <- factor(stat.test$ContinuousVariable, levels = ContVars, labels = sjlabelled::get_label(Data[ContVars], def.value = colnames(Data[ContVars])))
    stat.test$CategoricalVariable <- factor(stat.test$CategoricalVariable, levels = CatVars, labels = sjlabelled::get_label(Data[CatVars], def.value = colnames(Data[CatVars])))
  } else {
    stat.test$ContinuousVariable <- factor(stat.test$ContinuousVariable, levels = ContVars)
    stat.test$CategoricalVariable <- factor(stat.test$CategoricalVariable, levels = CatVars)
  }

  # Create scatter plot without FDR correction
  p <- ggplot(stat.test, aes(y = ContinuousVariable, x = CategoricalVariable, shape = `p<.05`, text = PlotText)) +
    geom_point(aes(size = logp, colour = p)) +
    theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    scale_color_gradientn(trans = "log", colours = rev(brewer.pal(9, "RdPu")), limits = c(0.0000001,1), oob = scales::squish, breaks = c(1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)) +
    guides(size = FALSE) +
    labs(subtitle = "No Multiple Comparison Correction") + xlab("") + ylab("")

  # Create scatter plot with FDR correction
  p_FDR <- ggplot(stat.test, aes(y = ContinuousVariable, x = CategoricalVariable, shape = p.adj.signif, text = PlotText)) +
    geom_point(aes(size = logp_FDR, colour = p.adj)) +
    theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    scale_color_gradientn(trans = "log", colours = rev(brewer.pal(9, "RdPu")), limits = c(0.0000001,1), oob = scales::squish, breaks = c(1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)) +
    guides(size = FALSE) +
    labs(subtitle = "FDR Correction") +
    scale_x_discrete(limits = levels(stat.test$CategoricalVariable)) +
    scale_y_discrete(limits = levels(stat.test$ContinuousVariable))

  # Return the plots and p-value table
  return(list(p = p, pvaltable = pvaltable, p_FDR = p_FDR, FCStats = FCStats))
}
