#' PlotMiningMatrix
#'
#' This function generates a matrix of statistical relationships between specified outcome and predictor variables
#' in a dataset. It includes visualizations for correlations, ANOVA results, and FDR corrections.
#'
#' @param Data A data frame containing the dataset to analyze.
#' @param OutcomeVars A vector of outcome variables to be analyzed.
#' @param PredictorVars A vector of predictor variables to be analyzed.
#' @param Covariates An optional vector of covariates to adjust the analysis (default is NULL).
#' @param Relabel Logical flag indicating whether to relabel the variables in the output (default is TRUE).
#' @param Parametric Logical flag indicating whether to use parametric methods (default is TRUE).
#' If FALSE, non-parametric methods will be used.
#'
#' @return A list containing the following elements:
#' \item{Unadjusted}{A list with the unadjusted p-value table and corresponding plot.}
#' \item{FDRCorrected}{A list with the FDR-adjusted p-value table and corresponding plot.}
#' \item{method}{The method used for correlation ("pearson" for parametric, "spearman" for non-parametric).}
#' \item{Relabel}{The value of the Relabel parameter.}
#' \item{Covariates}{The covariates used in the analysis, if any.}
#'
#' @import dplyr ggplot2 paletteer rstatix sjlabelled
#' @importFrom stats p.adjust
#' @importFrom magrittr %>%
#'
#'
#' @export
PlotMiningMatrix <- function(Data, OutcomeVars, PredictorVars, Covariates = NULL, Relabel = TRUE, Parametric = TRUE) {
  # Reverse order of OutcomeVars for the y axis
  OutcomeVars <- rev(OutcomeVars)
  # Determine correlation method based on Parametric flag
  method <- ifelse(Parametric, "pearson", "spearman")

  # Separate Outcome and Predictor Variables into categorical and numeric
  num_Outcomes <- getNumVars(Data %>% select(all_of(OutcomeVars)))
  cat_Outcomes <- getCatVars(Data %>% select(all_of(OutcomeVars)))

  num_Predictors <- getNumVars(Data %>% select(all_of(PredictorVars)))
  cat_Predictors <- getCatVars(Data %>% select(all_of(PredictorVars)))

  # Continuous Predictors and Continuous Outcomes: Correlation heatmap
  if(length(num_Outcomes) >0 & length(num_Predictors)>0){
    P_Correlations <- PlotCorrelationsHeatmap(Data, num_Predictors, num_Outcomes, covars = Covariates, Relabel = Relabel)
    P_Correlations <- P_Correlations$Unadjusted$plot$data %>%
      select(-P_adj) %>%
      mutate(Test = method, p = P)


    P_Correlations <- P_Correlations %>%
      dplyr::mutate(
        temp = XLabel,     # Create a temporary column
        XLabel = YLabel,   # Switch XLabel with YLabel
        YLabel = temp      # Assign temp to YLabel
      ) %>%
      dplyr::select(-temp)


  }else{
    P_Correlations<- NULL}




  # Categorical Predictors and Continuous Outcomes: ANOVA or t-test
  if(length(cat_Predictors)>0 & length(num_Outcomes)>0){
    P_Anova1 <- PlotAnovaRelationshipsMatrix(Data, cat_Predictors, num_Outcomes, Covariates = Covariates, Parametric = Parametric)
    P_Anova1 <- P_Anova1$Unadjusted$PvalTable %>%
      ungroup() %>%
      as.data.frame() %>%
      select(-p.adj, -p.adj.signif, -logp_FDR)
    P_Anova1$nPairs <- P_Anova1$DFd

    # Switch Order

    P_Anova1 <- P_Anova1 %>%
      dplyr::mutate(
        temp = XLabel,     # Create a temporary column
        XLabel = YLabel,   # Switch XLabel with YLabel
        YLabel = temp      # Assign temp to YLabel
      ) %>%
      dplyr::select(-temp)


  }else{
    P_Anova1 <- NULL
  }


  # Continuous Predictors and Categorical Outcomes: ANOVA or t-test
  if(length(cat_Outcomes)>0 & length(num_Predictors)>0){
    P_Anova2 <- PlotAnovaRelationshipsMatrix(Data,  cat_Outcomes, num_Predictors, Covariates = Covariates, Parametric = Parametric)
    P_Anova2 <- P_Anova2$Unadjusted$PvalTable %>%
      ungroup() %>%
      as.data.frame() %>%
      select(-p.adj, -p.adj.signif, -logp_FDR)
    P_Anova2$nPairs <- P_Anova2$DFd
  }else{
    P_Anova2 <- NULL
  }

  # Categorical Predictors and Categorical Outcomes: ChiSquared
  if(length(cat_Outcomes)>0 & length(cat_Predictors)>0){
    P_Chi <- PlotChiSqCovar(Data, cat_Predictors, cat_Outcomes, covars = Covariates)
    P_Chi <- P_Chi$p$data
    P_Chi$p <- P_Chi$pval
    P_Chi$p.adj <- P_Chi$pval.adj


  }else{
    # P_Anova2 <- NULL
    P_Chi <- NULL
  }


  # Combine Correlations and ANOVA results
  dfs <- list(P_Correlations, P_Anova1, P_Anova2, P_Chi)

  # Filter out NULL dataframes
  non_null_dfs <- dfs[!sapply(dfs, is.null)]

  # Merge the remaining dataframes
  P_Combined <- Reduce(function(x, y) full_join(x, y), non_null_dfs)

  P_Combined <- P_Combined %>%
    select(-any_of(c("stars", "stars_FDR", "CategoricalVariable", "ContinuousVariable",
           "p.signif", "PlotText", "ges", "Effect")))

  # Apply FDR Adjustment
  P_Combined$p_adj <- p.adjust(P_Combined$p, method = "fdr")

  # Add significance stars for p-values and FDR-adjusted p-values
  P_Combined <- P_Combined %>%
    rstatix::add_significance(p.col = "p", output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr")

  # Reorder factors and calculate log-transformed p-values
  P_Combined$stars_fdr <- factor(P_Combined$stars_fdr, levels = c("ns", "*", "**", "***", "****"))
  P_Combined$stars <- factor(P_Combined$stars, levels = c("ns", "*", "**", "***", "****"))
  P_Combined$logp <- -log10(P_Combined$p)
  P_Combined$logpfdr <- -log10(P_Combined$p_adj)

  # Reorder labels for plotting
  yorder <- sjlabelled::get_label(Data %>% select(as.character(PredictorVars)), def.value = PredictorVars)
  xorder <- sjlabelled::get_label(Data %>% select(as.character(OutcomeVars)), def.value = OutcomeVars)
  P_Combined$XLabel <- factor(P_Combined$XLabel, levels = xorder)
  P_Combined$YLabel <- factor(P_Combined$YLabel, levels = yorder)
  P_Combined<- P_Combined%>% filter(!is.na(p))
  # Create plot for unadjusted p-values
  p <- ggplot(P_Combined, aes(y = YLabel, x = XLabel)) +
    geom_point(aes(size = stars, colour = stars), shape = 18) +
    scale_color_manual(values = rev(paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5))) +
    guides(size = "none") +
    labs(subtitle = "No Multiple Comparison Correction") +
    xlab("") + ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())

  # Create plot for FDR-corrected p-values
  p_FDR <- ggplot(P_Combined, aes(y = YLabel, x = XLabel)) +
    geom_point(aes(size = stars_fdr, colour = stars_fdr), shape = 18) +
    scale_color_manual(values = rev(paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5))) +
    guides(size = "none") +
    labs(subtitle = "FDR Correction") +
    xlab("") + ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())

  # Prepare output list
  M <- list(PvalTable = P_Combined, plot = p)
  M_FDR <- list(PvalTable = P_Combined, plot = p_FDR)

  return(list(Unadjusted = M, FDRCorrected = M_FDR, method = method, Relabel = Relabel, Covariates = Covariates))
}
