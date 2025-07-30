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
#' @import dplyr ggplot2 rstatix sjlabelled paletteer
#' @importFrom stats p.adjust
#' @importFrom magrittr %>%
#' @export
PlotMiningMatrix <- function(Data,
                             OutcomeVars,
                             PredictorVars,
                             Covariates = NULL,
                             Relabel = TRUE,
                             Parametric = TRUE) {

  OutcomeVars <- rev(OutcomeVars)
  method <- ifelse(Parametric, "pearson", "spearman")

  # split numeric vs categorical
  num_Outcomes <- getNumVars(Data %>% select(all_of(OutcomeVars)))
  cat_Outcomes <- getCatVars(Data %>% select(all_of(OutcomeVars)))
  num_Predictors <- getNumVars(Data %>% select(all_of(PredictorVars)))
  cat_Predictors <- getCatVars(Data %>% select(all_of(PredictorVars)))

  # correlations
  if (length(num_Outcomes) > 0 && length(num_Predictors) > 0) {
    P_Correlations <- PlotCorrelationsHeatmap(
      Data, num_Predictors, num_Outcomes,
      covars = Covariates, Relabel = Relabel
    )$Unadjusted$plot$data %>%
      select(-P_adj) %>%
      mutate(Test = method, p = P) %>%
      rename(XLabel = YLabel, YLabel = XLabel)
  } else {
    P_Correlations <- NULL
  }

  # ANOVA for categorical predictors x numeric outcomes
  if (length(cat_Predictors) > 0 && length(num_Outcomes) > 0) {
    P_Anova1 <- PlotAnovaRelationshipsMatrix(
      Data, cat_Predictors, num_Outcomes,
      Covariates = Covariates, Parametric = Parametric
    )$Unadjusted$PvalTable %>%
      ungroup() %>% as.data.frame() %>%
      select(-p.adj, -p.adj.signif, -logp_FDR) %>%
      mutate(nPairs = DFd) %>%
      rename(XLabel = YLabel, YLabel = XLabel)
  } else {
    P_Anova1 <- NULL
  }

  # ANOVA for categorical outcomes x numeric predictors
  if (length(cat_Outcomes) > 0 && length(num_Predictors) > 0) {
    P_Anova2 <- PlotAnovaRelationshipsMatrix(
      Data, cat_Outcomes, num_Predictors,
      Covariates = Covariates, Parametric = Parametric
    )$Unadjusted$PvalTable %>%
      ungroup() %>% as.data.frame() %>%
      select(-p.adj, -p.adj.signif, -logp_FDR) %>%
      mutate(nPairs = DFd)
  } else {
    P_Anova2 <- NULL
  }

  # Chi-square for categorical x categorical
  if (length(cat_Outcomes) > 0 && length(cat_Predictors) > 0) {
    P_Chi <- PlotChiSqCovar(
      Data, cat_Predictors, cat_Outcomes,
      covars = Covariates
    )$p$data %>%
      mutate(p = pval, p.adj = pval.adj)
  } else {
    P_Chi <- NULL
  }

  # combine all non-null
  non_null_dfs <- list(P_Correlations, P_Anova1, P_Anova2, P_Chi) %>%
    keep(~ !is.null(.x))
  P_Combined <- Reduce(full_join, non_null_dfs)

  # drop unused columns
  P_Combined <- P_Combined %>%
    select(-any_of(c(
      "stars","stars_FDR","CategoricalVariable","ContinuousVariable",
      "p.signif","PlotText","ges","Effect"
    )))

  # FDR correction and significance annotation
  P_Combined <- P_Combined %>%
    mutate(p_adj = p.adjust(p, method = "fdr")) %>%
    rstatix::add_significance(p.col = "p",     output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr") %>%
    mutate( # max at 3
      stars     = ifelse(stars     == "****", "***", as.character(stars)),
      stars_fdr = ifelse(stars_fdr == "****", "***", as.character(stars_fdr))
    )
  # prepare factor levels for axes
  if (Relabel) {
    xorder <- sjlabelled::get_label(
      Data %>% select(all_of(OutcomeVars)),
      def.value = OutcomeVars
    )
    yorder <- sjlabelled::get_label(
      Data %>% select(all_of(PredictorVars)),
      def.value = PredictorVars
    )
  } else {
    xorder <- OutcomeVars
    yorder <- PredictorVars
  }
  xorder <- unique(xorder)
  yorder <- unique(yorder)

  # reorder and relabel
  P_Combined <- P_Combined %>%
    filter(!is.na(p)) %>%
    mutate(
      XLabel = factor(XLabel, levels = xorder),
      YLabel = factor(YLabel, levels = yorder)
    )

  # define significance levels so 'ns' is first/lightest
  sig_levels <- c("ns","*","**","***","****")

  # apply significance ordering and compute -log10(p)
  P_Combined <- P_Combined %>%
    mutate(
      stars     = factor(stars,     levels = sig_levels),
      stars_fdr = factor(stars_fdr, levels = sig_levels),
      logp      = -log10(p),
      logpfdr   = -log10(p_adj)
    )

  # plotting
  p <- ggplot(P_Combined, aes(y = YLabel, x = XLabel)) +
    geom_point(aes(size = stars, colour = stars), shape = 18) +
    scale_color_manual(
      values = paletteer::paletteer_c(
        "grDevices::Purple-Yellow", n = 5,
        direction = -1
      )
    ) +
    guides(size = "none") +
    labs(subtitle = "No Multiple Comparison Correction") +
    xlab("") + ylab("") + theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    )

  p_FDR <- ggplot(P_Combined, aes(y = YLabel, x = XLabel)) +
    geom_point(aes(size = stars_fdr, colour = stars_fdr), shape = 18) +
    scale_color_manual(
      values = paletteer::paletteer_c(
        "grDevices::Purple-Yellow", n = 5,
        direction = -1
      )
    ) +
    guides(size = "none") +
    labs(subtitle = "FDR Correction") +
    xlab("") + ylab("") + theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    )

  # return
  M     <- list(PvalTable = P_Combined, plot = p)
  M_FDR <- list(PvalTable = P_Combined, plot = p_FDR)
  return(list(
    Unadjusted   = M,
    FDRCorrected = M_FDR,
    method       = method,
    Relabel      = Relabel,
    Covariates   = Covariates
  ))
}
