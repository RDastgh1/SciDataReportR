#' Plot Chi-Square Test with Covariate
#'
#' This function conducts chi-square tests with covariates and visualizes the results.
#'
#' @param Data The dataset containing the variables.
#' @param xVars A character vector of the names of the x-axis categorical variables.
#' @param yVars A character vector of the names of the y-axis categorical variables.
#' @param covars A character vector of the names of covariate variables. Defaults to NULL.
#' @param Relabel Logical indicating whether to relabel the variables based on their labels in the dataset. Defaults to TRUE.
#' @return A list containing ggplot objects for visualizations and tables of p-values.
#' @export
PlotChiSqCovar <- function(Data, xVars, yVars, covars = NULL, Relabel = TRUE, Ordinal = TRUE) {
  if(is.null(yVars)){
    yVars = xVars
  }

  DataSubset <- Data[unique(c(xVars, yVars, covars))]
  DataSubset[unique(c(xVars, yVars))] <- lapply(DataSubset[unique(c(xVars, yVars))], factor)
  varOverlap <- yVars[yVars %in% xVars]
  DataSubset$rowID <- seq(1, nrow(DataSubset))

  # Remove NAs

  # Measure of association

  mData <- pivot_longer(DataSubset, cols = all_of(xVars))
  if(length(varOverlap) > 0){
    mData <- mData %>% left_join(DataSubset, by = "rowID")
  }
  mData <- pivot_longer(mData, cols = all_of(yVars), names_to = "Covariate", values_to = "Category")

  cOverlapV <- covars[covars %in% unique(c(xVars, yVars))]
  # Remove NAs

  if(length(cOverlapV) > 0){
    mData <- left_join(mData, DataSubset[c("rowID", cOverlapV)], by = "rowID")
  }

  mData <- na.omit(mData)

  nlevels.test <- mData %>% group_by(name, Covariate) %>%
    mutate(n = length(unique((Category))))

  mData <- mData[nlevels.test$n > 1,]

  nlevels.test <- mData %>% group_by(name, Covariate) %>%
    mutate(n = length(unique((value))))

  mData <- mData[nlevels.test$n > 1,]

  stat.test <- mData %>% group_by(name, Covariate) %>%
    summarise(pval = chisq.test(value, Category)$p.value) %>% add_significance() %>%
    adjust_pvalue(method = "fdr") %>% add_significance()

  stat.test$logp <- -log10(stat.test$pval)
  stat.test$`p<.05` <- stat.test$pval.signif
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns", "*", "**", "***", "****"))

  stat.test$logp_FDR <- -log10(stat.test$pval.adj)

  stat.test$pval.adj.signif <- factor(stat.test$pval.adj.signif, levels = c("ns", "*", "**", "***", "****"))

  # If there is a diagonal, remove it
  stat.test[stat.test$name == stat.test$Covariate, 3:ncol(stat.test)] <- NA

  if(Relabel == TRUE){
    Data<- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[stat.test$Covariate], def.value = colnames(Data[stat.test$Covariate])) %>% as.data.frame() %>% rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[stat.test$name], def.value = colnames(Data[stat.test$name])) %>% as.data.frame() %>% rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label

  } else {
    stat.test$XLabel <- stat.test$Covariate
    stat.test$YLabel <- stat.test$name
  }

  PlotText <- paste("</br>XVar Label:", stat.test$XLabel,
                    "</br>XVar:", stat.test$Covariate,
                    "</br>YVar Label:", stat.test$YLabel,
                    "</br>YVar:", stat.test$name,
                    "</br>P-Value: ", stat.test$pval, stat.test$pval.signif,
                    "</br>FDR-corrected P: ", stat.test$pval.adj, stat.test$pval.adj.signif)

  if(Relabel == TRUE){
    # Convert to factor and set labels
    stat.test$name <- factor(stat.test$name, levels = xVars, labels = sjlabelled::get_label(Data[xVars], def.value = colnames(Data[xVars])))
    stat.test$Covariate <- factor(stat.test$Covariate, levels = yVars, labels = sjlabelled::get_label(Data[yVars], def.value = colnames(Data[yVars])))
  } else {
    stat.test$name <- factor(stat.test$name, levels = xVars)
    stat.test$Covariate <- factor(stat.test$Covariate, levels = yVars)
  }

  stat.test$logp_lim <- if_else(stat.test$logp > 5, 5, stat.test$logp)
  stat.test$logp_FDR_lim <- if_else(stat.test$logp_FDR > 5, 5, stat.test$logp_FDR)

  p_g <- ggplot(stat.test, aes(y = name, x = Covariate, shape = `p<.05`, text = PlotText)) +
    geom_point(aes(size = logp_lim, colour = pval)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_shape_manual(values = c(7, 16, 17, 15, 18), breaks = c("ns", "*", "**", "***", "****"), drop = FALSE) +
    scale_color_gradientn(trans = "log", colours = rev(brewer.pal(9, "Oranges")[-1]), limits = c(0.0000001, 1), oob = scales::squish, breaks = c(1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)) +
    scale_size_continuous(breaks = seq(0, 10, 0.5), limits = c(0, 5)) + guides(size = FALSE) +
    labs(subtitle = "No Multiple Comparison Correction")

  p_g_FDR <- ggplot(stat.test, aes(y = name, x = Covariate, shape = pval.adj.signif, text = PlotText)) +
    geom_point(aes(size = logp_FDR_lim, colour = pval.adj)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_shape_manual(values = c(7, 16, 17, 15, 18), breaks = c("ns", "*", "**", "***", "****"), drop = FALSE) +
    scale_color_gradientn(trans = "log", colours = rev(brewer.pal(9, "Oranges")), limits = c(0.0000001, 1), oob = scales::squish, breaks = c(1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)) +
    scale_size_continuous(breaks = seq(0, 10, 0.5), limits = c(0, 5)) + guides(size = FALSE) +
    labs(subtitle = "FDR Correction")

  # Convert stat.test to table

  pvaltable <- stat.test %>% data.frame() %>% select(c(name, Covariate, pval)) %>% pivot_wider(names_from = name, values_from = pval) %>% filter(complete.cases(Covariate))

  pvaltable_FDR <- stat.test %>% data.frame() %>% select(c(name, Covariate, pval.adj)) %>% pivot_wider(names_from = name, values_from = pval.adj) %>% filter(complete.cases(Covariate))

  return(list(p = p_g, pvaltable = pvaltable, p_FDR = p_g_FDR, pvaltable_FDR = pvaltable_FDR))
}
