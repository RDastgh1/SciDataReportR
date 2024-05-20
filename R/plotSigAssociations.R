#' Plot Significant Associations
#'
#' Generate ggbetweenstats for significant correlations based on a previously generated anova matrix
#'
#' @param DataFrame The dataset used to generate the scatterplots.
#' @param AnovaMatrixObject The output of the PlotAnovaRelationshipsMatrix function.
#' @param PVar The name of the column used to filter for significance (default is "P").
#' @param Pthresh The significance threshold (default is 0.05).
#' @return A list of scatterplot objects for significant correlations.
#' @importFrom ggstatsplot ggscatterstats
#' @importFrom dplyr filter
#' @importFrom sjlabelled get_label
#' @export
plotSigAssociations <- function(DataFrame, AnovaMatrixObject, PVar = "p", Pthresh = 0.05) {
  DataFrame <- as.data.frame(DataFrame)
  Plot <- AnovaMatrixObject

  # Extract relevant data from the Plot object
  Relabeled <- Plot$Relabel
  SigVarCombos <- Plot$Unadjusted$plot$data %>% as.data.frame() %>% filter(!!sym(PVar) < Pthresh)

  # Initialize list to store scatterplot objects
  plist <- list()

  # Generate scatterplots for each significant correlation
  for (i in seq(1, nrow(SigVarCombos))) {
    categorical_var <- sym(SigVarCombos$CategoricalVariable[i])
    continuous_var <- sym(as.character(SigVarCombos$ContinuousVariable[i]))

    p <- ggbetweenstats(DataFrame, x = !!categorical_var, y = !!continuous_var)

    # Customize subtitle with correlation information
    # p$labels$subtitle = bquote(paste(widehat(italic("r")) == .(SigVarCombos$R[i] %>% round(5)), ", ",
    #                                  italic("p") == .(SigVarCombos$P[i] %>% round(5)), ", " ,
    #                                  italic("p_FDR") == .(SigVarCombos$P_adj[i] %>% round(5)), ", ",
    #                                  italic("n_pairs") == .(SigVarCombos$nPairs[i])))
    #
    # Check if adjusted for covariates and customize caption accordingly
    # covars <- Plot$Covariates
    # if (is.null(covars)) {
    #   p$labels$caption = "Not adjusted for any covariates"
    # } else {
    #   if (!Relabeled) {
    #     p$labels$caption <- paste("Adjusted for", covars)
    #   } else {
    #     covarlabels <- sjlabelled::get_label(DataFrame %>% select(covars))
    #     p$labels$caption <- paste("Adjusted for", covarlabels)
    #   }
    # }
    #
    # Change labels if Relabel is TRUE
    if (Relabeled) {
      p$labels$x <- sjlabelled::get_label(DataFrame %>% select(all_of(SigVarCombos$CategoricalVariable[i])))[[1]]
      p$labels$y <- sjlabelled::get_label(DataFrame %>% select(all_of(SigVarCombos$ContinuousVariable[i])))[[1]]
    }

    # Add scatterplot object to the list
    plist[[i]] <- p
  }

  # Return the list of scatterplot objects
  return(plist)
}
