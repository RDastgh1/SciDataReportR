#' Plot Significant Correlations
#'
#' Generate scatterplots for significant correlations based on a previously generated correlation heatmap.
#'
#' @param data The dataset used to generate the scatterplots.
#' @param CorrelationHeatmapObject The output of the PlotCorrelationsHeatmap function.
#' @param PVar The name of the column used to filter for significance (default is "P").
#' @param Pthresh The significance threshold (default is 0.05).
#' @return A list of scatterplot objects for significant correlations.
#' @importFrom ggstatsplot ggscatterstats
#' @importFrom dplyr filter
#' @importFrom sjlabelled get_label
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @examples
#' \donttest{
#' # Build a correlation heatmap, then plot the significant pairs
#' ch <- PlotCorrelationsHeatmap(
#'   mtcars,
#'   predictor_vars = c("mpg", "wt", "hp"),
#'   outcome_vars = c("mpg", "wt", "hp")
#' )
#'
#' plots <- plotSigCorrelations(mtcars, ch)
#'
#' # Display the first significant-correlation scatterplot
#' plots[[1]]
#' }
#' @export
plotSigCorrelations <- function(data,
    CorrelationHeatmapObject,
    PVar = "P",
    Pthresh = 0.05,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "plotSigCorrelations(DataFrame)", "plotSigCorrelations(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data


  Plot <- CorrelationHeatmapObject

  # Extract relevant data from the Plot object
  Relabeled <- Plot$Relabel
  SigVarCombos <- Plot$Unadjusted$plot$data %>% filter(!!sym(PVar) < Pthresh)

  # Initialize list to store scatterplot objects
  plist <- list()

  # Generate scatterplots for each significant correlation
  for (i in seq(1, nrow(SigVarCombos))) {
    p <- ggscatterstats(DataFrame, x = !!SigVarCombos$XVar[i], y = !!SigVarCombos$YVar[i],
                        xsidehistogram.args = list(fill = "#6C568CFF", color = "black", na.rm = TRUE),
                        ysidehistogram.args = list(fill = "#607345FF", color = "black", na.rm = TRUE)
    )

    # Customize subtitle with correlation information
    p$labels$subtitle = bquote(paste(widehat(italic("r")) == .(SigVarCombos$R[i] %>% round(5)), ", ",
                                     italic("p") == .(SigVarCombos$P[i] %>% round(5)), ", " ,
                                     italic("p_FDR") == .(SigVarCombos$P_adj[i] %>% round(5)), ", ",
                                     italic("n_pairs") == .(SigVarCombos$nPairs[i])))

    # Check if adjusted for covariates and customize caption accordingly
    covars <- Plot$Covariates
    if (is.null(covars)) {
      p$labels$caption = "Not adjusted for any covariates"
    } else {
      if (!Relabeled) {
        p$labels$caption <- paste("Adjusted for", paste(covars, collapse = ", "))
      } else {
        covarlabels <- sjlabelled::get_label(DataFrame %>% select(all_of(covars)), def.value = covars)
        p$labels$caption <- paste("Adjusted for", paste(covarlabels, collapse = ", "))
      }
    }

    # Change labels if Relabel is TRUE
    if (Relabeled) {
      p$labels$x <- sjlabelled::get_label(DataFrame %>% select(all_of(SigVarCombos$XVar[i])), def.value = SigVarCombos$XVar[i])[[1]]
      p$labels$y <- sjlabelled::get_label(DataFrame %>% select(all_of(SigVarCombos$YVar[i])), def.value = SigVarCombos$YVar[i])[[1]]
    }

    # Add scatterplot object to the list
    plist[[i]] <- p
  }

  # Return the list of scatterplot objects
  return(plist)
}
