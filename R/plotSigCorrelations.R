#' Plot Significant Correlations
#'
#' Generate scatterplots for significant correlations based on a previously generated correlation heatmap.
#'
#' @param data The dataset used to generate the scatterplots.
#' @param CorrelationHeatmapObject The output of the PlotCorrelationsHeatmap function.
#' @param PVar P-value column used to filter significance. Defaults to raw
#'   `"PValue"`; use `"PValueFDR"` for BH-adjusted selection.
#' @param Pthresh The significance threshold (default is 0.05).
#' @return A list of scatterplot objects for significant correlations.
#' @importFrom ggstatsplot ggscatterstats
#' @importFrom dplyr filter
#' @importFrom sjlabelled get_label
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @examples
#' \donttest{
#' # Build a correlation heatmap, then plot the significant pairs. Keep the
#' # predictor and outcome sets disjoint so no variable is correlated with
#' # itself.
#' ch <- PlotCorrelationsHeatmap(
#'   mtcars,
#'   predictor_vars = c("wt", "hp", "disp"),
#'   outcome_vars = c("mpg", "qsec")
#' )
#'
#' plots <- plotSigCorrelations(mtcars, ch)
#'
#' # Display the first significant-correlation scatterplot
#' plots[[1]]
#' }
#' @export
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
plotSigCorrelations <- function(data,
    CorrelationHeatmapObject,
    PVar = "PValue",
    Pthresh = 0.05,
    DataFrame = lifecycle::deprecated()) {
  .scidr_inference <- ScidrStartInference("plotSigCorrelations", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "plotSigCorrelations(DataFrame)", "plotSigCorrelations(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data


  Plot <- CorrelationHeatmapObject

  # Extract relevant data from the Plot object
  Relabeled <- Plot$Relabel
  SigVarCombos <- ScidrSignificantPairs(Plot, PVar, Pthresh, .scidr_inference, environment())

  # Initialize list to store scatterplot objects
  plist <- list()

  # Generate scatterplots for each significant correlation
  for (i in seq_len(nrow(SigVarCombos))) {
    p <- ScidrTryInference(ggscatterstats(DataFrame, x = !!SigVarCombos$XVar[i], y = !!SigVarCombos$YVar[i],
                        xsidehistogram.args = list(fill = "#6C568CFF", color = "black", na.rm = TRUE),
                        ysidehistogram.args = list(fill = "#607345FF", color = "black", na.rm = TRUE)
    ))
    if (is.null(p)) next

    # Customize subtitle with correlation information
    p$labels$subtitle <- bquote(paste(widehat(italic("r")) == .(SigVarCombos$R[i] %>% round(5)), ", ",
                                     italic("p") == .(SigVarCombos$P[i] %>% round(5)), ", " ,
                                     italic("p_FDR") == .(SigVarCombos$P_adj[i] %>% round(5)), ", ",
                                     italic("n_pairs") == .(SigVarCombos$nPairs[i])))

    # Check if adjusted for covariates and customize caption accordingly
    covars <- Plot$Covariates
    if (is.null(covars)) {
      p$labels$caption <- "Not adjusted for any covariates"
    } else {
      if (!Relabeled) {
        p$labels$caption <- paste("Adjusted for", paste(covars, collapse = ", "))
      } else {
        covarlabels <- sjlabelled::get_label(DataFrame %>% dplyr::select(dplyr::all_of(covars)), def.value = covars)
        p$labels$caption <- paste("Adjusted for", paste(covarlabels, collapse = ", "))
      }
    }

    # Change labels if Relabel is TRUE
    if (Relabeled) {
      p$labels$x <- sjlabelled::get_label(DataFrame %>% dplyr::select(dplyr::all_of(SigVarCombos$XVar[i])), def.value = SigVarCombos$XVar[i])[[1]]
      p$labels$y <- sjlabelled::get_label(DataFrame %>% dplyr::select(dplyr::all_of(SigVarCombos$YVar[i])), def.value = SigVarCombos$YVar[i])[[1]]
    }

    # Add scatterplot object to the list
    plist[[i]] <- p
  }

  # Return the list of scatterplot objects
  return(plist)
}
