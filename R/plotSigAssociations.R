#' Plot Significant Associations
#'
#' Generate ggbetweenstats for significant correlations based on a previously generated anova matrix
#'
#' @param data The dataset used to generate the scatterplots.
#' @param AnovaMatrixObject The output of the PlotAnovaRelationshipsMatrix function.
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
#' # Build an ANOVA relationships matrix, then plot the significant pairs
#' av <- PlotAnovaRelationshipsMatrix(
#'   mtcars,
#'   CatVars = c("cyl", "gear"),
#'   ContVars = c("mpg", "wt", "hp")
#' )
#'
#' plots <- plotSigAssociations(mtcars, av)
#'
#' # Display the first significant-association plot
#' plots[[1]]
#' }
#' @export
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
plotSigAssociations <- function(data,
    AnovaMatrixObject,
    PVar = "PValue",
    Pthresh = 0.05,
    DataFrame = lifecycle::deprecated()) {
  .scidr_inference <- ScidrStartInference("plotSigAssociations", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "plotSigAssociations(DataFrame)", "plotSigAssociations(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  DataFrame <- as.data.frame(DataFrame)
  Plot <- AnovaMatrixObject

  # Extract relevant data from the Plot object
  Relabeled <- Plot$Relabel
  SigVarCombos <- ScidrSignificantPairs(Plot, PVar, Pthresh, .scidr_inference, environment())

  # Initialize list to store scatterplot objects
  plist <- list()

  # Generate scatterplots for each significant correlation
  for (i in seq_len(nrow(SigVarCombos))) {
    categorical_var <- sym(SigVarCombos$CategoricalVariable[i])
    continuous_var <- sym(as.character(SigVarCombos$ContinuousVariable[i]))

    p <- ScidrTryInference(ggbetweenstats(DataFrame, x = !!categorical_var, y = !!continuous_var))
    if (is.null(p)) next

    # Customize subtitle with correlation information
    # p$labels$subtitle <- bquote(paste(widehat(italic("r")) == .(SigVarCombos$R[i] %>% round(5)), ", ",
    #                                  italic("p") == .(SigVarCombos$P[i] %>% round(5)), ", " ,
    #                                  italic("p_FDR") == .(SigVarCombos$P_adj[i] %>% round(5)), ", ",
    #                                  italic("n_pairs") == .(SigVarCombos$nPairs[i])))
    #
    # Check if adjusted for covariates and customize caption accordingly
    # covars <- Plot$Covariates
    # if (is.null(covars)) {
    #   p$labels$caption <- "Not adjusted for any covariates"
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
      p$labels$x <- sjlabelled::get_label(DataFrame %>% dplyr::select(dplyr::all_of(SigVarCombos$CategoricalVariable[i])))[[1]]
      p$labels$y <- sjlabelled::get_label(DataFrame %>% dplyr::select(dplyr::all_of(SigVarCombos$ContinuousVariable[i])))[[1]]
    }

    # Add scatterplot object to the list
    plist[[i]] <- p
  }

  # Return the list of scatterplot objects
  return(plist)
}
