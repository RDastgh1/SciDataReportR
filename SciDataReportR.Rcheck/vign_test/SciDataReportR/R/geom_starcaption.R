#' Add a Caption Explaining Star Annotations
#'
#' This function adds a caption to a ggplot explaining the meaning of star
#' annotations (*, **, ***). It is most commonly added to correlation heatmaps
#' produced by [PlotCorrelationsHeatmap()] or to downstream plots derived from
#' that function with [add_r_and_stars()].
#'
#' @return A `labs()` object that can be added to a ggplot, especially
#'   SciDataReportR heatmaps that use star annotations.
#'
#' @section Input requirements:
#' `geom_starcaption()` does not take an input object directly. Add it to a
#' ggplot with `+`, typically a plot returned inside the list created by
#' [PlotCorrelationsHeatmap()] or a plot returned by [add_r_and_stars()].
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Add the caption to a correlation heatmap whose tiles already show stars.
#' heatmap <- PlotCorrelationsHeatmap(
#'   data = df_Labelled,
#'   predictor_vars = c("Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10"),
#'   outcome_vars = c("MMP7", "TRAIL_R3", "Ferritin", "Fibrinogen", "MIF")
#' )
#' heatmap$Unadjusted$plot + geom_starcaption()
#' @export
geom_starcaption <- function() {
  # Caption text comes from the shared star scale so it cannot drift from the
  # thresholds the plotting code actually applies.
  caption <- ScidrStarCaptionText()

  # Use labs to add the caption to the plot
  labs(caption = caption)
}
