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
#' @export
geom_starcaption <- function() {
  # Create a caption with the desired text
  caption <- expression(paste("* = p < 0.05, ** = p < 0.01, *** = p < 0.001"))

  # Use labs to add the caption to the plot
  labs(caption = caption)
}
