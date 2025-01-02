#' Add a Caption Explaining Star Annotations
#'
#' This function adds a caption to a ggplot explaining the meaning of star annotations (*, **, ***).
#'
#' @return A `labs()` object that can be added to a ggplot, specifically designed for the SciDataReportR heatmaps
#' @export
geom_starcaption <- function() {
  # Create a caption with the desired text
  caption <- expression(paste("* = p < 0.05, ** = p < 0.01, *** = p < 0.001"))

  # Use labs to add the caption to the plot
  labs(caption = caption)
}
