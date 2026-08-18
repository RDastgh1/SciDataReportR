#' Add values to a biomarker performance heatmap
#'
#' Adds numeric cell labels to the heatmap returned by
#' [ScreenBiomarkerPerformance()]. This helper is intended for downstream
#' annotation of SciDataReportR biomarker heatmaps while keeping the default
#' heatmap uncluttered and hover-friendly.
#'
#' @param plot A biomarker heatmap ggplot, typically
#'   `ScreenBiomarkerPerformance(...)$Plots$Heatmap`.
#' @param value_var Column in `plot$data` used for labels. Default is
#'   `"HeatmapValue"`.
#' @param digits Number of decimal places. Default is `2`.
#' @param size Text size. Default is `3`.
#' @param color Text color. Default is `"black"`.
#' @return A ggplot with numeric values added to each non-missing heatmap cell.
#'
#' @examples
#' \dontrun{
#' screen$Plots$Heatmap %>% add_biomarker_values()
#' }
#' @export
add_biomarker_values <- function(
    plot,
    value_var = "HeatmapValue",
    digits = 2,
    size = 3,
    color = "black") {

  if (!inherits(plot, "ggplot")) {
    stop("plot must be a ggplot object.", call. = FALSE)
  }

  if (!is.data.frame(plot$data) || !value_var %in% names(plot$data)) {
    stop(
      "Column '", value_var,
      "' was not found in plot$data. Use a heatmap returned by ScreenBiomarkerPerformance().",
      call. = FALSE
    )
  }

  label_data <- plot$data %>%
    dplyr::filter(!is.na(.data[[value_var]])) %>%
    dplyr::mutate(
      .scidr_biomarker_label = sprintf(
        paste0("%.", digits, "f"),
        .data[[value_var]]
      )
    )

  plot +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(label = .data$.scidr_biomarker_label),
      inherit.aes = TRUE,
      size = size,
      color = color
    )
}
