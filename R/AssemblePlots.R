#' Assemble ggplot objects into a unified multi-panel figure
#'
#' Combine ggplot objects or lists of ggplot objects into a consistently styled
#' multi-panel layout with optional legend collection, global theming,
#' automatic layout estimation, and suggested figure dimensions for reporting
#' workflows.
#'
#' This function is designed as a centralized visualization orchestration layer
#' for SciDataReportR workflows. It standardizes layout behavior across
#' analytical outputs while preserving full ggplot compatibility.
#'
#' Internally, cowplot is used for layout assembly and legend extraction
#' because it provides more predictable spacing and legend behavior for
#' publication-style figures than patchwork.
#'
#' Supports:
#' - Named or unnamed plot lists
#' - Automatic row/column estimation
#' - Global themes
#' - Global ggplot layers/scales
#' - Shared legends
#' - Automatic NULL removal
#' - Suggested figure dimensions
#'
#' @param Plots A ggplot object or list of ggplot objects.
#' @param ncol Optional number of columns.
#' @param nrow Optional number of rows.
#' @param AutoLayout Logical; automatically estimate layout when
#'   `ncol` and `nrow` are not supplied.
#' @param RemoveNULL Logical; remove NULL plots before assembly.
#' @param CollectLegend Logical; combine legends into a shared legend.
#' @param LegendPosition Position of shared legend. One of
#'   `"top"`, `"bottom"`, `"left"`, `"right"`, or `"none"`.
#' @param LegendRelativeSize Relative size allocated to legend area.
#' @param Theme Optional global ggplot theme.
#' @param BaseFontSize Base font size applied globally.
#' @param GlobalTheme Optional additional theme applied globally.
#' @param GlobalLayers Optional list of ggplot layers/scales applied
#'   to all plots.
#' @param RemoveTitles Logical; remove plot titles globally.
#' @param UseNamesAsTitles Logical; use plot list names as titles when
#'   plot titles are missing.
#' @param Align Plot alignment passed to cowplot.
#' @param Axis Axis alignment passed to cowplot.
#' @param Labels Optional panel labels.
#' @param LabelSize Panel label font size.
#' @param SuggestedBaseWidth Base width per column used for suggested
#'   figure dimensions.
#' @param SuggestedBaseHeight Base height per row used for suggested
#'   figure dimensions.
#' @param ReturnMetadata Logical; if `TRUE`, returns plot plus metadata.
#'
#' @return
#' If `ReturnMetadata = FALSE`, returns a ggplot object.
#'
#' If `ReturnMetadata = TRUE`, returns a list containing:
#' \describe{
#'   \item{Plot}{Combined plot object}
#'   \item{nrow}{Estimated number of rows}
#'   \item{ncol}{Estimated number of columns}
#'   \item{SuggestedWidth}{Suggested figure width}
#'   \item{SuggestedHeight}{Suggested figure height}
#'   \item{NumPlots}{Number of plots}
#' }
#'
#' @examples
#' library(ggplot2)
#'
#' p1 <- ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point()
#'
#' p2 <- ggplot(mtcars, aes(hp, wt)) +
#'   geom_point()
#'
#' AssemblePlots(list(p1, p2))
#'
#' AssemblePlots(
#'   list(MPG = p1, Horsepower = p2),
#'   UseNamesAsTitles = TRUE,
#'   LegendPosition = "top"
#' )
#'
#' @export
AssemblePlots <- function(
    Plots,
    ncol = NULL,
    nrow = NULL,
    AutoLayout = TRUE,
    RemoveNULL = TRUE,
    CollectLegend = TRUE,
    LegendPosition = "bottom",
    LegendRelativeSize = 0.08,
    Theme = ggplot2::theme_minimal(),
    BaseFontSize = 12,
    GlobalTheme = NULL,
    GlobalLayers = NULL,
    RemoveTitles = FALSE,
    UseNamesAsTitles = FALSE,
    Align = "hv",
    Axis = "tblr",
    Labels = NULL,
    LabelSize = 14,
    SuggestedBaseWidth = 4,
    SuggestedBaseHeight = 4,
    ReturnMetadata = FALSE
) {

  # Validate inputs

  if (inherits(Plots, "ggplot")) {
    Plots <- list(Plots)
  }

  if (!is.list(Plots)) {
    stop("Plots must be a ggplot object or list of ggplot objects.")
  }

  valid_positions <- c(
    "top",
    "bottom",
    "left",
    "right",
    "none"
  )

  if (!LegendPosition %in% valid_positions) {
    stop(
      "LegendPosition must be one of: ",
      paste(valid_positions, collapse = ", ")
    )
  }

  # Flatten nested lists

  flatten_plots <- function(x) {
    unlist(x, recursive = FALSE)
  }

  Plots <- flatten_plots(Plots)

  # Remove NULL plots

  if (RemoveNULL) {

    Plots <- Plots[
      !vapply(Plots, is.null, logical(1))
    ]
  }

  if (length(Plots) == 0) {
    stop("No valid plots found.")
  }

  # Validate ggplot objects

  valid_plots <- vapply(
    Plots,
    function(x) inherits(x, "ggplot"),
    logical(1)
  )

  if (!all(valid_plots)) {
    stop("All elements in Plots must be ggplot objects.")
  }

  # Infer layout

  nplots <- length(Plots)

  if (AutoLayout &&
      is.null(ncol) &&
      is.null(nrow)) {

    if (nplots <= 2) {

      ncol <- nplots

    } else if (nplots <= 4) {

      ncol <- 2

    } else if (nplots <= 9) {

      ncol <- 3

    } else {

      ncol <- ceiling(sqrt(nplots))
    }

    nrow <- ceiling(nplots / ncol)
  }

  if (!is.null(ncol) && is.null(nrow)) {
    nrow <- ceiling(nplots / ncol)
  }

  if (!is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(nplots / nrow)
  }

  # Apply names as titles

  if (UseNamesAsTitles &&
      !is.null(names(Plots))) {

    for (i in seq_along(Plots)) {

      current_title <- Plots[[i]]$labels$title

      if (is.null(current_title) ||
          identical(current_title, "")) {

        if (!is.null(names(Plots)[i]) &&
            names(Plots)[i] != "") {

          Plots[[i]] <- Plots[[i]] +
            ggplot2::ggtitle(names(Plots)[i])
        }
      }
    }
  }

  # Remove titles

  if (RemoveTitles) {

    Plots <- lapply(
      Plots,
      function(p) {

        p +
          ggplot2::theme(
            plot.title = ggplot2::element_blank()
          )
      }
    )
  }

  # Apply base theme

  if (!is.null(Theme)) {

    base_theme <- Theme +
      ggplot2::theme(
        text = ggplot2::element_text(
          size = BaseFontSize
        )
      )

    Plots <- lapply(
      Plots,
      function(p) {
        p + base_theme
      }
    )
  }

  # Apply additional global theme

  if (!is.null(GlobalTheme)) {

    Plots <- lapply(
      Plots,
      function(p) {
        p + GlobalTheme
      }
    )
  }

  # Apply global layers/scales

  if (!is.null(GlobalLayers)) {

    for (layer in GlobalLayers) {

      Plots <- lapply(
        Plots,
        function(p) {
          p + layer
        }
      )
    }
  }

  # Extract shared legend

  legend <- NULL

  if (CollectLegend &&
      LegendPosition != "none") {

    legend_plot <- Plots[[1]] +
      ggplot2::theme(
        legend.position = LegendPosition,
        legend.box = ifelse(
          LegendPosition %in% c("top", "bottom"),
          "horizontal",
          "vertical"
        )
      )

    legend <- cowplot::get_legend(legend_plot)

    # Remove legends from individual plots

    Plots <- lapply(
      Plots,
      function(p) {
        p +
          ggplot2::theme(
            legend.position = "none"
          )
      }
    )
  }

  # Build plot grid

  plot_grid <- cowplot::plot_grid(
    plotlist = Plots,
    ncol = ncol,
    nrow = nrow,
    align = Align,
    axis = Axis,
    labels = Labels,
    label_size = LabelSize
  )

  # Combine legend and plots

  final_plot <- plot_grid

  if (!is.null(legend)) {

    if (LegendPosition == "top") {

      final_plot <- cowplot::plot_grid(
        legend,
        plot_grid,
        ncol = 1,
        rel_heights = c(
          LegendRelativeSize,
          1
        )
      )

    } else if (LegendPosition == "bottom") {

      final_plot <- cowplot::plot_grid(
        plot_grid,
        legend,
        ncol = 1,
        rel_heights = c(
          1,
          LegendRelativeSize
        )
      )

    } else if (LegendPosition == "left") {

      final_plot <- cowplot::plot_grid(
        legend,
        plot_grid,
        nrow = 1,
        rel_widths = c(
          LegendRelativeSize,
          1
        )
      )

    } else if (LegendPosition == "right") {

      final_plot <- cowplot::plot_grid(
        plot_grid,
        legend,
        nrow = 1,
        rel_widths = c(
          1,
          LegendRelativeSize
        )
      )
    }
  }

  # Estimate dimensions

  SuggestedWidth <- ncol * SuggestedBaseWidth
  SuggestedHeight <- nrow * SuggestedBaseHeight

  if (LegendPosition %in% c("top", "bottom")) {
    SuggestedHeight <- SuggestedHeight * 1.08
  }

  if (LegendPosition %in% c("left", "right")) {
    SuggestedWidth <- SuggestedWidth * 1.08
  }

  # Return result

  if (ReturnMetadata) {

    return(list(
      Plot = final_plot,
      nrow = nrow,
      ncol = ncol,
      SuggestedWidth = SuggestedWidth,
      SuggestedHeight = SuggestedHeight,
      NumPlots = nplots
    ))

  } else {

    return(final_plot)
  }
}
