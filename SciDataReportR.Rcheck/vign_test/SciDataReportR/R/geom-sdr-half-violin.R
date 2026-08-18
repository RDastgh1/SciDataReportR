# Internal half-violin geom for SciDataReportR split-violin workflows.
#
# Adapted from gghalves::geom_half_violin(), originally by Frederik
# Tiedemann. gghalves is MIT licensed:
# https://github.com/erocoar/gghalves
#
# This internal version intentionally keeps only the minimal behavior needed
# by PlotSplitViolin(): left/right half-violin polygons from ggplot2's
# y-density statistic. It is not exported as a general plotting API.

geom_sdr_half_violin <- function(mapping = NULL,
                                 data = NULL,
                                 stat = "ydensity",
                                 position = "identity",
                                 ...,
                                 side = "l",
                                 nudge = 0,
                                 trim = TRUE,
                                 scale = "area",
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSdrHalfViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      side = side,
      nudge = nudge,
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  )
}

GeomSdrHalfViolin <- ggplot2::ggproto(
  "GeomSdrHalfViolin",
  ggplot2::GeomViolin,

  setup_data = function(data, params) {
    x_data <- ggplot2::GeomBoxplot$setup_data(data, params)
    data$xmin <- x_data$xmin
    data$xmax <- x_data$xmax
    data
  },

  setup_params = function(data, params) {
    params$side <- match.arg(params$side, c("l", "r"))
    params
  },

  draw_group = function(data, panel_params, coord, side = "l", nudge = 0) {
    if (identical(side, "l")) {
      data$xminv <- data$x + data$violinwidth * (data$xmin - data$x) - nudge
      data$xmaxv <- data$x - nudge
    } else {
      data$xminv <- data$x + nudge
      data$xmaxv <- data$x + data$violinwidth * (data$xmax - data$x) + nudge
    }

    newdata <- rbind(
      transform(data, x = xminv)[order(data$y), ],
      transform(data, x = xmaxv)[order(data$y, decreasing = TRUE), ]
    )

    newdata <- rbind(newdata, newdata[1, ])

    ggplot2::GeomPolygon$draw_panel(newdata, panel_params, coord)
  }
)
