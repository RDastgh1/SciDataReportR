#' SciDataReportR qualitative color palette
#'
#' Return colors from the SciDataReportR qualitative color system. The palette
#' is designed for categorical scientific graphics and begins with the package
#' anchor colors, Navy and Orange.
#'
#' @param n Number of colors to return. If `NULL`, return the full palette.
#' @param names Logical indicating whether to preserve the color names.
#'
#' @return A character vector of hexadecimal color values.
#'
#' @examples
#' SciDataPalette(3)
#' SciDataPalette(8)
#' SciDataPalette()
#'
#' @export
SciDataPalette <- function(n = NULL, names = TRUE) {

  colors <- c(
    Navy          = "#0B1F5E",
    Orange        = "#E87400",
    Emerald       = "#16835F",
    Cranberry     = "#C2185B",
    Plum          = "#641B68",
    DeepTeal      = "#0C5F68",
    Burgundy      = "#841B37",
    ForestGreen   = "#285F3D",
    Mulberry      = "#9B315F",
    Aubergine     = "#49305C",
    Espresso      = "#593B2B",
    Blue          = "#1769D2",
    Gold          = "#F2B134",
    Violet        = "#7758A3",
    Coral         = "#F2645A",
    JewelBlue     = "#4D73B5",
    LeafGreen     = "#72B65D",
    Magenta       = "#D62976",
    Lime          = "#9EBB1F",
    BurntOrange   = "#E65D0A",
    Ochre         = "#B77B22",
    Mauve         = "#B07A9E",
    Rust          = "#B44D2E",
    SteelBlue     = "#6E8FA8",
    Olive         = "#687548",
    Sienna        = "#A85B32",
    Mint          = "#8FD8B8",
    Rose          = "#D9A0AE",
    PowderBlue    = "#B8D8E8",
    Periwinkle    = "#A6AFE8",
    Lavender      = "#C7B4D9",
    Sage          = "#B5C7A3",
    Peach         = "#F2B29A",
    Cashmere      = "#E7C6B2"
  )

  if (is.null(n)) {
    out <- colors
  } else {

    if (
      !is.numeric(n) ||
      length(n) != 1 ||
      is.na(n) ||
      n < 1 ||
      n %% 1 != 0
    ) {
      stop("`n` must be a single positive integer.")
    }

    if (n > length(colors)) {
      stop(
        "`n` cannot exceed ",
        length(colors),
        ", the number of colors in the SciDataReportR palette."
      )
    }

    out <- colors[seq_len(n)]
  }

  if (!names) {
    out <- unname(out)
  }

  out
}

#' Shade a set of colors toward black or white
#'
#' Used to extend the qualitative palette past its natural length. Shading in
#' HSV keeps each derived color recognisably related to the base color it came
#' from, which a linear ramp through the palette would not.
#'
#' @param colors Character vector of colors.
#' @param amount Numeric. Positive lightens, negative darkens, roughly as a
#'   proportion of the remaining headroom.
#' @return A character vector of hexadecimal colors, the same length as
#'   `colors`.
#' @noRd
.SciDataShade <- function(colors, amount) {

  hsv_values <- grDevices::rgb2hsv(grDevices::col2rgb(colors))

  hue <- hsv_values["h", ]
  saturation <- hsv_values["s", ]
  value <- hsv_values["v", ]

  if (amount >= 0) {
    value <- value + (1 - value) * amount
    # Desaturating alongside lightening keeps pastels from all collapsing
    # toward the same washed-out tone.
    saturation <- saturation * (1 - amount * 0.6)
  } else {
    value <- value * (1 + amount)
  }

  grDevices::hsv(
    h = hue,
    s = pmin(pmax(saturation, 0), 1),
    v = pmin(pmax(value, 0), 1)
  )
}

#' Palette values for an arbitrary number of categories
#'
#' The internal counterpart to [SciDataPalette()]. Where the exported function
#' errors past the length of the palette - correct behavior for an explicit
#' request - this one always returns exactly `n` colors, so a package plot can
#' never fail merely because a variable has many levels.
#'
#' The first 34 values are the palette verbatim. Beyond that, shaded cycles of
#' the full palette are appended. Values stay distinct to 170; past that they
#' recycle with a warning, on the basis that a figure with 170 categories has
#' already stopped communicating.
#'
#' @param n Number of colors required.
#' @return A character vector of `n` colors.
#' @noRd
.SciDataColorValues <- function(n) {

  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n %% 1 != 0) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }

  base_colors <- SciDataPalette(names = FALSE)

  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  }

  out <- base_colors

  for (amount in c(-0.35, 0.45, -0.6, 0.7)) {
    if (length(out) >= n) break
    out <- c(out, .SciDataShade(base_colors, amount))
  }

  if (length(out) < n) {
    warning(
      "More than ", length(out), " categories requested; palette colors are ",
      "being reused. Consider collapsing categories before plotting.",
      call. = FALSE
    )
    out <- rep_len(out, n)
  }

  out[seq_len(n)]
}

#' Assign palette colors to levels, holding some levels at fixed colors
#'
#' Levels named in `fixed` keep their color and consume no palette slot, so a
#' reserved level such as `Noise` or a missing-data category never shifts the
#' colors assigned to the real categories.
#'
#' @param levels Character vector of factor levels, in display order.
#' @param fixed Optional named character vector of level-to-color overrides.
#' @return A named character vector of colors, ordered to match `levels`.
#' @noRd
.SciDataNamedValues <- function(levels, fixed = NULL) {

  levels <- as.character(levels)
  reserved <- intersect(levels, names(fixed))
  free <- setdiff(levels, reserved)

  values <- character(0)

  if (length(free)) {
    values <- stats::setNames(.SciDataColorValues(length(free)), free)
  }

  if (length(reserved)) {
    values <- c(values, fixed[reserved])
  }

  values[levels]
}

# Manual scales carrying enough values that a package plot never hits the
# "Insufficient values in manual scale" error. 170 is the point past which
# .SciDataColorValues() starts recycling, so this is the largest count that
# stays warning-free.
.SciDataScaleValues <- function() .SciDataColorValues(170)

#' Default SciDataReportR discrete scales for package figures
#'
#' Internal equivalents of [scale_fill_SciData()] / [scale_color_SciData()],
#' sized so that high-cardinality categories degrade rather than error.
#'
#' @param ... Passed to the underlying ggplot2 scale.
#' @return A ggplot2 discrete scale.
#' @noRd
.SciDataFillScale <- function(...) {
  ggplot2::scale_fill_manual(values = .SciDataScaleValues(), ...)
}

#' @noRd
.SciDataColourScale <- function(...) {
  ggplot2::scale_colour_manual(values = .SciDataScaleValues(), ...)
}

#' Pick a legible text color for a given background
#'
#' Returns black or white per background color using the WCAG relative
#' luminance formula. Needed wherever a label is drawn on top of a filled
#' shape: the palette spans 1.4:1 to 14:1 against black, so a fixed text color
#' is illegible on part of it.
#'
#' @param fills Character vector of background colors.
#' @param dark,light Colors returned for light and dark backgrounds.
#' @return A character vector the same length as `fills`.
#' @noRd
.SciDataContrastText <- function(fills, dark = "black", light = "white") {

  channels <- grDevices::col2rgb(fills) / 255

  linear <- ifelse(
    channels <= 0.03928,
    channels / 12.92,
    ((channels + 0.055) / 1.055)^2.4
  )

  luminance <- 0.2126 * linear[1, ] + 0.7152 * linear[2, ] + 0.0722 * linear[3, ]

  ifelse(luminance > 0.179, dark, light)
}

#' SciDataReportR discrete color scale
#'
#' Apply the SciDataReportR qualitative palette to the `color` aesthetic in a
#' ggplot. Use this scale for categorical groups represented by points, lines,
#' outlines, or other color-mapped geometries.
#'
#' @param ... Additional arguments passed to `ggplot2::scale_color_manual()`.
#'
#' @return A ggplot2 discrete color scale.
#'
#' @examples
#' ggplot2::ggplot(
#'   iris,
#'   ggplot2::aes(
#'     x = Sepal.Length,
#'     y = Sepal.Width,
#'     color = Species
#'   )
#' ) +
#'   ggplot2::geom_point() +
#'   scale_color_SciData()
#'
#' @export
scale_color_SciData <- function(...) {

  ggplot2::scale_color_manual(
    values = SciDataPalette(names = FALSE),
    ...
  )
}

#' SciDataReportR discrete fill scale
#'
#' Apply the SciDataReportR qualitative palette to the `fill` aesthetic in a
#' ggplot. Use this scale for categorical groups represented by bars, boxes,
#' violins, areas, or other filled geometries.
#'
#' @param ... Additional arguments passed to `ggplot2::scale_fill_manual()`.
#'
#' @return A ggplot2 discrete fill scale.
#'
#' @examples
#' ggplot2::ggplot(
#'   iris,
#'   ggplot2::aes(
#'     x = Species,
#'     y = Sepal.Length,
#'     fill = Species
#'   )
#' ) +
#'   ggplot2::geom_boxplot() +
#'   scale_fill_SciData()
#'
#' @export
scale_fill_SciData <- function(...) {

  ggplot2::scale_fill_manual(
    values = SciDataPalette(names = FALSE),
    ...
  )
}
