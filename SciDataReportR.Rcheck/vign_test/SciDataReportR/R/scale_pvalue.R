# Internal constants define the evidence-aware mapping in one place so the
# scale can be refined later without changing the public API.
.pvalue_warp_breaks <- c(1, 0.10, 0.05, 0.01, 0.001, 1e-5, 1e-8)
.pvalue_warp_positions <- c(0, 0.16, 0.32, 0.56, 0.75, 0.90, 1)


# Convert raw p-values to threshold-aware positions between zero and one.
.warp_pvalues <- function(p_values) {
  evidence <- -log10(p_values)
  evidence_breaks <- -log10(.pvalue_warp_breaks)

  stats::approx(
    x = evidence_breaks,
    y = .pvalue_warp_positions,
    xout = evidence,
    rule = 2,
    ties = "ordered"
  )$y
}


# Convert threshold-aware positions back to raw p-values. The inverse is
# required so ggplot2 can construct a guide in the same coordinate system as
# the palette while continuing to display raw p-value labels.
.unwarp_pvalues <- function(positions) {
  evidence_breaks <- -log10(.pvalue_warp_breaks)

  evidence <- stats::approx(
    x = .pvalue_warp_positions,
    y = evidence_breaks,
    xout = positions,
    rule = 2,
    ties = "ordered"
  )$y

  10^(-evidence)
}


# Use a genuine scale transformation so both the palette and colorbar use the
# threshold-aware coordinates. A rescaler alone only warps palette lookup and
# leaves guide ticks compressed in raw p-value space.
.pvalue_transform <- function() {
  scales::new_transform(
    name = "pvalue_warp",
    transform = .warp_pvalues,
    inverse = .unwarp_pvalues,
    breaks = function(limits) .pvalue_warp_breaks,
    format = .format_pvalue_labels,
    domain = c(0, 1)
  )
}


# Reduce saturation only above p = 0.05 while preserving the underlying hue
# and restoring full saturation at the significance threshold.
.desaturate_pvalue_colors <- function(colors, positions) {
  threshold_position <- .warp_pvalues(0.05)

  saturation_multiplier <- rep(1, length(positions))
  nonsignificant <- !is.na(positions) & positions < threshold_position

  saturation_multiplier[nonsignificant] <-
    0.20 +
    0.80 * (
      positions[nonsignificant] / threshold_position
    )

  rgb_values <- grDevices::col2rgb(
    colors,
    alpha = TRUE
  )

  hsv_values <- grDevices::rgb2hsv(
    rgb_values[1:3, , drop = FALSE],
    maxColorValue = 255
  )

  hsv_values["s", ] <-
    hsv_values["s", ] * saturation_multiplier

  grDevices::hsv(
    h = hsv_values["h", ],
    s = hsv_values["s", ],
    v = hsv_values["v", ],
    alpha = rgb_values["alpha", ] / 255
  )
}


# Build the palette shared by the fill and color scales.
.pvalue_palette <- function(palette, direction) {
  palette <- match.arg(
    palette,
    choices = c("inferno", "viridis")
  )

  base_colors <- switch(
    palette,
    inferno = viridisLite::inferno(
      n = 256,
      direction = direction
    ),
    viridis = viridisLite::viridis(
      n = 256,
      direction = direction
    )
  )

  color_function <- grDevices::colorRampPalette(
    colors = base_colors,
    space = "Lab"
  )

  function(x) {
    missing_positions <- is.na(x)

    x <- pmin(
      pmax(x, 0),
      1
    )

    color_index <- round(x * 255) + 1
    colors <- color_function(256)[color_index]

    colors <- .desaturate_pvalue_colors(
      colors = colors,
      positions = x
    )

    colors[missing_positions] <- NA_character_
    colors
  }
}


# Format raw p-values rather than transformed evidence values in the legend.
.format_pvalue_labels <- function(x) {
  vapply(
    x,
    FUN.VALUE = character(1),
    function(value) {
      if (is.na(value)) {
        return(NA_character_)
      }

      if (value == 0) {
        return("0")
      }

      if (value >= 0.001) {
        return(
          format(
            value,
            scientific = FALSE,
            trim = TRUE,
            digits = 3
          )
        )
      }

      format(
        value,
        scientific = TRUE,
        trim = TRUE,
        digits = 1
      )
    }
  )
}


# Shared construction prevents color and fill scales from drifting apart.
.scale_pvalue <- function(
    aesthetics,
    palette,
    direction,
    name,
    breaks,
    labels,
    limits,
    na.value,
    guide,
    ...) {

  palette <- match.arg(
    palette,
    choices = c("inferno", "viridis")
  )

  if (
    !is.numeric(direction) ||
    length(direction) != 1 ||
    is.na(direction) ||
    !direction %in% c(-1, 1)
  ) {
    stop(
      "`direction` must be either -1 or 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(limits) ||
    length(limits) != 2 ||
    anyNA(limits) ||
    any(!is.finite(limits)) ||
    any(limits <= 0) ||
    limits[1] >= limits[2]
  ) {
    stop(
      paste0(
        "`limits` must be a numeric vector containing two finite, ",
        "positive values in increasing order. Received: ",
        paste(limits, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (is.null(guide)) {
    guide <- ggplot2::guide_colourbar(
      theme = ggplot2::theme(
        legend.key.height = grid::unit(70, "mm")
      )
    )
  }

  ggplot2::continuous_scale(
    aesthetics = aesthetics,
    palette = .pvalue_palette(
      palette = palette,
      direction = direction
    ),
    name = name,
    breaks = breaks,
    labels = labels,
    limits = limits,
    transform = .pvalue_transform(),
    oob = scales::squish,
    na.value = na.value,
    guide = guide,
    ...
  )
}


#' Apply an evidence-aware p-value color scale
#'
#' Maps raw p-values to either the Inferno or Viridis color palette using a
#' threshold-aware transformation. The transformation allocates additional
#' visual resolution around commonly interpreted p-value thresholds of 0.05,
#' 0.01, and 0.001. Values above 0.05 are progressively desaturated to reduce
#' their visual emphasis while retaining a continuous representation of the
#' underlying p-values. The colorbar uses the same warped coordinates, so its
#' raw p-value labels remain separated around these thresholds. Use this scale
#' when p-values are mapped to the `color` aesthetic.
#'
#' @param palette A character value specifying the color palette. Must be
#'   `"inferno"` or `"viridis"`. Defaults to `"inferno"`.
#' @param direction A numeric value controlling the direction of the palette.
#'   Use `-1`, the default, for darker colors to indicate smaller p-values and
#'   stronger statistical evidence. Use `1` to reverse the palette.
#' @param name The legend title. Defaults to `"P-value"`.
#' @param breaks A numeric vector of raw p-values to display as legend breaks.
#'   Defaults to commonly interpreted statistical thresholds and reference
#'   values.
#' @param labels A function or character vector used to label the legend
#'   breaks. By default, the legend displays raw p-values rather than
#'   transformed values.
#' @param limits A numeric vector containing the minimum and maximum p-values
#'   represented by the scale. Values outside these limits are squished to
#'   the nearest limit. Defaults to `c(1e-8, 1)`.
#' @param na.value The color assigned to missing p-values. Defaults to
#'   `"grey80"`.
#' @param guide A guide function or guide name. The default, `NULL`, uses a
#'   taller colorbar so the threshold-aware breaks remain legible. Supply
#'   `"colourbar"` for ggplot2's standard-sized colorbar or another guide to
#'   override this behavior.
#' @param ... Additional arguments passed to [ggplot2::continuous_scale()].
#'
#' @return A continuous ggplot2 color scale.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Four outcomes by six predictors: 24 standardized associations.
#' screen <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = c("tau", "p_tau", "AXL", "Ferritin"),
#'   predictor_vars = c(
#'     "age", "sex", "Adiponectin", "Cortisol", "Insulin", "Leptin"
#'   ),
#'   Standardize = TRUE
#' )
#'
#' # A volcano plot is the natural home for this scale: effect size on x,
#' # evidence on y, and the colour reinforcing where the thresholds fall so
#' # the eye is not left to judge distance along a log axis.
#' ggplot2::ggplot(
#'   screen$Results,
#'   ggplot2::aes(
#'     x = Estimate,
#'     y = -log10(PValue),
#'     color = PValue
#'   )
#' ) +
#'   ggplot2::geom_point(size = 2.5) +
#'   ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#'   ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
#'   scale_color_pvalue() +
#'   ggplot2::labs(
#'     title = "Standardized associations across four outcomes",
#'     subtitle = "Dashed line: p = 0.05. Dotted line: no effect.",
#'     x = "Estimate per SD", y = expression(-log[10](p))
#'   ) +
#'   ggplot2::theme_bw()
#' }
#'
#' @export
scale_color_pvalue <- function(
    palette = "inferno",
    direction = -1,
    name = "P-value",
    breaks = c(
      1,
      0.1,
      0.05,
      0.01,
      0.001,
      1e-5,
      1e-8
    ),
    labels = .format_pvalue_labels,
    limits = c(1e-8, 1),
    na.value = "grey80",
    guide = NULL,
    ...) {

  .scale_pvalue(
    aesthetics = "colour",
    palette = palette,
    direction = direction,
    name = name,
    breaks = breaks,
    labels = labels,
    limits = limits,
    na.value = na.value,
    guide = guide,
    ...
  )
}


#' Apply an evidence-aware p-value fill scale
#'
#' Maps raw p-values to either the Inferno or Viridis color palette using a
#' threshold-aware transformation. The transformation allocates additional
#' visual resolution around commonly interpreted p-value thresholds of 0.05,
#' 0.01, and 0.001. Values above 0.05 are progressively desaturated to reduce
#' their visual emphasis while retaining a continuous representation of the
#' underlying p-values. The colorbar uses the same warped coordinates, so its
#' raw p-value labels remain separated around these thresholds. Use this scale
#' when p-values are mapped to the `fill` aesthetic.
#'
#' @param palette A character value specifying the color palette. Must be
#'   `"inferno"` or `"viridis"`. Defaults to `"inferno"`.
#' @param direction A numeric value controlling the direction of the palette.
#'   Use `-1`, the default, for darker colors to indicate smaller p-values and
#'   stronger statistical evidence. Use `1` to reverse the palette.
#' @param name The legend title. Defaults to `"P-value"`.
#' @param breaks A numeric vector of raw p-values to display as legend breaks.
#'   Defaults to commonly interpreted statistical thresholds and reference
#'   values.
#' @param labels A function or character vector used to label the legend
#'   breaks. By default, the legend displays raw p-values rather than
#'   transformed values.
#' @param limits A numeric vector containing the minimum and maximum p-values
#'   represented by the scale. Values outside these limits are squished to
#'   the nearest limit. Defaults to `c(1e-8, 1)`.
#' @param na.value The fill color assigned to missing p-values. Defaults to
#'   `"grey80"`.
#' @param guide A guide function or guide name. The default, `NULL`, uses a
#'   taller colorbar so the threshold-aware breaks remain legible. Supply
#'   `"colourbar"` for ggplot2's standard-sized colorbar or another guide to
#'   override this behavior.
#' @param ... Additional arguments passed to [ggplot2::continuous_scale()].
#'
#' @return A continuous ggplot2 fill scale.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # A univariate screen: each biomarker against diagnosis, standardized so the
#' # effect sizes are comparable.
#' screen <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = "Diagnosis",
#'   predictor_vars = c(
#'     "age", "sex", "AXL", "Adiponectin", "Cortisol",
#'     "Ferritin", "Insulin", "Leptin", "tau", "p_tau"
#'   ),
#'   Standardize = TRUE
#' )
#'
#' # Bar length is the effect, fill is the evidence for it. Reading the two
#' # together is the point: a long pale bar is a large estimate nobody should
#' # rely on, and a short dark bar is a small effect that is real.
#' ggplot2::ggplot(
#'   screen$Results,
#'   ggplot2::aes(
#'     x = Estimate,
#'     y = stats::reorder(TermLabel, Estimate),
#'     fill = PValue
#'   )
#' ) +
#'   ggplot2::geom_col() +
#'   ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
#'   scale_fill_pvalue() +
#'   ggplot2::labs(
#'     title = "Association with impairment, per standard deviation",
#'     subtitle = "Dashed line: odds ratio of 1, no association",
#'     x = "Odds ratio per SD", y = NULL
#'   ) +
#'   ggplot2::theme_bw()
#' }
#'
#' @export
scale_fill_pvalue <- function(
    palette = "inferno",
    direction = -1,
    name = "P-value",
    breaks = c(
      1,
      0.1,
      0.05,
      0.01,
      0.001,
      1e-5,
      1e-8
    ),
    labels = .format_pvalue_labels,
    limits = c(1e-8, 1),
    na.value = "grey80",
    guide = NULL,
    ...) {

  .scale_pvalue(
    aesthetics = "fill",
    palette = palette,
    direction = direction,
    name = name,
    breaks = breaks,
    labels = labels,
    limits = limits,
    na.value = na.value,
    guide = guide,
    ...
  )
}
