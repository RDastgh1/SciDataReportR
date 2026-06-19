#' Plot cluster boxplots by variable
#'
#' Create grouped cluster boxplots where clusters are shown on the x-axis and
#' selected variables are shown as filled boxplots within each cluster. This is
#' useful for visualizing cognitive, clinical, biomarker, or domain profiles
#' across discovered clusters or subgroup solutions.
#'
#' Variables can optionally be z-scored across all participants before plotting.
#' Reference lines can be automatically added for z-score or T-score style
#' interpretation, or omitted entirely for custom overlays. By default,
#' reference lines are only added automatically when `Scale = TRUE`.
#'
#' Variable labels are used by default when available from a codebook or from
#' variable label attributes.
#'
#' @param Data A data frame.
#' @param ClusterVar Character string naming the cluster/grouping variable.
#' @param Variables Character vector of variable names to plot.
#' @param Codebook Optional codebook data frame with columns `Variable` and
#'   `Label`.
#' @param Scale Logical. If `TRUE`, variables are z-scored across all
#'   participants before plotting. Default is `FALSE`.
#' @param ScoreType Character. One of `"auto"`, `"z"`, `"t"`, or `"raw"`.
#'   Used for axis labeling and reference-line defaults.
#' @param ReferenceLines Character. One of `"auto"`, `"z"`, `"t"`, or
#'   `"none"`.
#'
#'   If `"auto"`, z-score reference lines are added only when `Scale = TRUE`.
#'
#'   If `"z"`, lines are added at `-1`, `-0.5`, `0`, `0.5`, and `1`.
#'
#'   If `"t"`, lines are added at `40`, `45`, `50`, `55`, and `60`.
#'
#'   If `"none"`, no reference lines are added so users can customize
#'   overlays manually using additional ggplot layers.
#' @param ClusterLabel Character. One of `"n_percent"`, `"n"`, or `"none"`.
#'
#'   `"n_percent"` adds cluster sample size and percent to the x-axis label.
#'
#'   `"n"` adds only sample size.
#'
#'   `"none"` uses only the cluster name/value.
#' @param Relabel Logical. If `TRUE`, variable labels are used when available.
#'   Labels are pulled first from `Codebook`, then from variable label
#'   attributes. Default is `TRUE`.
#' @param FillTitle Character string used as the fill legend title. Default is
#'   `"Test"`.
#' @param Palette Optional character vector of colors used for the fill scale.
#'   If `NULL`, a stable 20-color default palette is used. If more than 20
#'   variables are plotted, colors are interpolated from the default palette so
#'   the function does not fail.
#' @param YLabel Optional y-axis label. If `NULL`, an appropriate label is
#'   chosen automatically from `Scale` and `ScoreType`.
#' @param BoxplotWidth Numeric width passed to
#'   `ggplot2::geom_boxplot()`. Default is `0.75`.
#' @param OutlierSize Numeric size of outlier points. Default is `0.8`.
#' @param BaseSize Base font size for the plot theme. Default is `14`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' PlotClusterBoxplot(
#'   Data = mtcars,
#'   ClusterVar = "cyl",
#'   Variables = c("mpg", "disp", "hp"),
#'   Scale = TRUE,
#'   ReferenceLines = "z"
#' )
#'
#' PlotClusterBoxplot(
#'   Data = mtcars,
#'   ClusterVar = "cyl",
#'   Variables = c("mpg", "disp", "hp"),
#'   ClusterLabel = "n"
#' )
#'
#' @export
PlotClusterBoxplot <- function(
    Data,
    ClusterVar,
    Variables,
    Codebook = NULL,
    Scale = FALSE,
    ScoreType = c("auto", "z", "t", "raw"),
    ReferenceLines = c("auto", "z", "t", "none"),
    ClusterLabel = c("n_percent", "n", "none"),
    Relabel = TRUE,
    FillTitle = "Test",
    Palette = NULL,
    YLabel = NULL,
    BoxplotWidth = 0.75,
    OutlierSize = 0.8,
    BaseSize = 14
) {

  ScoreType <- match.arg(ScoreType)
  ReferenceLines <- match.arg(ReferenceLines)
  ClusterLabel <- match.arg(ClusterLabel)

  # Validate inputs

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }

  if (!is.character(ClusterVar) || length(ClusterVar) != 1) {
    stop("ClusterVar must be a single character string.")
  }

  if (!ClusterVar %in% names(Data)) {
    stop(
      "ClusterVar was not found in Data. Expected column: ",
      ClusterVar
    )
  }

  if (!is.character(Variables) || length(Variables) == 0) {
    stop("Variables must be a non-empty character vector.")
  }

  missing_vars <- setdiff(Variables, names(Data))

  if (length(missing_vars) > 0) {
    warning(
      "The following Variables were not found in Data and were removed: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  Variables <- intersect(Variables, names(Data))

  if (length(Variables) == 0) {
    stop("No valid Variables remain after checking against Data.")
  }

  numeric_vars <- Variables[vapply(Data[Variables], is.numeric, logical(1))]
  invalid_vars <- setdiff(Variables, numeric_vars)

  if (length(invalid_vars) > 0) {
    warning(
      "The following Variables were removed because they are not numeric: ",
      paste(invalid_vars, collapse = ", "),
      call. = FALSE
    )
  }

  Variables <- numeric_vars

  if (length(Variables) == 0) {
    stop("No numeric Variables remain after removing invalid variables.")
  }

  if (!is.null(Codebook)) {

    if (!is.data.frame(Codebook)) {
      stop("Codebook must be a data frame when supplied.")
    }

    required_codebook_cols <- c("Variable", "Label")
    missing_codebook_cols <- setdiff(required_codebook_cols, names(Codebook))

    if (length(missing_codebook_cols) > 0) {
      stop(
        "Codebook must contain columns Variable and Label. Missing: ",
        paste(missing_codebook_cols, collapse = ", ")
      )
    }
  }

  if (!is.logical(Scale) || length(Scale) != 1 || is.na(Scale)) {
    stop("Scale must be TRUE or FALSE.")
  }

  if (!is.logical(Relabel) || length(Relabel) != 1 || is.na(Relabel)) {
    stop("Relabel must be TRUE or FALSE.")
  }

  if (!is.character(FillTitle) || length(FillTitle) != 1) {
    stop("FillTitle must be a single character string.")
  }

  if (!is.null(Palette) && !is.character(Palette)) {
    stop("Palette must be NULL or a character vector of color values.")
  }

  if (!is.null(YLabel) && (!is.character(YLabel) || length(YLabel) != 1)) {
    stop("YLabel must be NULL or a single character string.")
  }

  if (!is.numeric(BoxplotWidth) ||
      length(BoxplotWidth) != 1 ||
      is.na(BoxplotWidth) ||
      BoxplotWidth <= 0) {
    stop("BoxplotWidth must be a single positive number.")
  }

  if (!is.numeric(OutlierSize) ||
      length(OutlierSize) != 1 ||
      is.na(OutlierSize) ||
      OutlierSize < 0) {
    stop("OutlierSize must be a single non-negative number.")
  }

  if (!is.numeric(BaseSize) ||
      length(BaseSize) != 1 ||
      is.na(BaseSize) ||
      BaseSize <= 0) {
    stop("BaseSize must be a single positive number.")
  }

  # Prepare data

  plot_data <- Data %>%
    dplyr::select(dplyr::all_of(c(ClusterVar, Variables))) %>%
    dplyr::mutate(
      .ClusterPlot = .data[[ClusterVar]]
    )

  if (Scale) {
    plot_data <- plot_data %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(Variables),
          ~ as.numeric(scale(.x))
        )
      )
  }

  if (!is.factor(plot_data$.ClusterPlot)) {

    cluster_levels <- plot_data$.ClusterPlot %>%
      unique() %>%
      stats::na.omit() %>%
      sort()

    plot_data <- plot_data %>%
      dplyr::mutate(
        .ClusterPlot = factor(
          .data$.ClusterPlot,
          levels = cluster_levels
        )
      )
  }

  cluster_counts <- plot_data %>%
    dplyr::filter(!is.na(.data$.ClusterPlot)) %>%
    dplyr::count(.ClusterPlot, name = "n") %>%
    dplyr::mutate(
      Percent = 100 * .data$n / sum(.data$n),
      ClusterBase = as.character(.data$.ClusterPlot)
    )

  if (ClusterLabel == "n_percent") {

    cluster_counts <- cluster_counts %>%
      dplyr::mutate(
        ClusterDisplay = paste0(
          .data$ClusterBase,
          " (n = ",
          .data$n,
          ", ",
          sprintf("%.0f", .data$Percent),
          "%)"
        )
      )

  } else if (ClusterLabel == "n") {

    cluster_counts <- cluster_counts %>%
      dplyr::mutate(
        ClusterDisplay = paste0(
          .data$ClusterBase,
          " (n = ",
          .data$n,
          ")"
        )
      )

  } else {

    cluster_counts <- cluster_counts %>%
      dplyr::mutate(
        ClusterDisplay = .data$ClusterBase
      )
  }

  cluster_label_lookup <- stats::setNames(
    cluster_counts$ClusterDisplay,
    cluster_counts$ClusterBase
  )

  long_data <- plot_data %>%
    dplyr::select(dplyr::all_of(c(".ClusterPlot", Variables))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(Variables),
      names_to = "Variable",
      values_to = "Score"
    )

  # Apply labels

  label_lookup <- stats::setNames(Variables, Variables)

  if (Relabel) {

    if (!is.null(Codebook)) {

      codebook_labels <- Codebook %>%
        dplyr::filter(.data$Variable %in% Variables) %>%
        dplyr::mutate(
          Label = dplyr::if_else(
            is.na(.data$Label) | .data$Label == "",
            .data$Variable,
            as.character(.data$Label)
          )
        ) %>%
        dplyr::select(dplyr::all_of(c("Variable", "Label")))

      label_lookup[codebook_labels$Variable] <- codebook_labels$Label
    }

    for (this_var in Variables) {

      if (label_lookup[[this_var]] == this_var) {

        attr_label <- attr(Data[[this_var]], "label", exact = TRUE)

        if (!is.null(attr_label) &&
            length(attr_label) == 1 &&
            !is.na(attr_label) &&
            attr_label != "") {

          label_lookup[[this_var]] <- as.character(attr_label)
        }
      }
    }
  }

  long_data <- long_data %>%
    dplyr::mutate(
      Label = unname(label_lookup[.data$Variable]),
      Label = factor(
        .data$Label,
        levels = unname(label_lookup[Variables])
      )
    )

  # Resolve score type

  if (ScoreType == "auto") {

    if (Scale) {

      ScoreType <- "z"

    } else {

      all_values <- long_data$Score

      value_mean <- mean(all_values, na.rm = TRUE)
      value_sd <- stats::sd(all_values, na.rm = TRUE)

      ScoreType <- dplyr::case_when(
        is.finite(value_mean) &&
          is.finite(value_sd) &&
          abs(value_mean) < 0.5 &&
          value_sd > 0.5 &&
          value_sd < 2 ~ "z",

        is.finite(value_mean) &&
          is.finite(value_sd) &&
          value_mean > 35 &&
          value_mean < 65 &&
          value_sd > 5 &&
          value_sd < 15 ~ "t",

        TRUE ~ "raw"
      )
    }
  }

  # Resolve reference lines

  if (ReferenceLines == "auto") {

    if (Scale) {
      ReferenceLines <- "z"
    } else {
      ReferenceLines <- "none"
    }
  }

  # Resolve y-axis label

  if (is.null(YLabel)) {

    YLabel <- dplyr::case_when(
      Scale ~ "Scaled score",
      ScoreType == "z" ~ "Z-score",
      ScoreType == "t" ~ "T-score",
      TRUE ~ "Score"
    )
  }

  # Resolve fill palette

  fill_levels <- levels(long_data$Label)
  n_fill_levels <- length(fill_levels)

  default_palette <- c(
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF"
  )

  if (is.null(Palette)) {

    if (n_fill_levels <= length(default_palette)) {

      Palette <- default_palette[seq_len(n_fill_levels)]

    } else {

      Palette <- grDevices::colorRampPalette(default_palette)(
        n_fill_levels
      )
    }

  } else {

    Palette <- rep(Palette, length.out = n_fill_levels)
  }

  names(Palette) <- fill_levels

  # Build plot

  p <- ggplot2::ggplot(
    long_data,
    ggplot2::aes(
      x = .data$.ClusterPlot,
      y = .data$Score,
      fill = .data$Label
    )
  ) +
    ggplot2::geom_boxplot(
      width = BoxplotWidth,
      outlier.size = OutlierSize,
      outlier.alpha = 0.8
    ) +
    ggplot2::scale_x_discrete(
      labels = cluster_label_lookup,
      drop = FALSE,
      na.translate = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = Palette,
      name = FillTitle
    ) +
    ggplot2::labs(
      x = "Cluster",
      y = YLabel
    ) +
    ggplot2::theme_bw(base_size = BaseSize) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_line(
        linewidth = 0.2,
        color = "grey92"
      ),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.3,
        color = "grey88"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(0.5, "cm")
    )

  # Add reference lines

  if (ReferenceLines == "z") {

    p <- p +
      ggplot2::geom_hline(
        yintercept = c(-0.5, 0.5),
        linetype = "dashed",
        linewidth = 0.45,
        color = "grey35"
      ) +
      ggplot2::geom_hline(
        yintercept = c(-1, 1),
        linetype = "dashed",
        linewidth = 0.8,
        color = "grey20"
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        linewidth = 0.45,
        color = "grey50"
      )
  }

  if (ReferenceLines == "t") {

    p <- p +
      ggplot2::geom_hline(
        yintercept = c(45, 55),
        linetype = "dashed",
        linewidth = 0.45,
        color = "grey35"
      ) +
      ggplot2::geom_hline(
        yintercept = c(40, 60),
        linetype = "dashed",
        linewidth = 0.8,
        color = "grey20"
      ) +
      ggplot2::geom_hline(
        yintercept = 50,
        linetype = "dotted",
        linewidth = 0.45,
        color = "grey50"
      )
  }

  return(p)
}
