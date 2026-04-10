#' Plot Continuous Distributions
#'
#' Creates rain-cloud plots (half-violin + box/median + scatter) for one or
#' more continuous variables, with optional group-wise colouring.
#'
#' @param DataFrame A data frame containing the variables to be plotted.
#' @param Variables Character vector of column names to plot.
#' @param Fill Optional column name for grouping.
#' @param Relabel Logical; use variable labels when available.
#' @param FacetLabelStyle One of "both", "label_only", "variable_only", "auto".
#' @param ncol Number of columns in the facet grid.
#' @param Ordinal Logical; include labelled-ordinal variables as numeric.
#'
#' @return A ggplot object.
#' @export
PlotContinuousDistributions <- function(
    DataFrame,
    Variables = NULL,
    Fill      = NULL,
    Relabel   = TRUE,
    FacetLabelStyle = c("both", "label_only", "variable_only", "auto"),
    ncol      = 3,
    Ordinal   = TRUE
) {

  FacetLabelStyle <- match.arg(FacetLabelStyle)

  # Select Variables
  if (is.null(Variables)) {
    Variables <- getNumVars(DataFrame, Ordinal = FALSE)
    if (Ordinal) Variables <- getNumVars(DataFrame, Ordinal = TRUE)
  }

  # Convert and Relabel
  if (Ordinal) {
    OriginalLabels <- sjlabelled::get_label(DataFrame, def.value = colnames(DataFrame))
    DataFrame      <- ConvertOrdinalToNumeric(DataFrame, Variables)
    DataFrame[Variables] <- lapply(DataFrame[Variables], as.numeric)

    for (v in Variables) {
      DataFrame[[v]] <- sjlabelled::set_label(DataFrame[[v]], OriginalLabels[[v]])
    }
  }

  # Build facet labels
  var_labels <- sjlabelled::get_label(DataFrame[Variables], def.value = Variables)

  build_label <- function(var, label) {

    if (!Relabel) return(var)

    if (is.null(label) || is.na(label) || label == "") {
      return(var)
    }

    switch(
      FacetLabelStyle,
      both = paste0(var, "\n", label),
      label_only = label,
      variable_only = var,
      auto = if (identical(var, label)) var else label
    )
  }

  facetlabels <- mapply(build_label, Variables, var_labels, USE.NAMES = FALSE)

  # Pivot Longer
  long_vars <- c(Variables, Fill)

  ContData  <- DataFrame %>%
    dplyr::select(dplyr::all_of(long_vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(Variables)) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(Mean = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup()

  ContData$name <- factor(
    ContData$name,
    levels = Variables,
    labels = facetlabels
  )

  # Fill handling
  if (is.null(Fill)) {
    ContData[[".fill_const"]] <- "1"
    fill_var <- ".fill_const"
    legend   <- FALSE
  } else {
    fill_var <- Fill
    legend   <- TRUE
  }

  # Plot
  p <- ggplot2::ggplot(
    ContData,
    ggplot2::aes(
      y = value,
      x = 1,
      fill   = .data[[fill_var]],
      colour = .data[[fill_var]]
    )
  ) +
    ggrain::geom_rain(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ name, scales = "free", ncol = ncol) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  # Legend handling
  if (!legend) {
    p <- p +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::scale_fill_manual(values = "#6EC259") +
      ggplot2::scale_colour_manual(values = "#6EC259")
  }

  return(p)
}
