#' Plot Continuous Distributions
#'
#' Creates rain-cloud plots (half-violin + box/median + scatter) for one or
#' more continuous variables, with optional group-wise colouring.
#'
#' @param DataFrame A data frame containing the variables to be plotted.
#' @param Variables Character vector of column names to plot.  If NULL, all
#'        numeric variables (or numeric + ordinal if Ordinal = TRUE) are used.
#' @param Fill      *Optional.*  Single column name to map to the `fill`/`colour`
#'        aesthetic (e.g., a factor such as sex, treatment group, etc.).
#' @param Relabel   Logical; if TRUE, facet labels use variable labels when
#'        available (sjlabelled).
#' @param ncol      Number of columns in the facet grid.
#' @param Ordinal   Logical; include labelled-ordinal variables as numeric?
#' @return          A `ggplot` object.
#'
#' @importFrom dplyr mutate pivot_longer select all_of group_by ungroup
#' @importFrom ggplot2 ggplot aes theme_bw guides coord_flip facet_wrap
#'         scale_fill_manual scale_colour_manual element_blank
#' @suggests ggrain, sjlabelled, rlang
#' @export
PlotContinuousDistributions <- function(
    DataFrame,
    Variables = NULL,
    Fill      = NULL,
    Relabel   = TRUE,
    ncol      = 3,
    Ordinal   = TRUE
) {

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

  facetlabels <- if (isTRUE(Relabel)) {
    createFacetLabels(dplyr::select(DataFrame, dplyr::all_of(Variables)))
  } else {
    Variables
  }

#Pivot Longer

  long_vars <- c(Variables, Fill)
  ContData  <- DataFrame %>%
    dplyr::select(dplyr::all_of(long_vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(Variables)) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(Mean = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup()

  ContData$name <- factor(ContData$name,
                          levels = Variables,
                          labels = facetlabels)


# Set up Fill/Color aesthetics

  if (is.null(Fill)) {
    # Single colour fallback
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
    ggplot2::aes(y = value, x = 1,
                 fill   = .data[[fill_var]],
                 colour = .data[[fill_var]])
  ) +
    ggrain::geom_rain(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ name, scales = "free", ncol = ncol) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  #  Legend & colour handling
  if (!legend) {
    p <- p +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::scale_fill_manual(values = "#6EC259") +
      ggplot2::scale_colour_manual(values = "#6EC259")
  }

  return(p)
}
