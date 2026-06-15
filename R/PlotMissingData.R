#' Plot Missing Data
#'
#' Visualize missing data patterns with variables as rows and observations as
#' columns. Optional hover variables can be included to facilitate quality
#' control workflows when converting the plot to an interactive Plotly figure.
#'
#' @param DataFrame A data frame.
#' @param Variables Character vector of variables to visualize. If NULL,
#'   all columns except HoverVars are used.
#' @param HoverVars Optional character vector of columns to include in hover
#'   text. Useful for participant IDs, visit names, dates, sites, etc.
#' @param Relabel Logical. If TRUE, variable labels are used when available.
#' @param show_perc Logical. If TRUE, overall missingness percentages are shown
#'   in the legend.
#' @param show_perc_var Logical. If TRUE, variable-specific missingness
#'   percentages are appended to y-axis labels.
#' @param cluster Logical. If TRUE, variables are clustered by missingness
#'   pattern.
#'
#' @return A ggplot object.
#'
#' @examples
#' PlotMissingData(
#'   mtcars,
#'   HoverVars = "cyl"
#' )
#'
#' @export
PlotMissingData <- function(
    DataFrame,
    Variables = NULL,
    HoverVars = NULL,
    Relabel = TRUE,
    show_perc = TRUE,
    show_perc_var = TRUE,
    cluster = FALSE
) {

  # Validate inputs

  if (!is.data.frame(DataFrame)) {
    stop("DataFrame must be a data.frame.")
  }

  if (!is.null(HoverVars)) {

    missing_hover <- setdiff(HoverVars, names(DataFrame))

    if (length(missing_hover) > 0) {
      stop(
        "HoverVars not found in DataFrame: ",
        paste(missing_hover, collapse = ", ")
      )
    }
  }

  # Determine variables

  if (is.null(Variables)) {

    Variables <- setdiff(
      names(DataFrame),
      HoverVars
    )

  } else {

    missing_vars <- setdiff(Variables, names(DataFrame))

    if (length(missing_vars) > 0) {
      stop(
        "Variables not found in DataFrame: ",
        paste(missing_vars, collapse = ", ")
      )
    }
  }

  analysis_df <- DataFrame %>%
    dplyr::select(dplyr::all_of(Variables))

  original_vars <- names(analysis_df)

  # Apply labels

  if (Relabel) {

    analysis_df <- ReplaceMissingLabels(analysis_df)

    display_vars <- sjlabelled::get_label(analysis_df)

    names(display_vars) <- NULL

    display_vars <- make.unique(display_vars)

    names(display_vars) <- original_vars

    names(analysis_df) <- display_vars

  } else {

    display_vars <- original_vars
    names(display_vars) <- original_vars

  }

  # Missing percentages

  miss_pct <- colMeans(is.na(analysis_df)) * 100

  names(miss_pct) <- names(analysis_df)

  # Build hover text

  if (!is.null(HoverVars)) {

    hover_text <- apply(
      DataFrame[, HoverVars, drop = FALSE],
      1,
      function(x) {

        paste(
          paste(
            HoverVars,
            as.character(x),
            sep = ": "
          ),
          collapse = "<br>"
        )

      }
    )

  } else {

    hover_text <- paste0(
      "Row: ",
      seq_len(nrow(DataFrame))
    )

  }

  # Prepare data

  missing_df <- analysis_df %>%
    mutate(
      row_num = dplyr::row_number(),
      hover_text = hover_text
    ) %>%
    mutate(
      across(
        .cols = all_of(names(analysis_df)),
        .fns = is.na
      )
    )

  plot_data <- missing_df %>%
    tidyr::pivot_longer(
      cols = -c(row_num, hover_text),
      names_to = "variable",
      values_to = "is_missing"
    )

  plot_data$variable <- factor(
    plot_data$variable,
    levels = rev(names(analysis_df))
  )

  # Cluster variables

  if (cluster) {

    miss_matrix <- missing_df %>%
      select(-row_num, -hover_text) %>%
      as.matrix()

    hc <- hclust(dist(t(miss_matrix)))

    var_order <- colnames(miss_matrix)[hc$order]

    plot_data$variable <- factor(
      plot_data$variable,
      levels = rev(var_order)
    )
  }

  # Add variable percentages

  if (show_perc_var) {

    current_levels <- levels(plot_data$variable)

    pct_labels <- paste0(
      current_levels,
      " (",
      round(miss_pct[current_levels], 1),
      "%)"
    )

    levels(plot_data$variable) <- pct_labels
  }

  # Build hover field

  plot_data$hover_display <- paste0(
    plot_data$hover_text,
    "<br>Variable: ",
    as.character(plot_data$variable),
    "<br>Status: ",
    ifelse(
      plot_data$is_missing,
      "Missing",
      "Present"
    )
  )

  # Create plot

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = row_num,
      y = variable,
      fill = is_missing,
      text = hover_display
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(
      name = "",
      values = c(
        "FALSE" = "grey80",
        "TRUE" = "grey20"
      ),
      labels = if (show_perc) {

        c(
          "FALSE" = paste0(
            "Present (",
            round(
              100 - mean(plot_data$is_missing) * 100,
              1
            ),
            "%)"
          ),
          "TRUE" = paste0(
            "Missing (",
            round(
              mean(plot_data$is_missing) * 100,
              1
            ),
            "%)"
          )
        )

      } else {

        c(
          "FALSE" = "Present",
          "TRUE" = "Missing"
        )

      }
    ) +
    ggplot2::scale_x_continuous(
      name = "Observations",
      expand = c(0, 0),
      position = "top"
    ) +
    ggplot2::scale_y_discrete(
      name = "Variables",
      expand = c(0, 0)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(reverse = TRUE)
    )

  return(p)

}
