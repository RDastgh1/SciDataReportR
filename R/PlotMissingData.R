#' Plot Missing Data
#'
#' Visualize missing data patterns with variables as rows and observations as
#' columns. Optional hover variables can be included to facilitate quality
#' control workflows when converting the plot to an interactive Plotly figure.
#'
#' @param data A data frame.
#' @param variables Character vector of variables to visualize. If NULL,
#'   all columns except `HoverVars`, `x_var`, and `facet_by` are used.
#' @param HoverVars Optional character vector of columns to include in hover
#'   text. Useful for participant IDs, visit names, dates, sites, etc.
#' @param x_var Optional single column name to use for the x-axis. Numeric and
#'   date variables retain their original scale; categorical variables use a
#'   discrete axis. Missing x values are displayed as `"Missing"`.
#' @param facet_by Optional single column name used to create missingness
#'   panels. Missing facet values are displayed in a `"Missing"` panel.
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
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # The revalued data includes missingness defined in the codebook
#' vars <- c("age", "AXL", "Angiotensinogen", "BMP_6", "IL_6",
#'           "Fetuin_A", "NT_proBNP", "ENA_78")
#'
#' PlotMissingData(
#'   Labelled,
#'   variables = vars,
#'   HoverVars = "Diagnosis"
#' )
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @export
PlotMissingData <- function(data,
    variables = NULL,
    HoverVars = NULL,
    x_var = NULL,
    facet_by = NULL,
    Relabel = TRUE,
    show_perc = TRUE,
    show_perc_var = TRUE,
    cluster = FALSE,
    DataFrame = lifecycle::deprecated(),
    Variables = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMissingData(DataFrame)", "PlotMissingData(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "PlotMissingData(Variables)", "PlotMissingData(variables)")
    variables <- Variables
  }
  Variables <- variables


  # Validate inputs

  if (!is.data.frame(DataFrame)) {
    stop("DataFrame must be a data.frame.")
  }

  ValidateSingleControlVar <- function(var, arg_name) {
    if (is.null(var)) return(NULL)
    if (!is.character(var) || length(var) != 1L || is.na(var) || !nzchar(var)) {
      stop(arg_name, " must be a single, non-missing column name.")
    }
    if (!var %in% names(DataFrame)) {
      stop(arg_name, " was not found in DataFrame: ", var)
    }
    var
  }

  x_var <- ValidateSingleControlVar(x_var, "x_var")
  facet_by <- ValidateSingleControlVar(facet_by, "facet_by")

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
      unique(c(HoverVars, x_var, facet_by))
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
      hover_text = hover_text,
      .x_value = if (is.null(x_var)) row_num else DataFrame[[x_var]],
      .facet_value = if (is.null(facet_by)) "All observations" else DataFrame[[facet_by]]
    ) %>%
    mutate(
      across(
        .cols = all_of(names(analysis_df)),
        .fns = is.na
      )
    )

  plot_data <- missing_df %>%
    tidyr::pivot_longer(
      cols = -c(row_num, hover_text, .x_value, .facet_value),
      names_to = "variable",
      values_to = "is_missing"
    )

  plot_data$.facet_value <- as.character(plot_data$.facet_value)
  plot_data$.facet_value[is.na(plot_data$.facet_value) | !nzchar(plot_data$.facet_value)] <- "Missing"

  if (is.null(x_var)) {
    plot_data$.x_plot <- plot_data$row_num
    x_scale <- ggplot2::scale_x_continuous(
      name = "Observations",
      expand = c(0, 0),
      position = "top"
    )
  } else {
    x_values <- DataFrame[[x_var]]
    x_label <- sjlabelled::get_label(x_values, def.value = x_var)
    is_continuous_x <- is.numeric(x_values) || inherits(x_values, c("Date", "POSIXct", "POSIXt"))

    if (is_continuous_x && !anyNA(x_values)) {
      plot_data$.x_plot <- plot_data$.x_value
      x_scale <- if (inherits(x_values, "Date")) {
        ggplot2::scale_x_date(name = x_label, expand = c(0, 0), position = "top")
      } else if (inherits(x_values, c("POSIXct", "POSIXt"))) {
        ggplot2::scale_x_datetime(name = x_label, expand = c(0, 0), position = "top")
      } else {
        ggplot2::scale_x_continuous(name = x_label, expand = c(0, 0), position = "top")
      }
    } else {
      x_display <- as.character(plot_data$.x_value)
      x_display[is.na(x_display) | !nzchar(x_display)] <- "Missing"
      x_levels <- as.character(x_values)
      x_levels[is.na(x_levels) | !nzchar(x_levels)] <- "Missing"
      plot_data$.x_plot <- factor(x_display, levels = unique(x_levels))
      x_scale <- ggplot2::scale_x_discrete(name = x_label, expand = c(0, 0), position = "top")
    }
  }

  plot_data$variable <- factor(
    plot_data$variable,
    levels = rev(names(analysis_df))
  )

  # Cluster variables

  if (cluster) {

    miss_matrix <- missing_df %>%
      select(-row_num, -hover_text, -.x_value, -.facet_value) %>%
      as.matrix()

    hc <- hclust(dist(t(miss_matrix)))

    var_order <- colnames(miss_matrix)[hc$order]

    plot_data$variable <- factor(
      plot_data$variable,
      levels = rev(var_order)
    )
  }

  # Add variable percentages. Faceted panels receive subgroup-specific labels.
  plot_data$.variable_display <- as.character(plot_data$variable)
  if (show_perc_var) {
    if (is.null(facet_by)) {
      current_levels <- levels(plot_data$variable)
      pct_labels <- paste0(current_levels, " (", round(miss_pct[current_levels], 1), "%)")
      plot_data$.variable_display <- pct_labels[match(plot_data$.variable_display, current_levels)]
    } else {
      plot_data <- plot_data %>%
        dplyr::group_by(.facet_value, variable) %>%
        dplyr::mutate(
          .variable_display = paste0(as.character(variable), " (", round(mean(is_missing) * 100, 1), "%)")
        ) %>%
        dplyr::ungroup()
    }
  }

  variable_order <- unique(plot_data$.variable_display)
  plot_data$.variable_display <- factor(plot_data$.variable_display, levels = rev(variable_order))

  # Build hover field

  plot_data$hover_display <- paste0(
    plot_data$hover_text,
    "<br>Variable: ",
    as.character(plot_data$.variable_display),
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
      x = .x_plot,
      y = .variable_display,
      fill = is_missing,
      text = hover_display
    )
  ) +
    ggplot2::geom_tile() +
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
    x_scale +
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

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.facet_value),
      scales = "free_y"
    )
  }

  return(p)

}
