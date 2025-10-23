#' Create Kynurenine-Tryptophan Pathway Plot
#'
#' Creates a pathway diagram for the kynurenine-tryptophan metabolic pathway with
#' color-coded fold changes or correlations and significance indicators.
#'
#' @param results_table Data frame with columns: Metabolite, p_value, p_adj, and either "% Change" or "correlation"
#' @param title Character string for plot title
#' @param value_type Character string: "auto", "fold_change", or "correlation"
#' @param metabolite_mapping Named character vector mapping results table names to standard names.
#'   For example: c("N'-Formylkynurenine" = "N-Formylkynurenine", "Quinolinic Acid(log10)" = "Quinolinic Acid")
#' @param use_fdr Logical: if TRUE uses FDR-adjusted p-values (p_adj) for significance, if FALSE uses raw p-values. Default is FALSE.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with raw p-values
#' plot <- CreatePathwayPlot_KT(results, "My Pathway")
#'
#' # Use FDR-adjusted p-values for significance
#' plot <- CreatePathwayPlot_KT(results, "My Pathway", use_fdr = TRUE)
#'
#' # With custom metabolite name mapping
#' name_map <- c(
#'   "N'-Formylkynurenine" = "N-Formylkynurenine",
#'   "Quinolinic Acid(log10)" = "Quinolinic Acid",
#'   "3-OH-kynurenine" = "3-Hydroxykynurenine"
#' )
#' plot <- CreatePathwayPlot_KT(results, "My Pathway", metabolite_mapping = name_map)
#' }
CreatePathwayPlot_KT <- function(results_table,
                                 title = "",
                                 value_type = "auto",
                                 metabolite_mapping = NULL,
                                 use_fdr = FALSE) {

  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)

  # Auto-detect value type
  if (value_type == "auto") {
    if ("% Change" %in% names(results_table)) {
      value_type <- "fold_change"
      value_col <- "% Change"
      scale_limits <- c(-10, 10)
      scale_name <- "% Change"
    } else if ("correlation" %in% names(results_table)) {
      value_type <- "correlation"
      value_col <- "correlation"
      scale_limits <- c(-1, 1)
      scale_name <- "Correlation (r)"
    } else {
      stop("Could not auto-detect value type. Results table must contain either '% Change' or 'correlation' column.")
    }
  } else if (value_type == "fold_change") {
    value_col <- "% Change"
    scale_limits <- c(-10, 10)
    scale_name <- "% Change"
  } else if (value_type == "correlation") {
    value_col <- "correlation"
    scale_limits <- c(-1, 1)
    scale_name <- "Correlation (r)"
  }

  # Define metabolite positions
  pathway_positions <- data.frame(
    Metabolite = c("Tryptophan",
                   "Serotonin",
                   "N-Formylkynurenine",
                   "Kynurenine",
                   "Kynurenic Acid",
                   "3-Hydroxykynurenine",
                   "Anthranilic Acid",
                   "Xanthurenic Acid",
                   "3-Hydroxyanthranilic acid",
                   "Quinolinic Acid"),
    x = c(6, 3, 9, 9, 5, 9, 13, 9, 13, 13),
    y = c(10, 8, 8, 6, 4, 4, 4, 2, 2, 0),
    width = c(2.8, 2.8, 3.4, 2.8, 2.8, 3.4, 2.8, 2.8, 3.8, 2.8),
    height = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7),
    stringsAsFactors = FALSE
  )

  # Handle different metabolite names
  results_table_adj <- results_table

  # Default mapping if none provided
  if (is.null(metabolite_mapping)) {
    metabolite_mapping <- c(
      "N'-Formylkynurenine" = "N-Formylkynurenine"
    )
  }

  # Apply metabolite name mapping
  for (old_name in names(metabolite_mapping)) {
    if (old_name %in% results_table_adj$Metabolite) {
      results_table_adj$Metabolite[results_table_adj$Metabolite == old_name] <- metabolite_mapping[old_name]
    }
  }

  # Merge with results
  plot_data <- pathway_positions %>%
    dplyr::left_join(results_table_adj, by = "Metabolite")

  # Replace NA values with 0 for color mapping and ensure no missing data
  if (value_col %in% names(plot_data)) {
    plot_data[[value_col]][is.na(plot_data[[value_col]])] <- 0
  }
  plot_data$p_value[is.na(plot_data$p_value)] <- 1  # Set NA p-values to 1
  plot_data$p_adj[is.na(plot_data$p_adj)] <- 1  # Set NA p-values to 1

  # Choose which p-value to use based on use_fdr parameter
  if (use_fdr) {
    # Check if p_adj column exists
    if (!"p_adj" %in% names(plot_data)) {
      warning("p_adj column not found in results table. Using p_value instead.")
      p_col <- "p_value"
    } else {
      p_col <- "p_adj"
    }
  } else {
    # Check if p_value column exists
    if (!"p_value" %in% names(plot_data)) {
      stop("p_value column not found in results table.")
    }
    p_col <- "p_value"
  }

  # Add plot aesthetics
  plot_data <- plot_data %>%
    dplyr::mutate(
      color_value = if (value_col %in% names(.)) .data[[value_col]] else 0,
      color_value = pmax(pmin(color_value, scale_limits[2]), scale_limits[1]),
      border_width = ifelse(!is.na(.data[[p_col]]) & .data[[p_col]] < 0.05, 2, 0.8),
      border_color = "black"
    )

  # Remove any rows with NA positions (shouldn't happen but just in case)
  plot_data <- plot_data[complete.cases(plot_data[c("x", "y")]), ]

  # Create the plot
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    # Add arrows manually with proper positioning
    # Tryptophan to Serotonin
    annotate("segment", x = 6, y = 9.65, xend = 3, yend = 8.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Tryptophan to N-Formylkynurenine
    annotate("segment", x = 6, y = 9.65, xend = 9, yend = 8.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # N-Formylkynurenine to Kynurenine - STRAIGHT DOWN
    annotate("segment", x = 9, y = 7.65, xend = 9, yend = 6.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Kynurenine to Kynurenic Acid
    annotate("segment", x = 9, y = 5.65, xend = 5, yend = 4.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Kynurenine to 3-Hydroxykynurenine - STRAIGHT DOWN
    annotate("segment", x = 9, y = 5.65, xend = 9, yend = 4.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Kynurenine to Anthranilic Acid
    annotate("segment", x = 9, y = 5.65, xend = 13, yend = 4.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # 3-Hydroxykynurenine to Xanthurenic Acid - STRAIGHT DOWN
    annotate("segment", x = 9, y = 3.65, xend = 9, yend = 2.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # 3-Hydroxykynurenine to 3-Hydroxyanthranilic acid
    annotate("segment", x = 9, y = 3.65, xend = 13, yend = 2.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Anthranilic Acid to 3-Hydroxyanthranilic acid - STRAIGHT DOWN
    annotate("segment", x = 13, y = 3.65, xend = 13, yend = 2.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # 3-Hydroxyanthranilic acid to Quinolinic Acid - STRAIGHT DOWN
    annotate("segment", x = 13, y = 1.65, xend = 13, yend = 0.35,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black", linewidth = 0.8) +
    # Add metabolite boxes
    geom_tile(aes(fill = color_value, width = width, height = height),
              color = plot_data$border_color,
              linewidth = plot_data$border_width) +
    # Add metabolite labels
    geom_text(aes(label = Metabolite),
              size = 3, fontface = "plain") +
    # Color scale
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = scale_limits,
                         name = scale_name,
                         na.value = "grey90") +
    # Theme
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 1)
    ) +
    # Set axis limits
    xlim(1, 15) +
    ylim(-1, 11) +
    # Add title
    ggtitle(title) +
    # Ensure aspect ratio
    coord_fixed(ratio = 1.2)

  return(p)
}
