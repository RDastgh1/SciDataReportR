#' Plot Missing Data
#'
#' A custom function to plot missing data with variables as rows and observations as columns.
#'
#' @param DataFrame The dataframe containing missing data.
#' @param Variables A character vector specifying the variables to be included in the plot. If NULL (default), all variables in the dataframe will be used.
#' @param Relabel Logical indicating whether to relabel the columns based on their labels. Default is TRUE.
#' @param show_perc Logical indicating whether to show percentage labels in legend. Default is TRUE.
#' @param show_perc_var Logical indicating whether to show percentage of missing for each variable. Default is TRUE.
#' @param cluster Logical indicating whether to cluster rows by missingness pattern. Default is FALSE.
#' @return A ggplot object visualizing missing data with variables as rows.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
#'
PlotMissingData <- function(DataFrame,
                            Variables = NULL,
                            Relabel = TRUE,
                            show_perc = TRUE,
                            show_perc_var = TRUE,
                            cluster = FALSE){

  # If Variables argument is not provided, use all columns of the DataFrame
  if(is.null(Variables)){
    Variables = colnames(DataFrame)
  }

  # Subset DataFrame to include only specified Variables
  DataFrame <- DataFrame[Variables]

  # Store original variable names
  original_vars <- colnames(DataFrame)

  # Relabel columns based on their labels if Relabel is TRUE
  if(Relabel){
    DataFrame <- ReplaceMissingLabels(DataFrame)
    labels <- get_label(DataFrame)
    names(labels) <- NULL
    labels <- make.unique(labels)
    colnames(DataFrame) <- labels
    display_vars <- labels
  } else {
    display_vars <- original_vars
  }

  # Calculate missing percentages for each variable
  miss_pct <- colMeans(is.na(DataFrame)) * 100

  # Create a binary missing data frame first (before pivoting)
  # This avoids type conflicts
  missing_df <- DataFrame %>%
    mutate(across(everything(), is.na)) %>%
    mutate(row_num = row_number())

  # Now pivot the binary data
  plot_data <- missing_df %>%
    pivot_longer(cols = -row_num,
                 names_to = "variable",
                 values_to = "is_missing") %>%
    mutate(variable = factor(variable, levels = rev(display_vars))) # Reverse for correct orientation

  # If clustering is requested, reorder variables by missingness pattern
  if(cluster){
    # Create a binary matrix of missingness
    miss_matrix <- as.matrix(missing_df %>% select(-row_num))
    # Perform hierarchical clustering on variables (transposed)
    hc <- hclust(dist(t(miss_matrix)))
    # Reorder variables based on clustering
    var_order <- display_vars[hc$order]
    plot_data$variable <- factor(plot_data$variable, levels = rev(var_order))
  }

  # Create variable labels with percentages if requested
  if(show_perc_var){
    var_labels <- paste0(display_vars, " (", round(miss_pct, 1), "%)")
    names(var_labels) <- display_vars

    # Update factor levels with percentage labels
    current_levels <- levels(plot_data$variable)
    new_levels <- var_labels[current_levels]
    levels(plot_data$variable) <- new_levels
  }

  # Create the plot
  p <- ggplot(plot_data, aes(x = row_num, y = variable, fill = is_missing)) +
    geom_raster() +
    scale_fill_manual(
      name = "",
      values = c("FALSE" = "grey80", "TRUE" = "grey20"),
      labels = if(show_perc) {
        c("FALSE" = paste0("Present (", round(100 - mean(plot_data$is_missing) * 100, 1), "%)"),
          "TRUE" = paste0("Missing (", round(mean(plot_data$is_missing) * 100, 1), "%)"))
      } else {
        c("FALSE" = "Present", "TRUE" = "Missing")
      }
    ) +
    scale_x_continuous(
      name = "Observations",
      expand = c(0, 0),
      position = "top"
    ) +
    scale_y_discrete(
      name = "Variables",
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(fill = guide_legend(reverse = TRUE))

  return(p)
}
