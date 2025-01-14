#' Create a Forest Plot from Univariate Regression Tables
#'
#' This function generates a forest plot from a list of formatted univariate regression tables.
#'
#' @param UnivariateRegressionTables A list containing formatted regression tables and styling information. Expected structure:
#'   - `FormattedTable$tbls`: A list of tables, each containing a `table_body` dataframe.
#'   - `LargeTable$table_styling$header`: A dataframe with a `label` and `spanning_header` column for headers.
#' @param pSize Numeric. Size of the points in the plot. Default is 2.
#' @return A ggplot object representing the forest plot.
#' @examples
#' # Example usage:
#' # plotForestFromTable(UnivariateRegressionTables)
#' @export
plotForestFromTable <- function(UnivariateRegressionTables, pSize = 2) {

  # Extract tables and headers
  list_Tables <- UnivariateRegressionTables$FormattedTable$tbls
  title_Tables <- UnivariateRegressionTables$LargeTable$table_styling$header %>%
    filter(label == "tbl_id1") %>%
    pull(spanning_header)

  # Combine all tables into a single dataframe
  for (i in seq_along(list_Tables)) {
    d <- list_Tables[[i]]$table_body
    d$Header <- title_Tables[i]
    if (i == 1) {
      df_Combined <- d
    } else {
      df_Combined <- rbind(df_Combined, d)
    }
  }

  # Mark significant results
  df_Combined$sig <- df_Combined$p.value < 0.05

  # Update labels for categorical variables
  df_Combined <- df_Combined %>%
    filter(!is.na(estimate)) %>%
    mutate(
      var_label = if_else(
        var_type == "categorical",
        paste0(var_label, " : ", label), # Combine var_label and label
        var_label # Keep var_label unchanged for other cases
      )
    )

  # Reverse the order of var_label for plotting
  df_Combined$var_label <- factor(df_Combined$var_label, levels = rev(unique(df_Combined$var_label)))

  # Create the forest plot
  p <- df_Combined %>%
    ggplot(aes(x = estimate, y = var_label, color = sig)) +
    geom_point(size = pSize) +  # Plot the point for the estimate
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +  # Add the error bars
    facet_wrap(~Header, nrow = 1) +  # Facet by Header
    geom_vline(xintercept = 0, linetype = "dashed") +  # Add a dashed line at 0
    labs(
      x = "Estimate",
      y = "Variable"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank()
    ) +
    scale_color_manual(values = c("darkgrey", "black"))+ theme(legend.position="none")

  return(p)
}
