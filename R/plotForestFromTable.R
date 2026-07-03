#' Create a Forest Plot from Univariate Regression Tables
#'
#' This function generates a forest plot from a list of formatted univariate regression tables.
#'
#' @param UnivariateRegressionTables A list containing formatted regression tables and styling information. Expected structure:
#'   - `FormattedTable$tbls`: A list of tables, each containing a `table_body` dataframe.
#'   - `LargeTable$table_styling$header`: A dataframe with a `label` and `spanning_header` column for headers.
#' @param pSize Numeric. Size of the points in the plot. Default is 2.
#' @param Flip Logical. If `FALSE`, outcomes are facets and predictors/terms
#'   are rows. If `TRUE`, predictors/terms are facets and outcomes are rows.
#' @return A ggplot object representing the forest plot.
#' @examples
#' # Example usage:
#' # plotForestFromTable(UnivariateRegressionTables)
#' @export
plotForestFromTable <- function(UnivariateRegressionTables, pSize = 2, Flip = FALSE) {

  if (!is.logical(Flip) || length(Flip) != 1 || is.na(Flip)) {
    stop("Flip must be TRUE or FALSE.")
  }

  # Extract tables and headers
  list_Tables <- UnivariateRegressionTables$FormattedTable$tbls
  title_Tables <- names(list_Tables)
  if (is.null(title_Tables) || any(is.na(title_Tables)) || any(title_Tables == "")) {
    title_Tables <- paste0("Outcome ", seq_along(list_Tables))
  }
  label_Tables <- title_Tables
  spanner_labels <- ScidrForestPlotSpannerLabels(
    UnivariateRegressionTables$FormattedTable,
    length(list_Tables)
  )
  if (all(!is.na(spanner_labels)) && all(spanner_labels != "")) {
    label_Tables <- spanner_labels
  }
  outcome_metadata <- UnivariateRegressionTables$Metadata$Outcomes
  if (!is.null(outcome_metadata) &&
      all(c("Outcome", "OutcomeLabel") %in% names(outcome_metadata))) {
    outcome_label_lookup <- stats::setNames(
      as.character(outcome_metadata$OutcomeLabel),
      as.character(outcome_metadata$Outcome)
    )
    matched_labels <- outcome_label_lookup[title_Tables]
    use_metadata_label <- !is.na(matched_labels) &
      matched_labels != "" &
      (matched_labels != title_Tables | label_Tables == title_Tables)
    label_Tables <- ifelse(
      use_metadata_label,
      matched_labels,
      label_Tables
    )
  }

  # Combine all tables into a single dataframe
  for (i in seq_along(list_Tables)) {
    d <- list_Tables[[i]]$table_body
    d$Header <- title_Tables[i]
    d$HeaderLabel <- label_Tables[i]
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
      TermLabel = if_else(
        !is.na(label) & label != var_label,
        paste0(var_label, " : ", label), # Combine var_label and label
        var_label # Keep var_label unchanged for other cases
      ),
      OutcomeLabel = HeaderLabel
    )

  # Reverse the row label order for plotting
  if (Flip) {
    df_Combined$PlotRow <- factor(
      df_Combined$OutcomeLabel,
      levels = rev(label_Tables)
    )
    df_Combined$PlotFacet <- factor(
      df_Combined$TermLabel,
      levels = unique(df_Combined$TermLabel)
    )
    y_label <- "Outcome"
  } else {
    df_Combined$PlotRow <- factor(
      df_Combined$TermLabel,
      levels = rev(unique(df_Combined$TermLabel))
    )
    df_Combined$PlotFacet <- factor(
      df_Combined$OutcomeLabel,
      levels = label_Tables
    )
    y_label <- "Variable"
  }

  # Create the forest plot
  p <- df_Combined %>%
    ggplot(aes(x = estimate, y = PlotRow, color = sig)) +
    geom_point(size = pSize) +  # Plot the point for the estimate
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.2) +  # Add the error bars
    facet_wrap(~PlotFacet, nrow = 1) +  # Facet by outcome or term
    geom_vline(xintercept = 0, linetype = "dashed") +  # Add a dashed line at 0
    labs(
      x = "Estimate",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank()
    ) +
    scale_color_manual(values = c("darkgrey", "black"))+ theme(legend.position="none")

  return(p)
}

ScidrForestPlotSpannerLabels <- function(gtsummary_table, n_tables) {
  out <- rep(NA_character_, n_tables)
  spanning_header <- gtsummary_table$table_styling$spanning_header
  if (is.null(spanning_header) ||
      !all(c("column", "spanning_header") %in% names(spanning_header))) {
    return(out)
  }

  table_index <- suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", spanning_header$column)))
  valid_rows <- !is.na(table_index) &
    table_index >= 1 &
    table_index <= n_tables &
    !is.na(spanning_header$spanning_header) &
    spanning_header$spanning_header != ""

  if (!any(valid_rows)) {
    return(out)
  }

  label_data <- data.frame(
    TableIndex = table_index[valid_rows],
    Label = as.character(spanning_header$spanning_header[valid_rows]),
    stringsAsFactors = FALSE
  )
  label_data <- label_data[!duplicated(label_data$TableIndex), , drop = FALSE]
  out[label_data$TableIndex] <- label_data$Label
  out
}
