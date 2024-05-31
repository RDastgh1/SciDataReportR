#' Plot Categorical Distributions
#'
#' This function creates plots to visualize the distributions of categorical variables in a dataframe.
#'
#' @param DataFrame The dataframe containing the variables to be plotted.
#' @param Variables Optional. A character vector specifying the names of the categorical variables to be plotted. If NULL, categorical variables are automatically detected.
#' @param Relabel Logical. If TRUE, missing labels in the dataframe are replaced with column names as labels for plotting.
#' @return A plot visualizing the distributions of categorical variables.
#' @importFrom sjlabelled get_label
#' @importFrom dplyr select
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom magrittr "%>%"
#' @export

PlotCategoricalDistributions <- function(DataFrame, Variables = NULL, Relabel = TRUE, Ordinal = TRUE){
  if (is.null(Variables)) {
    Variables <- getCatVars(DataFrame)
  }

  if (Relabel == TRUE) {
    DataFrame <- ReplaceMissingLabels(DataFrame)
    facetlabels <- sjlabelled::get_label(DataFrame %>% dplyr::select(all_of(Variables)))
    names(facetlabels) <- NULL
  } else {
    facetlabels <- Variables
  }

  df <- inspectdf::inspect_cat(DataFrame %>% dplyr::select(all_of(Variables)), include_int = TRUE)
  d <- df %>%
    inspectdf::show_plot(col_palette = 1)

  d <- d + scale_x_discrete(limits = Variables, labels = facetlabels)

  return(d)
}
