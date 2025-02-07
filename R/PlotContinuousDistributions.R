#' Plot Continuous Distributions
#'
#' This function creates raincloud plots to visualize the distributions of continuous variables in a data frame.
#'
#' @param DataFrame A data frame containing the variables to be plotted.
#' @param Variables A character vector specifying the names of the variables to be plotted. If NULL, all numeric variables in DataFrame will be used.
#' @param Relabel Logical, whether to relabel the variables using their labels if available. Default is TRUE.
#' @param ncol Number of columns in the facet grid. Default is 3.
#' @param Ordinal Logical, indicating whether ordinal variables should be included.
#' @return A raincloud plot visualizing the distributions of continuous variables.
#'
#' @importFrom dplyr mutate pivot_longer select
#' @importFrom ggplot2 ggplot theme_bw guides coord_flip facet_wrap scale_fill_manual element_blank
#' @suggests ggrain, sjlabelled
#' @export

PlotContinuousDistributions <- function (DataFrame, Variables = NULL, Relabel = TRUE, ncol = 3, Ordinal = T) {
  # If Variables argument is NULL, use all numeric variables
  if (is.null(Variables)) {
    Variables <- getNumVars(DataFrame, Ordinal = F)
    if (Ordinal) {
      Variables <- getNumVars(DataFrame, Ordinal = T)
    }
  }
  if (Ordinal) {
    OriginalLabels <- sjlabelled::get_label(DataFrame, def.value = colnames(DataFrame))
    DataFrame <- ConvertOrdinalToNumeric(DataFrame, Variables)
    DataFrame[Variables] <- lapply(DataFrame[Variables], as.numeric)
    for (var in Variables) {
      DataFrame[[var]] <- set_label(DataFrame[[var]], OriginalLabels[[var]])
    }
  }
  if (Relabel == TRUE) {
    facetlabels <- createFacetLabels(DataFrame %>% select(all_of(Variables)))
  } else {
    facetlabels <- Variables
  }
  ContData <- DataFrame %>% pivot_longer(cols = all_of(Variables)) %>%
    group_by(name) %>% mutate(Mean = mean(value, na.rm = TRUE)) %>%
    ungroup()
  ContData$name <- factor(ContData$name, levels = Variables, labels = facetlabels)

  p <- ggplot(ContData, aes(y = value, x = 1, fill = "1")) +
    ggrain::geom_rain(alpha = 0.5) + theme_bw() + guides(fill = "none", color = "none") +
    coord_flip() + facet_wrap(~name, scales = "free", ncol = ncol) +
    scale_fill_manual(values = c("#6EC259")) +
    theme(axis.title.y = element_blank())

  return(p)
}
