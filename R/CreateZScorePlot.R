#' Create a Z-Score Plot with Statistical Significance
#'
#' This function generates a Z-score plot for the provided data, highlighting the mean values, standard errors, and p-values for significance testing between groups.
#'
#' @param Data A data frame containing the target variable and variables of interest.
#' @param TargetVar A string specifying the name of the target variable (grouping variable).
#' @param Variables A vector of strings specifying the names of the variables to be plotted.
#' @param VariableCategories An optional vector of categories corresponding to the variables.
#' @param Relabel A logical value indicating whether to relabel the variables. Default is TRUE.
#' @param sort A logical value indicating whether to sort the variables by p-value and category. Default is TRUE.
#' @param RemoveXAxisLabels A logical value indicating whether to remove X-axis labels. Default is TRUE.
#' @return A ggplot object representing the Z-score plot.
#' @import dplyr
#' @import ggplot2
#' @import plotly
#' @import rstatix
#' @import paletteer
#' @import plotrix
#' @importFrom tidyr pivot_longer
#' @export
CreateZScorePlot <- function (Data, TargetVar, Variables, VariableCategories = NULL,
                              Relabel = TRUE, sort = TRUE, RemoveXAxisLabels = TRUE, Ordinal = TRUE)
  roupMeans <- as.data.frame(GroupMeans)
if (sort) {
  GroupMeans <- GroupMeans[order(GroupMeans$Category, GroupMeans$pval),
  ]
  GroupMeans$variable <- factor(GroupMeans$variable, levels = unique(GroupMeans$variable))
  GroupMeans$Label <- factor(GroupMeans$Label, levels = unique(GroupMeans$Label))
}
