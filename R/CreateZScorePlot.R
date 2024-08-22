#' Create a Z-Score Plot with Statistical Significance
#'
#' This function generates a Z-score plot to compare multiple variables across
#' different groups. It offers options for parametric or non-parametric tests,
#' ordinal variable conversion, and custom labeling. Significant p-values and
#' FDR-adjusted p-values are highlighted on the plot.
#'
#'
#' @param Data A dataframe containing the data to be analyzed.
#' @param TargetVar A string specifying the column name of the grouping variable.
#' @param Variables A vector of strings specifying the column names of the variables to be analyzed.
#' @param VariableCategories An optional vector categorizing the variables.
#' @param Relabel Logical; if TRUE, variables will be relabeled using their labels from the dataframe.
#' @param sort Logical; if TRUE, variables will be sorted by category and p-value.
#' @param RemoveXAxisLabels Logical; if TRUE, X-axis labels will be removed.
#' @param Ordinal Logical; if TRUE, ordinal variables will be converted to numeric.
#' @param Parametric Logical; if TRUE, parametric tests (t-test/ANOVA) will be used; otherwise, non-parametric tests (Wilcoxon/Kruskal-Wallis) will be used.
#' @param SigP_YCoord Numeric; the y-coordinate for marking significant p-values.
#' @param SigFDR_YCoord Numeric; the y-coordinate for marking significant FDR-adjusted p-values.
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
                              Relabel = TRUE, sort = TRUE, RemoveXAxisLabels = TRUE, Ordinal = TRUE,
                              Parametric = TRUE, SigP_YCoord = 1.5, SigFDR_YCoord = 1.6){
  if (Ordinal) {
    Data <- ConvertOrdinalToNumeric(Data, Variables)
  }

  classcolors <- c(paletteer::paletteer_d("calecopal::superbloom2"),
                   paletteer::paletteer_d("calecopal::vermillion"), paletteer::paletteer_d("fishualize::Antennarius_commerson"),
                   paletteer::paletteer_d("fishualize::Bodianus_rufus"))
  scaledData <- Data[c(TargetVar, Variables)]
  scaledData[Variables] <- scale(scaledData[Variables])
  colnames(scaledData)[1] <- "Group"
  n_groups <- length(unique(scaledData$Group))
  melted <- tidyr::pivot_longer(scaledData, cols = all_of(Variables),
                                names_to = "variable", values_to = "value")
  melted$variable <- factor(melted$variable, levels = unique(melted$variable))

  # Add contingency for when one group doesn't have much data
  melted$variable <- factor(melted$variable, levels = unique(melted$variable))
  melted <- melted[!is.na(melted$value),]
  melted <- melted %>% group_by(Group, variable)%>% filter(n() > 2) %>% ungroup() %>%
    group_by(variable) %>%
    filter(n_distinct(Group) > 1) %>% filter(n() > 2) %>%# filter out variables with few observations
    ungroup()

  # Perform statistical tests based on the number of groups and Parametric flag
  if (Parametric) {
    if (n_groups == 2) {
      stat.test <- melted %>%
        dplyr::group_by(variable) %>%
        rstatix::t_test(value ~ Group, var.equal = TRUE) %>%
        rstatix::adjust_pvalue(method = "BH") %>%
        rstatix::add_significance()
    } else {
      stat.test <- melted %>%
        dplyr::group_by(variable) %>%
        rstatix::anova_test(value ~ Group) %>%
        rstatix::adjust_pvalue(method = "BH") %>%
        rstatix::add_significance()
    }
  } else {
    # Non-parametric tests
    if (n_groups == 2) {
      stat.test <- melted %>%
        dplyr::group_by(variable) %>%
        rstatix::wilcox_test(value ~ Group, var.equal = TRUE) %>%
        rstatix::adjust_pvalue(method = "BH") %>%
        rstatix::add_significance()
    } else {
      stat.test <- melted %>%
        dplyr::group_by(variable) %>%
        rstatix::kruskal_test(value ~ Group) %>%
        rstatix::adjust_pvalue(method = "BH") %>%
        rstatix::add_significance()
    }
  }


  GroupMeans <- melted %>% dplyr::group_by(variable, Group) %>%
    dplyr::summarize(mean = mean(value, na.rm = TRUE), stderror = plotrix::std.error(value),
                     std = sd(value), n = n(), .groups = "drop")
  GroupMeans$pval <- stat.test$p[match(GroupMeans$variable,
                                       stat.test$variable)]
  GroupMeans$FDRpval <- stat.test$p.adj[match(GroupMeans$variable,
                                              stat.test$variable)]
  if (Relabel) {
    DataMergeTable <- data.frame(variable = colnames(Data),
                                 Label = sjlabelled::get_label(Data, def.value = colnames(Data)) %>%
                                   unname())
    GroupMeans <- left_join(GroupMeans, DataMergeTable, by = "variable")
  }else {
    GroupMeans$Label <- GroupMeans$variable
  }
  if (is.null(VariableCategories)) {
    GroupMeans$Category <- NA
  }else {
    GroupMeans$Category <- VariableCategories[match(GroupMeans$variable,
                                                    Variables)]
  }
  GroupMeans <- as.data.frame(GroupMeans)
  if (sort) {
    GroupMeans <- GroupMeans[order(GroupMeans$Category, GroupMeans$pval), ]
    GroupMeans$variable <- factor(GroupMeans$variable, levels = unique(GroupMeans$variable))
    GroupMeans$Label <- factor(GroupMeans$Label, levels = unique(GroupMeans$Label))
  }else {
    GroupMeans$variable <- factor(GroupMeans$variable, levels = unique(GroupMeans$variable))
    GroupMeans$Label <- factor(GroupMeans$Label, levels = unique(GroupMeans$Label))
  }
  pvaldata <- data.frame(variable = levels(GroupMeans$variable),
                         Label <- levels(GroupMeans$Label))
  pvaldata <- dplyr::right_join(pvaldata, stat.test, by = "variable")
  pvaldata$pvalline <- ifelse(pvaldata$p < 0.05, SigP_YCoord, NaN)
  pvaldata$FDRline <- ifelse(pvaldata$p.adj < 0.05, SigFDR_YCoord, NaN)
  pvaldata$variable <- factor(pvaldata$variable, levels = levels(GroupMeans$variable))
  pvaldata$Category <- ifelse(is.null(VariableCategories),
                              NA, VariableCategories)
  GroupMeans$Category[is.na(GroupMeans$Category)] <- "Unknown"
  GroupMeans <- GroupMeans %>% dplyr::mutate(Text = paste0("</br> Variable: ",
                                                           variable, "</br> Label: ", Label, "</br> p-value: ",
                                                           round(pval, 4), "</br> Group: ", Group, "</br> FDR: ",
                                                           round(FDRpval, 4), "</br> Category: ", Category))
  pvaldata <- pvaldata %>% dplyr::mutate(Text = paste0("</br> Variable: ",
                                                       variable, "</br> Label: ", Label, "</br> p-value: ",
                                                       round(p, 4), "</br> FDR: ", round(p.adj, 4), "</br> Category: ",
                                                       Category))
  pZ <- GroupMeans %>% ggplot(aes(x = Label, text = Text)) +
    geom_point(aes(y = mean, shape = Group, color = Category)) +
    geom_errorbar(aes(ymin = mean - stderror, ymax = mean +
                        stderror, color = Category), alpha = 0.5) + theme_minimal() +
    geom_point(data = pvaldata, aes(y = pvalline), color = "blue") +
    geom_point(data = pvaldata, aes(y = FDRline), color = "green") +
    guides(shape = FALSE, linetype = FALSE) + scale_y_continuous(limits = c(-2,
                                                                            2)) + ylab("Z-Score") + scale_color_manual(values = classcolors)
  if (RemoveXAxisLabels) {
    pZ <- pZ + theme(axis.text.x = element_blank())
  }
  else {
    pZ <- pZ + theme(axis.text.x = element_text(angle = 45,
                                                vjust = 1, hjust = 1))
  }
  return(pZ)
}
