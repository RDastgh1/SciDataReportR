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
CreateZScorePlot <- function(Data, TargetVar, Variables, VariableCategories = NULL, Relabel = TRUE, sort = TRUE, RemoveXAxisLabels = TRUE, Ordinal = TRUE) {

  if(Ordinal){

    #Then Convert to Numeric

    Data <- ConvertOrdinalToNumeric(Data, Variables)
    #DataFrame[Variables] <- lapply(DataFrame[Variables], as.numeric)
  }


# Relabeling (if required)
if (Relabel) {
  TargetVar <- sjlabelled::get_label(Data[TargetVar], def.value = colnames(Data[TargetVar]))
  colnames(Data)<- sjlabelled::get_label(Data %>% select(all_of(Variables)),
                                         def.value = colnames(Data %>% select(all_of(Variables))))
}

# Define colors for different categories
classcolors <- c(
  paletteer::paletteer_d("calecopal::superbloom2"),
  paletteer::paletteer_d("calecopal::vermillion"),
  paletteer::paletteer_d("fishualize::Antennarius_commerson"),
  paletteer::paletteer_d("fishualize::Bodianus_rufus")
)

# Scale the data
scaledData <- Data[c(TargetVar, Variables)]
scaledData[Variables] <- scale(scaledData[Variables])
colnames(scaledData)[1] <- "Group"

# Determine the number of groups
n_groups <- length(unique(scaledData$Group))

# Melt the data for plotting
melted <- tidyr::pivot_longer(scaledData, cols = Variables, names_to = "variable", values_to = "value")
melted$variable <- factor(melted$variable, levels = unique(melted$variable))

# Perform statistical tests
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

# Calculate group means and errors
GroupMeans <- melted %>%
  dplyr::group_by(variable, Group) %>%
  dplyr::summarize(
    mean = mean(value, na.rm = TRUE),
    stderror = plotrix::std.error(value),
    std = sd(value),
    n = n(),
    .groups = 'drop'
  )

GroupMeans$pval <- stat.test$p[match(GroupMeans$variable, stat.test$variable)]
GroupMeans$FDRpval <- stat.test$p.adj[match(GroupMeans$variable, stat.test$variable)]
GroupMeans$Category <- ifelse(is.null(VariableCategories), NA, VariableCategories[match(GroupMeans$variable, Variables)])
GroupMeans <- as.data.frame(GroupMeans)

# Sort the data if required
if (sort) {
  GroupMeans <- GroupMeans[order(GroupMeans$Category, GroupMeans$pval),]
  GroupMeans$variable <- factor(GroupMeans$variable, levels = unique(GroupMeans$variable))
}

# Add p-values
pvaldata <- data.frame(variable = levels(GroupMeans$variable))
pvaldata <- dplyr::right_join(pvaldata, stat.test, by = "variable")
pvaldata$pvalline <- ifelse(pvaldata$p < 0.05, 1.5, NaN)
pvaldata$FDRline <- ifelse(pvaldata$p.adj < 0.05, 1.6, NaN)
pvaldata$variable <- factor(pvaldata$variable, levels = levels(GroupMeans$variable))
pvaldata$Category <- ifelse(is.null(VariableCategories), NA, VariableCategories)

# Add unknown category if missing
GroupMeans$Category[is.na(GroupMeans$Category)] <- "Unknown"

# Add text annotations
GroupMeans <- GroupMeans %>%
  dplyr::mutate(Text = paste0("</br> Variable: ", variable,
                              "</br> p-value: ", round(pval, 4),
                              "</br> Group: ", Group,
                              "</br> FDR: ", round(FDRpval, 4),
                              "</br> Category: ", Category))

pvaldata <- pvaldata %>%
  dplyr::mutate(Text = paste0("</br> Variable: ", variable,
                              "</br> p-value: ", round(p, 4),
                              "</br> FDR: ", round(p.adj, 4),
                              "</br> Category: ", Category))

# Create the plot
pZ <- GroupMeans %>%
  ggplot(aes(x = variable, text = Text)) +
  geom_point(aes(y = mean, shape = Group, color = Category)) +
  geom_errorbar(aes(ymin = mean - stderror, ymax = mean + stderror, color = Category), alpha = 0.5) +
  theme_minimal() +
  geom_point(data = pvaldata, aes(y = pvalline), color = "blue") +
  geom_point(data = pvaldata, aes(y = FDRline), color = "green") +
  guides(shape = FALSE, linetype = FALSE) +
  scale_y_continuous(limits = c(-2, 2)) +
  ylab("Z-Score") +
  scale_color_manual(values = classcolors)

# Remove X-axis labels if specified
if (RemoveXAxisLabels) {
  pZ <- pZ + theme(axis.text.x = element_blank())
} else {
  pZ <- pZ + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

return(pZ)

# Use the line below to convert to plotly for an interactive plot
# ggplotly(pZ, tooltip = "text")
}
