#' Plot P-Value Comparisons
#'
#' This function generates a plot comparing p-values for different variables
#' across two or more groups.
#'
#' @param Data Data frame containing the variables to compare.
#' @param GroupVariable Character string specifying the name of the column in
#'   `Data` that contains the group labels.
#' @param Variables Character vector specifying the names of the columns in
#'   `Data` to include in the comparison. If `NULL`, all columns except
#'   `GroupVariable` are included.
#' @param VariableCategories Character vector specifying the categories for
#'   each variable. If `NULL`, no categories are used.
#' @param Relabel Logical indicating whether to replace missing labels with
#'   the column names.
#'
#' @return A ggplot object displaying the p-value comparisons.
#'
#' @import arsenal
#' @import ggplot2
#' @import ggrepel
#' @import gtools
#' @import sjlabelled
#' @import dplyr
#' @export
PlotPValueComparisons <- function(Data, GroupVariable, Variables = NULL, VariableCategories = NULL, Relabel = TRUE) {
  # Validate and prepare Variables
  if (is.null(Variables)) {
    Variables <- Data %>% select(-all_of(GroupVariable)) %>% colnames()
  }

  # Preserve original GroupVariable and create a temporary column
  Data$GroupVariablex <- Data[[GroupVariable]]
  Data <- Data %>% select(-all_of(GroupVariable))

  # Select relevant variables
  tData <- Data %>% select(GroupVariablex, all_of(Variables))

  # Remove factor variables with zero levels
  l <- tData %>% summarise_if(is.factor, ~ nlevels(factor(.))) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[l < 2])))

  # Remove numeric variables with zero or missing standard deviation
  l <- tData %>% summarise_if(is.numeric, ~ sd(., na.rm = TRUE)) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[is.na(l) | l == 0])))

  # Perform statistical tests using arsenal::tableby
  tab1 <- arsenal::tableby(GroupVariablex ~ ., data = tData)
  pvaltable <- arsenal::tests(tab1)

  # Calculate log-transformed p-values and significance levels
  pvaltable$logp <- -log10(pvaltable$p.value)
  pvaltable$Sig <- gtools::stars.pval(pvaltable$p.value)
  pvaltable$Sig[pvaltable$Sig %in% c(".", "+", " ")] <- "ns"
  pvaltable$Sig <- factor(pvaltable$Sig, levels = c("ns", "*", "**", "***", "****", "*****"))

  # Adjust p-values for multiple testing
  pvaltable$p.adj <- p.adjust(pvaltable$p.value, method = "fdr")

  # Optionally relabel variables
  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    VarLabels <- sjlabelled::get_label(Data[as.character(pvaltable$Variable)],
                                       def.value = pvaltable$Variable) %>%
      as.data.frame() %>% rownames_to_column()
    pvaltable$VarLabel <- VarLabels$.
  } else {
    pvaltable$VarLabel <- pvaltable$Variable
  }
  pvaltable$VarLabel <- factor(pvaltable$VarLabel, levels = unique(pvaltable$VarLabel))

  # Generate the plot
  if (is.null(VariableCategories)) {
    p <- pvaltable %>% ggplot(aes(y = Variable, x = logp, shape = Sig))
  } else {
    pvaltable$Category <- VariableCategories[match(pvaltable$Variable, Variables)]
    p <- pvaltable %>% ggplot(aes(y = Variable, x = logp, color = Category, shape = Sig))
  }

  # Add points, text labels for significant variables, and customize scales
  p <- p +
    geom_point() +
    ggrepel::geom_text_repel(data = subset(pvaltable, p.adj < 0.1),
                             aes(label = Variable)) +
    scale_y_discrete(limits = rev(pvaltable$Variable)) +
    xlab("-log10(p-value)") +
    scale_shape_manual(values = c(21, 16, 17, 15, 18, 8),
                       breaks = c("ns", "*", "**", "***", "****"), drop = FALSE)

  return(p)
}
