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
  # If Variables is not specified, include all columns except GroupVariable
  if (is.null(Variables)) {
    Variables <- Data %>% select(-all_of(GroupVariable)) %>% colnames()
  }

  # Create a new column for the group variable
  Data$GroupVariablex <- Data[[GroupVariable]]

  # Remove the original group variable column
  Data <- Data %>% select(-all_of(GroupVariable))

  # Remove zero-length variables
  tData <- Data %>% select(GroupVariablex, all_of(Variables))
  l <- tData %>%
    summarise_if(is.factor, funs(nlevels(factor(.)))) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[l < 2])))
  l <- tData %>%
    summarise_if(is.numeric, funs(sd(., na.rm = T))) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[is.na(l) | l == 0])))

  # Perform statistical tests
  tab1 <- arsenal::tableby(GroupVariablex ~ ., data = tData)
  pvaltable <- arsenal::tests(tab1)
  pvaltable$logp <- -log10(pvaltable$p.value)

  # Add significance stars
  pvaltable$Sig <- gtools::stars.pval(pvaltable$p.value)
  pvaltable$Sig[pvaltable$Sig == "." | pvaltable$Sig == "+" | pvaltable$Sig == " "] <- "ns"
  pvaltable$Sig <- factor(pvaltable$Sig, levels = c("ns", "*", "**", "***", "****", "*****"))

  # Adjust p-values for multiple comparisons
  pvaltable$p.adj <- p.adjust(pvaltable$p.value, method = "fdr")

  # Relabel variables if necessary
  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    VarLabels <- sjlabelled::get_label(Data[pvaltable$Variable], def.value = colnames(Data[pvaltable$Variable])) %>%
      as.data.frame() %>% rownames_to_column()
    pvaltable$VarLabel <- VarLabels$.
  } else {
    pvaltable$VarLabel <- pvaltable$Variable
  }
  pvaltable$VarLabel <- factor(pvaltable$VarLabel, levels = unique(pvaltable$VarLabel))

  # Create the plot
  if (is.null(VariableCategories)) {
    p <- pvaltable %>% ggplot(aes(y = VarLabel, x = logp, shape = Sig))
  } else {
    pvaltable$Category <- VariableCategories[match(pvaltable$Variable, Variables)]
    p <- pvaltable %>% ggplot(aes(y = Variable, x = logp, color = Category, shape = Sig))
  }
  p <- p + geom_point() +
    ggrepel::geom_text_repel(data = subset(pvaltable, `p.adj` < 0.1),
                             aes(label = VarLabel)) +
    scale_y_discrete(limits = rev(pvaltable$VarLabel)) +
    xlab("-log10(pval)") +
    scale_shape_manual(values = c(21, 16, 17, 15, 18, 8), breaks = c("ns", "*", "**", "***", "****"), drop = F)

  return(p)
}
