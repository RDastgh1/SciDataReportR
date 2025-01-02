#' Plot Point-Biserial Correlations Between Binary and Continuous Variables
#'
#' This function calculates point-biserial correlations between binary categorical variables and continuous variables,
#' creates a summary table, and generates heatmap-style plots to visualize the correlations.
#'
#' @param Data A dataframe containing the data.
#' @param CatVars A vector of column names representing binary categorical variables.
#' @param ContVars A vector of column names representing continuous variables.
#' @param Covariates A vector of column names representing covariates (optional).
#' @param Relabel A logical value indicating whether to relabel variables using variable labels (default: TRUE).
#' @param Ordinal A logical value indicating whether to treat ordinal variables as numeric (default: TRUE).
#' @return A list containing two elements:
#'   - `Unadjusted`: A list with the unadjusted p-value table and corresponding plot.
#'   - `FDRCorrected`: A list with the FDR-adjusted p-value table and corresponding plot.
#'   - `method`: The statistical method used ("R_pb").
#'   - `Relabel`: Logical indicating whether relabeling was applied.
#'   - `Covariates`: The covariates used.
#' @export
PlotPointCorrelationsHeatmap <- function (Data, CatVars, ContVars, Covariates = NULL, Relabel = TRUE, Ordinal = TRUE) {
  # Create a subset of the data
  DataSubset <- Data[c(CatVars, ContVars, if (!is.null(Covariates)) Covariates)]

  DataSubset <- ReplaceMissingLabels(DataSubset)  # Replace missing labels in the data

  # Function to check if variables are binary
  checkBinaryVars <- function(vars) {
    sapply(Data[vars], function(col) length(unique(na.omit(col))) <= 2)
  }

  # Verify that all categorical variables are binary
  non_binary_vars <- c(CatVars)[!checkBinaryVars(c(CatVars))]
  if (length(non_binary_vars) > 0) {
    stop("The following variables are not binary (they do not have exactly two unique values): ",
         paste(non_binary_vars, collapse = ", "))
  }

  # Pivot the data for analysis
  mData <- pivot_longer(DataSubset, cols = all_of(CatVars),
                        names_to = "CategoricalVariable", values_to = "CategoricalValue")
  mData <- pivot_longer(mData, cols = all_of(ContVars), names_to = "ContinuousVariable",
                        values_to = "ContinuousValue")

  # Convert variables to factors
  mData$ContinuousVariable <- as.factor(mData$ContinuousVariable)
  mData$CategoricalValue <- as.factor(mData$CategoricalValue)

  # Remove missing values
  mData <- na.omit(mData)

  # Filter out categorical variables with only one level
  nlevels.test <- mData %>% group_by(ContinuousVariable, CategoricalVariable) %>%
    mutate(n = length(unique(CategoricalValue)))
  mData <- mData[nlevels.test$n > 1, ]

  # Filter out continuous variables with zero variance
  nvar.test <- mData %>% group_by(ContinuousVariable) %>% mutate(v = var(ContinuousValue))
  mData <- mData[nvar.test$v > 0, ]

  # Calculate point-biserial correlation and number of valid pairs
  stat.test <- mData %>%
    group_by(ContinuousVariable, CategoricalVariable) %>%
    summarise(
      nPairs = sum(!is.na(CategoricalValue) & !is.na(ContinuousValue)),
      correlation = cor(ContinuousValue, as.numeric(factor(CategoricalValue)), use = "complete.obs"),
      p_value = cor.test(ContinuousValue, as.numeric(factor(CategoricalValue)))$p.value,
      .groups = "drop"
    )

  # Adjust p-values for multiple comparisons
  stat.test$p.adj <- p.adjust(stat.test$p_value, method = "fdr")
  stat.test <- stat.test %>% rstatix::add_significance() %>%
    rstatix::add_significance("p_value", output.col = "p<.05")
  stat.test$`p<.05`[stat.test$`p<.05` == "ns"] <- ""
  stat.test$p.adj.signif[stat.test$p.adj.signif == "ns"] <- ""
  stat.test$`p<.05` <- factor(stat.test$`p<.05`, levels = c("ns", "*", "**", "***", "****"),
                              labels = c("ns", "*", "**", "***", "***"))
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif,
                                   levels = c("ns", "*", "**", "***", "****"),
                                   labels = c("ns", "*", "**", "***", "***"))
  stat.test$test <- "point biserial correlation"

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[stat.test$CategoricalVariable],
                                     def.value = colnames(Data[stat.test$CategoricalVariable])) %>%
      as.data.frame() %>% rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")
    ylabels <- sjlabelled::get_label(Data[stat.test$ContinuousVariable],
                                     def.value = colnames(Data[stat.test$ContinuousVariable])) %>%
      as.data.frame() %>% rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")
    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  PlotText <- paste("</br>Cat Var Label:", stat.test$XLabel,
                    "</br> Cat Var:", stat.test$CategoricalVariable, "</br> Cont Var Label:",
                    stat.test$YLabel, "</br> Cont Var:", stat.test$ContinuousVariable,
                    "</br> P-Value: ", stat.test$p_value, stat.test$`p<.05`, "</br> FDR-corrected P: ",
                    stat.test$p.adj, stat.test$p.adj.signif, "</br> nPairs:", stat.test$nPairs)

  # Create unadjusted plot
  p <- ggplot(stat.test, aes(y = YLabel, x = XLabel,
                             fill = correlation, text = PlotText)) +
    geom_tile() + geom_text(aes(label = `p<.05`),
                            color = "black") +
    scale_fill_gradient2(limits = c(-1, 1), low = scales::muted("#FFA500"),
                         high = scales::muted("#008080")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = expression(r[pb]))

  # Create FDR-corrected plot
  p_FDR <- ggplot(stat.test, aes(y = YLabel, x = XLabel,
                                 fill = correlation, text = PlotText)) +
    geom_tile() + geom_text(aes(label = `p.adj.signif`),
                            color = "black") +
    scale_fill_gradient2(limits = c(-1, 1), low = scales::muted("#FFA500"),
                         high = scales::muted("#008080")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = expression(r[pb]))

  # Return results
  M <- list(PvalTable = stat.test, plot = p)
  M_FDR <- list(PvalTable = stat.test, plot = p_FDR)
  return(list(Unadjusted = M, FDRCorrected = M_FDR, method = "R_pb",
              Relabel = Relabel, Covariates = Covariates))
}
