#' Plot and Summarize Group Statistics
#'
#' This function creates plots and summary tables to compare two groups within a dataset.
#'
#' @param Data A data frame containing the data to be analyzed.
#' @param Variables A character vector of variable names to be included in the analysis.
#' @param VariableCategories A data frame containing variable categories (optional).
#' @param impClust A character string representing the important cluster.
#' @param normalClust A character string representing the normal cluster.
#' @return A list containing a plot, a summary table, and a p-value table.
#' @export
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import VIM
#' @import fastDummies
#' @import gtsummary
#' @import gtools
#' @import paletteer
#' @import stringr

Plot2GroupStats <- function(Data, Variables, VariableCategories = NULL, impClust, normalClust) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(VIM)
  library(fastDummies)
  library(gtsummary)
  library(gtools)
  library(paletteer)
  library(stringr)

  # Filter and prepare data
  tData <- Data %>%
    select(Cluster, all_of(Variables)) %>%
    filter(Cluster %in% c(normalClust, impClust))
  tData$Cluster <- factor(tData$Cluster, levels = c(normalClust, impClust))

  # Remove variables with over 80% missing data
  aggr_plot_all <- VIM::aggr(tData, col = c('navyblue', 'red'), numbers = TRUE, sortVars = TRUE,
                             cex.axis = .7, gap = 3, ylab = c("Histogram of missing data", "Pattern"),
                             bars = TRUE, labels = TRUE, plot = FALSE)
  thresh <- 0.75 * nrow(tData)
  remove_ind <- aggr_plot_all$missings$Count > thresh
  removed_75pmissing <- colnames(tData[, remove_ind])
  tData <- select(tData, -all_of(removed_75pmissing))

  # Remove factor variables with less than 2 levels
  l <- tData %>%
    summarise_if(is.factor, funs(nlevels(factor(.)))) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[l < 2])))

  # Remove numeric variables with zero standard deviation
  l <- tData %>%
    summarise_if(is.numeric, funs(sd(., na.rm = TRUE))) %>% as.list()
  tData <- tData %>% select(-all_of(names(l[is.na(l) | l == 0])))

  # Generate summary table
  gtabp <- tData %>%
    tbl_summary(
      by = Cluster,
      type = list(where(is.numeric) ~ "continuous"),
      statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)")
    ) %>%
    add_n() %>%
    add_p(test = list(all_continuous() ~ "aov", all_categorical() ~ "chisq.test")) %>%
    bold_p() %>%
    modify_header(label = "**Cluster**") %>%
    bold_labels()

  # Remove factor variables with more than 10 levels
  pred <- function(x) is.factor(x) & nlevels(factor(x)) > 10
  varsToRemove <- colnames(tData[sapply(tData, pred)])
  tData <- tData %>% select(-all_of(varsToRemove))

  # Create dummy variables if necessary
  if (sum(lapply(tData, class) == "factor") > 1) {
    d <- fastDummies::dummy_cols(tData %>% select(-Cluster), remove_selected_columns = TRUE, remove_first_dummy = TRUE)
    d$Cluster <- tData$Cluster

    l <- d %>%
      summarise_if(is.factor, funs(nlevels(factor(.)))) %>% as.list()
    d <- d %>% select(-all_of(names(l[l < 2])))

    newVars <- colnames(d)[colnames(d) %notin% colnames(tData)]
    d[newVars] <- lapply(d[newVars], as.factor)

    if (!is.null(VariableCategories)) {
      for (var in newVars) {
        oldvar <- str_split(var, "_")[[1]][1]
        newline <- VariableCategories %>% filter(Variable == oldvar)
        newline$Variable <- var
        VariableCategories <- rbind(VariableCategories, newline)
      }
    }
  } else {
    d <- tData
  }

  # Generate dummy summary table
  gtabd <- d %>%
    tbl_summary(
      by = Cluster,
      type = list(where(is.numeric) ~ "continuous", where(is.factor) ~ "dichotomous"),
      statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)", all_dichotomous() ~ "{p}%")
    ) %>%
    add_n() %>%
    add_difference()

  # Prepare p-value table
  pvaltable <- gtabd$table_body %>% filter(row_type == "label")
  pvaltable <- pvaltable[!is.na(pvaltable$p.value), ]
  pvaltable$Sig <- gtools::stars.pval(pvaltable$p.value)
  pvaltable$p.adj <- p.adjust(pvaltable$p.value, method = "fdr")
  pvaltable$logp <- log10(pvaltable$p.value) * sign(pvaltable$estimate)

  # Create plot
  if (is.null(VariableCategories)) {
    p <- pvaltable %>% ggplot(aes(y = variable, x = logp, color = Sig, shape = Sig))
  } else {
    pvaltable$Category <- VariableCategories$Category[match(pvaltable$variable, VariableCategories$Variable)]
    p <- pvaltable %>% ggplot(aes(y = variable, x = logp, color = Category, shape = Sig)) +
      xlab("-log10(pval)") +
      geom_vline(xintercept = 0) +
      theme(plot.title = element_text(hjust = 1))
  }

  p <- p +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.5) +
    geom_point() +
    geom_text_repel(data = subset(pvaltable, `p.adj` < 0.05), aes(label = variable)) +
    scale_y_discrete(limits = rev(pvaltable$variable)) +
    paletteer::scale_color_paletteer_d("fishualize::Scarus_quoyi") +
    xlab("-log10(pval)") +
    geom_vline(xintercept = 0) +
    theme_bw() +
    ggtitle(paste("---> Higher in", impClust)) +
    theme(plot.title = element_text(hjust = 1))

  return(list(p = p, ptable = gtabp, pvaltable = pvaltable))
}
