#' Plot Associations
#'
#' This function generates scatter plots or box plots to visualize the relationship between two variables.
#'
#' @param data The data frame containing the variables of interest.
#' @param Var1 The name of the first variable.
#' @param Var2 The name of the second variable.
#' @param Ordinal Logical, indicating whether ordinal variables should be included.
#' @return A ggplot object representing the relationship between the variables.
#' @import ggplot2 ggstatsplot dplyr
#'
#' @references
#' This function wraps \pkg{ggstatsplot}. Please cite:
#'
#' Patil, I. (2021). Visualizations with statistical details: The 'ggstatsplot'
#' approach. \emph{Journal of Open Source Software}, 6(61), 3167.
#' \doi{10.21105/joss.03167}
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @examples
#' data(SampleData)
#'
#' # Two categorical variables (grouped bar chart)
#' PlotAssociations(SampleData, "Diagnosis", "Genotype")
#'
#' # Two continuous variables (scatter plot with correlation)
#' PlotAssociations(SampleData, "age", "AXL")
#'
#' # One continuous and one categorical variable (box/violin plot)
#' PlotAssociations(SampleData, "Diagnosis", "AXL")
#' @export
PlotAssociations <- function(data,
    Var1,
    Var2,
    Ordinal = FALSE,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "PlotAssociations(DataFrame)", "PlotAssociations(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  TestFrame <- na.omit(DataFrame[c(Var1, Var2)])
  if (nrow(TestFrame) == 0) {
    return(ggplot(DataFrame, aes_string(x = Var1, y = Var2)))
  }
  T1 <- class(DataFrame[[Var1]])
  T2 <- class(DataFrame[[Var2]])
  if (length(DataFrame[[Var1]] %>% na.omit) == 0) {
    T1 <- "logical"
  }
  if (length(DataFrame[[Var2]] %>% na.omit) == 0) {
    T2 <- "logical"
  }
  if ((T1 == "numeric") & (T2 == "numeric")) {
    type <- "NumNum"
  }else if (T1 %in% c("character", "logical", "factor") & T2 %in%
            c("character", "logical", "factor")) {
    type <- "CatCat"
  }else {
    type <- "NumCat"
    if (T1 == "numeric") {
      NumVar <- Var1
      CatVar <- Var2
    } else {
      NumVar <- Var2
      CatVar <- Var1
    }
  }
  if (type == "NumNum") {
    p <- ggscatterstats(data = DataFrame, x = !!Var1, y = !!Var2,
                        bf.message = FALSE)
    pval <- extract_stats(p)[["subtitle_data"]][["p.value"]][1]
    rval <- extract_stats(p)[["subtitle_data"]][["estimate"]][1]
    if (is.na(pval)) {
      p <- p + theme(panel.background = element_rect(fill = "#ff6347"))
    }
    else {
      if (pval < 0.05) {
        p <- p + theme(panel.background = element_rect(fill = "lightblue"))
      }
      if (rval == 1) {
        p <- p + theme(panel.background = element_rect(fill = "#D8BFD8"))
      }
    }
  }
  else if (type == "CatCat") {
    DataFrame[[Var1]] <- DataFrame[[Var1]] %>% addNA()
    DataFrame[[Var2]] <- DataFrame[[Var2]] %>% addNA()
    p <- ggbarstats(data = DataFrame, x = !!Var1, y = !!Var2,
                    bf.message = FALSE, label = "both")
    s <- p$labels$subtitle %>% as.character()
    if (length(s) == 0) {
      p <- p + theme(panel.background = element_rect(fill = "#ff6347"))
    }
    else {
      pval <- extract_stats(p)[["subtitle_data"]][["p.value"]][1]
      rval <- extract_stats(p)[["subtitle_data"]][["estimate"]][1]
      if (pval < 0.05) {
        p <- p + theme(panel.background = element_rect(fill = "lightblue"))
      }
      if (rval == 1) {
        p <- p + theme(panel.background = element_rect(fill = "#D8BFD8"))
      }
    }
  }
  else if (type == "NumCat") {
    DataFrame <- DataFrame[!is.na(DataFrame[[NumVar]]), ]
    DataFrame[[CatVar]] <- DataFrame[[CatVar]] %>% as.factor()
    p <- ggstatsplot::ggbetweenstats(data = DataFrame, x = !!CatVar, y = !!NumVar,
                                     plot.type = "box", bf.message = FALSE, pairwise.comparisons = FALSE)
    s <- p$labels$subtitle %>% as.character()
    if (length(s) == 0) {
      p <- p + theme(panel.background = element_rect(fill = "#ff6347"))
    }
    else {
      pval <- extract_stats(p)[["subtitle_data"]][["p.value"]][1]
      if (pval < 0.05) {
        p <- p + theme(panel.background = element_rect(fill = "lightblue"))
      }
    }
  }
  p <- p + theme(text = element_text(size = 10))
  return(p)
}
