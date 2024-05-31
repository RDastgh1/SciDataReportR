#' Plot Associations
#'
#' This function generates scatter plots or box plots to visualize the relationship between two variables.
#'
#' @param DataFrame The data frame containing the variables of interest.
#' @param Var1 The name of the first variable.
#' @param Var2 The name of the second variable.
#' @return A ggplot object representing the relationship between the variables.
#' @import ggplot2 ggstatsplot dplyr
#' @export
PlotAssociations <- function (DataFrame, Var1, Var2)
{
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
    s <- p$labels$subtitle %>% as.character()
    pval <- parse_number(get_pval(s))
    rval <- parse_number(get_rval(s))
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
      pval <- parse_number(get_pval(s))
      rval <- parse_number(get_rval(s))
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
      pval <- parse_number(get_pval(s))
      if (pval < 0.05) {
        p <- p + theme(panel.background = element_rect(fill = "lightblue"))
      }
    }
  }
  p <- p + theme(text = element_text(size = 10))
  return(p)
}
