
## Todo: Add relabel options
## Todo: Add Covariate Option
#' Plot Categorical Interaction Effects Matrix
#'
#' This function calculates and visualizes the interaction effects between categorical variables.
#'
#' @param Data The dataset containing the variables.
#' @param xVars A character vector of the names of the x-axis categorical variables.
#' @param yVars A character vector of the names of the y-axis categorical variables. Defaults to NULL,
#'        in which case it takes the same values as xVars.
#' @param xVarLabels A character vector of labels for the x-axis variables. Defaults to NULL,
#'        in which case it takes the same values as xVars.
#' @param yVarLabels A character vector of labels for the y-axis variables. Defaults to NULL,
#'        in which case it takes the same values as yVars.
#' @param interVar The name of the interaction variable.
#' @return A list containing matrices of interaction coefficients, p-values, ggplot objects for visualizations,
#'         and tables of FDR-corrected p-values.
#' @export
PlotCatInteractionEffectsMatrix <- function(Data, xVars, yVars = NULL, xVarLabels = NULL, yVarLabels = NULL, interVar ) {
  if(is.null(yVars)){
    yVars = xVars
  }
  if(is.null(xVarLabels)){
    xVarLabels = xVars
  }
  if(is.null(yVarLabels)){
    yVarLabels = yVars
  }

  r_P <- r_C <- r_S <- matrix(0, nrow = length(xVars), ncol = length(yVars))

  for (i in seq_len(length(xVars))) {
    for (j in seq_len(length(yVars))) {
      xVar <- xVars[i]
      yVar <- yVars[j]

      if(is.character(Data[[xVar]])){
        Data[[xVar]] <- as.factor(Data[[xVar]])
      }

      if(is.character(Data[[yVar]])){
        Data[[yVar]] <- as.factor(Data[[yVar]])
      }
      hadError <- FALSE

      tryCatch({
        m <- lm(get(yVar) ~ get(xVar)*get(interVar), data = Data)
      }, error = function(err){
        hadError <<- TRUE
      }, interrupt = function(){
        hadError <<- FALSE
      })

      if(hadError){
        interC <- NA
        interP <- NA
        interS <- NA
      } else {
        m <- lm(get(yVar) ~ get(xVar)*get(interVar), data = Data)
        d <- summary(m)
        dd <- as.data.frame(d$coefficients)
        interC <- dd$Estimate %>% tail(1)
        interP <- dd$`Pr(>|t|)` %>% tail(1)
        interS <- sum(sign(dd$Estimate))
      }

      r_P[i,j] <- interP
      r_C[i,j] <- interC
      r_S[i,j] <- interS
    }
  }

  rownames(r_C) <- xVarLabels
  colnames(r_C) <- yVarLabels
  rownames(r_P) <- xVarLabels
  colnames(r_P) <- yVarLabels
  rownames(r_S) <- xVarLabels
  colnames(r_S) <- yVarLabels

  m_r_C <- r_C %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "C")
  m_r_P <- r_P %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "P")
  m_r_S <- r_S %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "S")

  m_G <- left_join(m_r_C, m_r_P) %>% left_join(m_r_S)

  m_G$sign <- sign(m_G$S) %>% factor(levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign[m_G$P > 0.05] <- "ns"
  m_G$sig <- gtools::stars.pval(m_G$P)
  m_G$sig[m_G$sig %in% c(".", "+", " ")] <- ""
  m_G$sig <- paste(m_G$sign, m_G$sig) %>% factor(levels = c("+ ***", "+ **", "+ *","ns ", "- *", "- **", "- ***"))

  ## FDR correction for r_P
  m_G$P_FDR <- p.adjust(m_G$P, method = "fdr")
  m_G$sign_FDR <- sign(m_G$S) %>% factor(levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign_FDR[m_G$P_FDR > 0.05] <- "ns"
  m_G$sig_FDR <- gtools::stars.pval(m_G$P_FDR)
  m_G$sig_FDR[m_G$sig_FDR %in% c(".", "+", " ")] <- ""
  m_G$sig_FDR <- paste(m_G$sign_FDR, m_G$sig_FDR) %>% factor(levels = c("+ ***", "+ **", "+ *","ns ", "- *", "- **", "- ***"))

  p <-  m_G %>%
    ggplot(aes(x = X, y = Y, fill = sig)) +
    geom_tile() +
    scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white", "lightblue2", "steelblue3", "blue")), drop = F) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "No Multiple Comparison Correction")

  p_FDR <-  m_G %>%
    ggplot(aes(x = X, y = Y, fill = sig_FDR)) +
    geom_tile() +
    scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white", "lightblue2", "steelblue3", "blue")), drop = F) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "FDR Correction")

  pvaltable_FDR <- m_G %>%
    data.frame() %>%
    select(X, Y, P_FDR) %>%
    pivot_wider(names_from = X, values_from = P_FDR)

  return(list(C = r_C, pvals = r_P, p = p, p_FDR = p_FDR, pvals_FDR = pvaltable_FDR))
}
