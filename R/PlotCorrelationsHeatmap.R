#' Plot Correlations Heatmap
#'
#' This function calculates correlations between variables and plots them as a heatmap.
#'
#' @param Data The dataset containing the variables.
#' @param xVars A character vector of the names of the x-axis variables.
#' @param yVars A character vector of the names of the y-axis variables. Defaults to NULL.
#' @param covars A character vector of the names of covariate variables. Defaults to NULL.
#' @param FS The font size for the plot. Defaults to 3.
#' @param method The correlation method. Defaults to "pearson".
#' @param Relabel A logical indicating whether to relabel variables. Defaults to TRUE.
#' @param Ordinal Logical, indicating whether ordinal variables should be included.
#' @return A list containing matrices, ggplot objects for visualizations, and details of the method used.
#' @export

PlotCorrelationsHeatmap <- function(Data, xVars = NULL, yVars = NULL, covars = NULL, method = "pearson", Relabel = TRUE, Ordinal = FALSE) {
  removediag <- FALSE
  if (is.null(yVars)) {
    yVars <- xVars
    removediag <- TRUE
  }

  Variables <- c(xVars, yVars)
  Data <- ReplaceMissingLabels(Data)

  if (is.null(xVars)) {
    xVars <- getNumVars(Data, Ordinal = F)

    if(Ordinal){
      xVars <- getNumVars(Data, Ordinal = T)
    }
  }
  # Then Convert to Numeric
  if(Ordinal){
    Variables <- c(xVars, yVars)
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)

    # if any of the variables in x or y are ordinal, convert them to numeric
    #allvars <- c(xVars, yVars)
    #Data <- ConvertOrdinalToNumeric(Data, allvars)
  }











  if (length(covars) > 0) {
    char_vars <- sapply(Data, is.character)
    Data[char_vars] <- lapply(Data[char_vars], as.factor)
    factor_vars <- sapply(Data, is.factor)
    Data[factor_vars] <- lapply(Data[factor_vars], as.numeric)
    num_vars <- sapply(Data, is.numeric)
    cdata <- Data[covars]
  }else{
    cdata<- data.frame()
  }



  cdata <- scale(cdata)

  xdata <- Data[xVars]
  ydata <- Data[yVars]

  MC <- matrix(nrow = ncol(xdata), ncol = ncol(ydata))
  rownames(MC) <- colnames(xdata)
  colnames(MC) <- colnames(ydata)
  MP <- MC
  MN <- MC

  for (xi in 1:ncol(xdata)) {
    for (yi in 1:ncol(ydata)) {
      x <- xdata[, xi]
      y <- ydata[, yi]
      if (nrow(cdata) > 0) {
        d <- cbind(x, y, cdata)
      } else {
        d <- cbind(x, y)
      }

      d <- na.omit(d)
      hadError <- FALSE

      MN[xi, yi] <- nrow(d)

      #default output if we never successfully compute
      o <- list(estimate = NA_real_, p.value = NA_real_)

      # only try the test if we have ≥3 rows
      if (nrow(d) >= 3) {
        tryCatch({
          if (nrow(cdata) > 0) {
            tmp <- ppcor::pcor.test(d[,1], d[,2],
                                    d[, seq(3, 2 + length(covars))],
                                    method = method)
            o$estimate <- unname(tmp$estimate)
            o$p.value  <- tmp$p.value
          } else {
            tmp <- stats::cor.test(d[,1], d[,2], method = method)
            o$estimate <- unname(tmp$estimate)
            o$p.value  <- tmp$p.value
          }
        }, error = function(e) {
          message("Correlation failed for ", xi, ",", yi, ": ", e$message)
          # o stays as NA
        })
      }


      if (hadError) {
        MC[xi, yi] <- NA
        MP[xi, yi] <- NA
      } else {
        MC[xi, yi] <- o$estimate
        MP[xi, yi] <- o$p.value
      }
    }
  }

  pvals <- MP

  if (removediag) {
    diag(pvals) <- NaN
    diag(MC) <- NaN
  }

  M <- list()
  M$r <- MC
  M$p <- pvals
  M$npairs <- MN


  # if npairs <2, set r and p to NA
  M$r[M$npairs<2]<- NaN
  M$p[M$npairs<2]<- NaN
  M$r[M$npairs<2]<- NaN
  MP[M$npairs<2]<- NaN


  M_FDR <- M
  padjusted <- p.adjust(MP, method = "fdr")
  M_FDR$p <- matrix(padjusted, ncol = length(yVars))

  pvals_FDR <- M_FDR$p
  #pvals_FDR[is.na(pvals_FDR)] <- 1


  if (removediag) {
    diag(M_FDR$p) <- NaN
    diag(M_FDR$r) <- NaN
  }

  colnames(M$r) <- yVars
  rownames(M$r) <- xVars
  colnames(M$p) <- yVars
  rownames(M$p) <- xVars

  colnames(M_FDR$r) <- yVars
  rownames(M_FDR$r) <- xVars
  colnames(M_FDR$p) <- yVars
  rownames(M_FDR$p) <- xVars

  plot.data_R <- pivot_longer(as.data.frame(M$r) %>% rownames_to_column(var = "XVar"), cols = colnames(M$r), names_to = "YVar", values_to = "R")
  plot.data_P <- pivot_longer(as.data.frame(M$p) %>% rownames_to_column(var = "XVar"), cols = colnames(M$r), names_to = "YVar", values_to = "P")
  plot.data_P_adj <- pivot_longer(as.data.frame(M_FDR$p) %>% rownames_to_column(var = "XVar"), cols = colnames(M$r), names_to = "YVar", values_to = "P_adj")
  plot.data_npairs <- pivot_longer(as.data.frame(M$npairs) %>% rownames_to_column(var = "XVar"), cols = colnames(M$npairs), names_to = "YVar", values_to = "nPairs")

  plot.data <- suppressMessages(left_join(plot.data_R, plot.data_P) %>%
                                  left_join(plot.data_P_adj) %>%
                                  left_join(plot.data_npairs))

  plot.data$stars <- cut(plot.data$P, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))
  plot.data$stars_FDR <- cut(plot.data$P_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[as.character(plot.data$XVar)], def.value = plot.data$XVar) %>%
      as.data.frame() %>% rownames_to_column()

    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(plot.data$YVar)], def.value = plot.data$YVar) %>%
      as.data.frame() %>% rownames_to_column()

    colnames(ylabels) <- c("Variable", "label")

    plot.data$XLabel <- xlabels$label
    plot.data$YLabel <- ylabels$label
  } else {
    plot.data$XLabel <- plot.data$XVar
    plot.data$YLabel <- plot.data$YVar
  }

  plot.data$XLabel <- factor(plot.data$XLabel, ordered = FALSE, levels = rev(unique(plot.data$XLabel)))
  plot.data$YLabel <- factor(plot.data$YLabel, ordered = FALSE, levels = unique(plot.data$YLabel))

  PlotText <- paste("YVar", plot.data$YVar,
                    "</br> XVAR: ", plot.data$XVar,
                    "</br> P-Value: ", plot.data$P, plot.data$stars,
                    "</br> FDR-corrected P: ", plot.data$P_adj, plot.data$stars_FDR)

  P <- plot.data %>%
    ggplot(aes(y = XLabel, x = YLabel, fill = R, text = PlotText)) +
    geom_tile() +
    geom_text(aes(label = stars), color = "black") +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  P_FDR <- plot.data %>%
    ggplot(aes(y = XLabel, x = YLabel, fill = R, text = PlotText)) +
    geom_tile() +
    geom_text(aes(label = stars_FDR), color = "black") +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  M$plot <- P
  M_FDR$plot <- P_FDR

  return(list(Unadjusted = M, FDRCorrected = M_FDR, method = method, Relabel = Relabel, Covariates = covars))
}
