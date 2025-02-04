#' Plot Phi Heatmap
#'
#' This function generates a heatmap of the Phi coefficient for the association between binary variables
#' in the provided dataset. The Phi coefficient measures the association between two binary variables,
#' and the function also performs chi-square tests for each pair of variables.
#'
#' @param Data A dataframe containing the variables of interest.
#' @param xVars A vector of variable names for the x-axis (predictor variables).
#' @param yVars A vector of variable names for the y-axis (outcome variables).
#' @param Relabel A logical indicating whether to relabel the variables with their labels (default is TRUE).
#' @param Ordinal A logical indicating whether to consider ordinal variables when determining binary variables (default is TRUE).
#'
#' @return A list containing the following components:
#'   \item{Unadjusted}{A list of results containing the Phi coefficients, p-values, and number of pairs for each variable combination.}
#'   \item{FDRCorrected}{A list containing FDR-corrected p-values.}
#'   \item{method}{A string indicating that the method used is "Phi".}
#'   \item{Relabel}{The relabel parameter.}
#'   \item{Covariates}{Currently set to NULL, as covariates are not implemented in this function.
#' @export
PlotPhiHeatmap <- function(Data, xVars = NULL, yVars = NULL, Relabel = TRUE, Ordinal = TRUE) {
  # Check if all variables in xVars and yVars are binary (i.e., have two or fewer unique values)

  checkBinaryVars <- function(vars) {
    sapply(Data[vars], function(col) length(unique(na.omit(col))) <= 2)
  }

  # If any variable is not binary, stop and return an error message
  non_binary_vars <- c(xVars, yVars)[!checkBinaryVars(c(xVars, yVars))]
  if (length(non_binary_vars) > 0) {
    stop("The following variables are not binary (they do not have exactly two unique values): ", paste(non_binary_vars, collapse = ", "))
  }

  # Initialize variables
  removediag <- FALSE

  # If yVars is NULL, set yVars equal to xVars and remove diagonal
  if (is.null(yVars)) {
    yVars <- xVars
    removediag <- TRUE
  }

  # Define Variables
  Variables <- c(xVars, yVars)
  Data <- ReplaceMissingLabels(Data)  # Replace missing labels in the data

  # If xVars is not provided, determine binary variables
  if (is.null(xVars)) {
    xVars <- getBinaryVars(Data, Ordinal = TRUE)
    if (!Ordinal) {
      xVars <- getBinaryVars(Data, Ordinal = FALSE)
    }
  }

  # Subset the data for xVars and yVars
  xdata <- Data[xVars]
  ydata <- Data[yVars]

  # Initialize matrices to store Phi coefficients, p-values, and pair counts
  MC <- matrix(nrow = ncol(xdata), ncol = ncol(ydata))
  rownames(MC) <- colnames(xdata)
  colnames(MC) <- colnames(ydata)
  MP <- MC  # p-values matrix
  MN <- MC  # pair counts matrix

  # Loop through all combinations of xVars and yVars to compute Phi and Chi-square tests
  for (xi in 1:ncol(xdata)) {
    for (yi in 1:ncol(ydata)) {
      phi <- NA
      chi_p <- NA
      x <- xdata[, xi]
      y <- ydata[, yi]
      d <- cbind(x, y)

      d <- na.omit(d)  # Remove missing data
      hadError <- FALSE
      MN[xi, yi] <- nrow(d)

      tryCatch({
        if (nrow(d) > 0) {
          contingency_table <- table(d[, 1], d[, 2])
          phi <- suppressWarnings(cor(as.numeric(d[, 1]), as.numeric(d[, 2])))  # Calculate Phi coefficient
          chi_test <- chisq.test(contingency_table)  # Perform chi-square test
          chi_p <- chi_test$p.value
        }
        else {
          phi = NA
        }
      }, error = function(e) {
        hadError <- TRUE  # If there's an error, mark it
      }, interrupt = function() {
        hadError <- FALSE
      })

      # Store results in the matrices
      if (hadError) {
        MC[xi, yi] <- NA
        MP[xi, yi] <- NA
      }
      else {
        MC[xi, yi] <- phi
        MP[xi, yi] <- chi_p
      }
    }
  }

  # Apply FDR correction to p-values
  pvals <- MP
  if (removediag) {
    diag(pvals) <- NaN  # Remove diagonal entries if specified
    diag(MC) <- NaN
  }

  # Store unadjusted results
  M <- list()
  M$phi <- MC
  M$p <- pvals
  M$npairs <- MN

  # Adjust p-values using FDR method
  M_FDR <- M
  padjusted <- p.adjust(MP, method = "fdr")
  M_FDR$p <- matrix(padjusted, ncol = length(yVars))
  pvals_FDR <- M_FDR$p
  if (removediag) {
    diag(M_FDR$p) <- NaN
    diag(M_FDR$phi) <- NaN
  }

  # Set row and column names for the results
  colnames(M$phi) <- yVars
  rownames(M$phi) <- xVars
  colnames(M$p) <- yVars
  rownames(M$p) <- xVars
  colnames(M_FDR$phi) <- yVars
  rownames(M_FDR$phi) <- xVars
  colnames(M_FDR$p) <- yVars
  rownames(M_FDR$p) <- xVars

  # Reshape the results for plotting
  plot.data_R <- pivot_longer(as.data.frame(M$phi) %>% rownames_to_column(var = "XVar"),
                              cols = colnames(M$phi), names_to = "YVar", values_to = "Phi")
  plot.data_P <- pivot_longer(as.data.frame(M$p) %>% rownames_to_column(var = "XVar"),
                              cols = colnames(M$p), names_to = "YVar", values_to = "P")
  plot.data_P_adj <- pivot_longer(as.data.frame(M_FDR$p) %>%
                                    rownames_to_column(var = "XVar"), cols = colnames(M$p),
                                  names_to = "YVar", values_to = "P_adj")
  plot.data_npairs <- pivot_longer(as.data.frame(M$npairs) %>%
                                     rownames_to_column(var = "XVar"), cols = colnames(M$npairs),
                                   names_to = "YVar", values_to = "nPairs")
  plot.data <- suppressMessages(left_join(plot.data_R, plot.data_P) %>% left_join(plot.data_P_adj) %>%
    left_join(plot.data_npairs))

  # Assign stars based on p-value thresholds
  plot.data$stars <- cut(plot.data$P, breaks = c(-Inf, 0.001,
                                                 0.01, 0.05, Inf), label = c("***", "**", "*", ""))
  plot.data$stars_FDR <- cut(plot.data$P_adj, breaks = c(-Inf,
                                                         0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

  # Relabel variables if requested
  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[plot.data$XVar],
                                     def.value = colnames(Data[plot.data$XVar])) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")
    ylabels <- sjlabelled::get_label(Data[plot.data$YVar],
                                     def.value = colnames(Data[plot.data$YVar])) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")
    plot.data$XLabel <- xlabels$label
    plot.data$YLabel <- ylabels$label
  }
  else {
    plot.data$XLabel <- plot.data$XVar
    plot.data$YLabel <- plot.data$YVar
  }

  # Factorize the labels for plotting
  plot.data$XLabel <- factor(plot.data$XLabel, ordered = FALSE,
                             levels = unique(plot.data$XLabel))
  plot.data$YLabel <- factor(plot.data$YLabel, ordered = FALSE,
                             levels = rev(unique(plot.data$YLabel)))

  PlotText <- paste("YVar", plot.data$YVar, "</br> XVAR: ",
                    plot.data$XVar, "</br> P-Value: ", plot.data$P, plot.data$stars,
                    "</br> FDR-corrected P: ", plot.data$P_adj, plot.data$stars_FDR)
  P <- plot.data %>% ggplot(aes(y = XLabel, x = YLabel, fill = Phi,
                                text = PlotText)) + geom_tile() +
    geom_text(aes(label = stars), color = "black") +
    scale_fill_gradient2(limits = c(-1,  1), low = scales::muted("purple"),  high = scales::muted("green")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),   axis.text.y = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  P_FDR <- plot.data %>% ggplot(aes(y = XLabel, x = YLabel,
                                    fill = Phi, text = PlotText)) + geom_tile() +
    geom_text(aes(label = stars_FDR),
              color = "black") +
    scale_fill_gradient2(limits = c(-1,  1), low = scales::muted("purple"),  high = scales::muted("green")) +

    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  M$plot <- P
  M_FDR$plot <- P_FDR
  BinaryMapping <- createBinaryMapping(Data, unique(c(xVars, yVars)))
  return(list(Unadjusted = M, FDRCorrected = M_FDR, method = "Phi",
              Relabel = Relabel, Covariates = NULL, BinaryMapping = BinaryMapping))
  }
