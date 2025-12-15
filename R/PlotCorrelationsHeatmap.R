#' Plot Correlations Heatmap
#'
#' This function calculates correlations between variables and plots them as a heatmap.
#'
#' @param Data The dataset containing the variables.
#' @param xVars A character vector of the names of the x-axis variables.
#' @param yVars A character vector of the names of the y-axis variables. Defaults to NULL.
#' @param covars A character vector of the names of covariate variables. Defaults to NULL.
#' @param method The correlation method. Defaults to "pearson".
#' @param Relabel A logical indicating whether to relabel variables. Defaults to TRUE.
#' @param Ordinal Logical, indicating whether ordinal variables should be included.
#' @return A list containing matrices, ggplot objects for visualizations, and details of the method used.
#' @export
PlotCorrelationsHeatmap <- function(
    Data,
    xVars = NULL,
    yVars = NULL,
    covars = NULL,
    method = "pearson",
    Relabel = TRUE,
    Ordinal = FALSE
) {

  removediag <- FALSE
  Data <- ReplaceMissingLabels(Data)

  # Determine variables
  if (is.null(xVars)) {
    xVars <- getNumVars(Data, Ordinal = isTRUE(Ordinal))
  }

  if (is.null(yVars)) {
    yVars <- xVars
    removediag <- TRUE
  }

  # If ordinal requested, convert ordinal vars to numeric (both x and y)
  if (isTRUE(Ordinal)) {
    Variables <- unique(c(xVars, yVars))
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # Prepare covariates data
  if (!is.null(covars) && length(covars) > 0) {
    # make character to factor, factor to numeric for covars (and only covars)
    cdata <- Data[, covars, drop = FALSE]
    char_vars <- vapply(cdata, is.character, logical(1))
    if (any(char_vars)) cdata[char_vars] <- lapply(cdata[char_vars], as.factor)

    factor_vars <- vapply(cdata, is.factor, logical(1))
    if (any(factor_vars)) cdata[factor_vars] <- lapply(cdata[factor_vars], as.numeric)

    # scale covars
    cdata <- scale(cdata)
  } else {
    cdata <- data.frame()
  }

  # X and Y data
  xdata <- Data[, xVars, drop = FALSE]
  ydata <- Data[, yVars, drop = FALSE]

  # Initialize matrices
  MC <- matrix(NA_real_, nrow = ncol(xdata), ncol = ncol(ydata),
               dimnames = list(colnames(xdata), colnames(ydata)))
  MP <- MC
  MN <- MC

  # Compute correlations
  for (xi in seq_len(ncol(xdata))) {
    for (yi in seq_len(ncol(ydata))) {

      x <- xdata[[xi]]
      y <- ydata[[yi]]

      d <- if (nrow(cdata) > 0) cbind(x, y, cdata) else cbind(x, y)
      d <- stats::na.omit(d)

      MN[xi, yi] <- nrow(d)

      # Default output
      est <- NA_real_
      pval <- NA_real_

      # Only compute if enough rows
      if (nrow(d) >= 3) {
        if (nrow(cdata) > 0) {
          # Partial correlation
          tmp <- tryCatch(
            ppcor::pcor.test(
              d[, 1], d[, 2],
              d[, seq(3, 2 + ncol(cdata)), drop = FALSE],
              method = method
            ),
            error = function(e) NULL
          )

          if (!is.null(tmp)) {
            est <- unname(tmp$estimate)
            pval <- tmp$p.value
          }
        } else {
          # Standard correlation
          tmp <- tryCatch(
            stats::cor.test(d[, 1], d[, 2], method = method),
            error = function(e) NULL
          )

          if (!is.null(tmp)) {
            est <- unname(tmp$estimate)
            pval <- tmp$p.value
          }
        }
      }

      MC[xi, yi] <- est
      MP[xi, yi] <- pval
    }
  }

  # Remove diagonal if symmetric heatmap
  if (removediag) {
    diag(MC) <- NaN
    diag(MP) <- NaN
    diag(MN) <- NaN
  }

  # If npairs < 2, blank values
  MC[MN < 2] <- NaN
  MP[MN < 2] <- NaN

  # FDR correction (vectorize then reshape to yVars length)
  padj <- stats::p.adjust(as.vector(MP), method = "fdr")
  M_FDR_p <- matrix(padj, nrow = nrow(MP), ncol = ncol(MP),
                    dimnames = dimnames(MP))

  if (removediag) {
    diag(M_FDR_p) <- NaN
  }

  # Build long plot data (include npairs)
  plot.data_R <- tidyr::pivot_longer(
    tibble::rownames_to_column(as.data.frame(MC), var = "XVar"),
    cols = colnames(MC), names_to = "YVar", values_to = "R"
  )

  plot.data_P <- tidyr::pivot_longer(
    tibble::rownames_to_column(as.data.frame(MP), var = "XVar"),
    cols = colnames(MP), names_to = "YVar", values_to = "P"
  )

  plot.data_P_adj <- tidyr::pivot_longer(
    tibble::rownames_to_column(as.data.frame(M_FDR_p), var = "XVar"),
    cols = colnames(M_FDR_p), names_to = "YVar", values_to = "P_adj"
  )

  plot.data_npairs <- tidyr::pivot_longer(
    tibble::rownames_to_column(as.data.frame(MN), var = "XVar"),
    cols = colnames(MN), names_to = "YVar", values_to = "nPairs"
  )

  plot.data <- suppressMessages(
    dplyr::left_join(plot.data_R, plot.data_P, by = c("XVar", "YVar")) |>
      dplyr::left_join(plot.data_P_adj, by = c("XVar", "YVar")) |>
      dplyr::left_join(plot.data_npairs, by = c("XVar", "YVar"))
  )

  # Significance stars
  plot.data$stars <- cut(
    plot.data$P,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
  )

  plot.data$stars_FDR <- cut(
    plot.data$P_adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
  )

  # Labels (relabel if desired)
  if (isTRUE(Relabel)) {
    Data <- ReplaceMissingLabels(Data)

    xlabels <- sjlabelled::get_label(
      Data[as.character(plot.data$XVar)],
      def.value = plot.data$XVar
    ) |>
      as.data.frame() |>
      tibble::rownames_to_column()

    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(
      Data[as.character(plot.data$YVar)],
      def.value = plot.data$YVar
    ) |>
      as.data.frame() |>
      tibble::rownames_to_column()

    colnames(ylabels) <- c("Variable", "label")

    plot.data$XLabel <- xlabels$label
    plot.data$YLabel <- ylabels$label
  } else {
    plot.data$XLabel <- plot.data$XVar
    plot.data$YLabel <- plot.data$YVar
  }

  plot.data$XLabel <- factor(plot.data$XLabel, ordered = FALSE, levels = rev(unique(plot.data$XLabel)))
  plot.data$YLabel <- factor(plot.data$YLabel, ordered = FALSE, levels = unique(plot.data$YLabel))

  # Plotly-friendly hover formatting helpers
  fmt_num <- function(x, digits = 3) {
    ifelse(is.na(x) | is.nan(x), "NA", formatC(x, format = "f", digits = digits))
  }
  fmt_p <- function(x) {
    ifelse(is.na(x) | is.nan(x), "NA", format.pval(x, digits = 3, eps = 1e-3))
  }

  # Hover text includes r, p, FDR p, and nPairs
  plot.data$PlotText <- paste0(
    "<b>Y:</b> ", plot.data$YVar,
    "<br><b>X:</b> ", plot.data$XVar,
    "<br><b>r:</b> ", fmt_num(plot.data$R, digits = 3),
    "<br><b>p:</b> ", fmt_p(plot.data$P), " ", plot.data$stars,
    "<br><b>FDR p:</b> ", fmt_p(plot.data$P_adj), " ", plot.data$stars_FDR,
    "<br><b>nPairs:</b> ", ifelse(is.na(plot.data$nPairs) | is.nan(plot.data$nPairs), "NA", plot.data$nPairs)
  )

  # Plots
  P <- ggplot2::ggplot(plot.data, ggplot2::aes(y = XLabel, x = YLabel, fill = R, text = PlotText)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = stars), color = "black") +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  P_FDR <- ggplot2::ggplot(plot.data, ggplot2::aes(y = XLabel, x = YLabel, fill = R, text = PlotText)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = stars_FDR), color = "black") +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Return objects in the same structure you already use
  M <- list(r = MC, p = MP, npairs = MN, plot = P)
  M_FDR <- list(r = MC, p = M_FDR_p, npairs = MN, plot = P_FDR)

  list(
    Unadjusted = M,
    FDRCorrected = M_FDR,
    method = method,
    Relabel = Relabel,
    Covariates = covars
  )
}
