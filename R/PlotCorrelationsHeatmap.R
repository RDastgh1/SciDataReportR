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
    Ordinal = FALSE,
    min_n = 3,
    eps = 1e-12
) {
  Data <- ReplaceMissingLabels(Data)

  # ---- determine vars ----
  if (is.null(xVars)) {
    xVars <- getNumVars(Data, Ordinal = isTRUE(Ordinal))
  }
  if (is.null(yVars)) {
    yVars <- xVars
    removediag <- TRUE
  } else {
    removediag <- FALSE
  }

  xVars <- intersect(as.character(xVars), names(Data))
  yVars <- intersect(as.character(yVars), names(Data))

  covars <- if (!is.null(covars) && length(covars) > 0) intersect(as.character(covars), names(Data)) else character(0)

  # ---- early exit if nothing to correlate ----
  if (length(xVars) == 0 || length(yVars) == 0) {
    empty_plot_data <- data.frame(
      XVar = character(0), YVar = character(0),
      R = numeric(0), P = numeric(0), P_adj = numeric(0),
      nPairs = integer(0),
      stars = character(0), stars_FDR = character(0),
      XLabel = character(0), YLabel = character(0),
      PlotText = character(0),
      stringsAsFactors = FALSE
    )

    P0 <- ggplot2::ggplot(empty_plot_data, ggplot2::aes(x = YLabel, y = XLabel)) +
      ggplot2::geom_blank() +
      ggplot2::theme_bw() +
      ggplot2::labs(subtitle = "No numeric variables available for correlation") +
      ggplot2::xlab("") + ggplot2::ylab("")

    z <- matrix(numeric(0), nrow = 0, ncol = 0)
    M0 <- list(r = z, p = z, npairs = z, plot = P0)

    return(list(
      Unadjusted = M0,
      FDRCorrected = M0,
      method = method,
      Relabel = Relabel,
      Covariates = covars
    ))
  }

  # ---- ordinal handling ----
  if (isTRUE(Ordinal)) {
    Variables <- unique(c(xVars, yVars))
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # ---- prepare covariates data ----
  if (length(covars) > 0) {
    cdata <- Data[, covars, drop = FALSE]

    # characters to factor; factors to numeric codes; then scale
    char_vars <- vapply(cdata, is.character, logical(1))
    if (any(char_vars)) cdata[char_vars] <- lapply(cdata[char_vars], as.factor)

    factor_vars <- vapply(cdata, is.factor, logical(1))
    if (any(factor_vars)) cdata[factor_vars] <- lapply(cdata[factor_vars], as.numeric)

    # numeric coercion safety
    cdata <- as.data.frame(lapply(cdata, function(z) suppressWarnings(as.numeric(as.character(z)))))
    cdata <- as.data.frame(scale(cdata))
  } else {
    cdata <- data.frame()
  }

  xdata <- Data[, xVars, drop = FALSE]
  ydata <- Data[, yVars, drop = FALSE]

  # ---- coerce x/y to numeric safely (drop truly non-numeric later at pair level) ----
  xdata <- as.data.frame(lapply(xdata, function(z) {
    if (is.numeric(z) || is.integer(z)) return(as.numeric(z))
    suppressWarnings(as.numeric(as.character(z)))
  }))
  ydata <- as.data.frame(lapply(ydata, function(z) {
    if (is.numeric(z) || is.integer(z)) return(as.numeric(z))
    suppressWarnings(as.numeric(as.character(z)))
  }))

  # ---- init matrices ----
  MC <- matrix(NA_real_, nrow = ncol(xdata), ncol = ncol(ydata),
               dimnames = list(colnames(xdata), colnames(ydata)))
  MP <- MC
  MN <- MC

  # ---- compute correlations ----
  for (xi in seq_len(ncol(xdata))) {
    for (yi in seq_len(ncol(ydata))) {

      x <- xdata[[xi]]
      y <- ydata[[yi]]

      d <- if (nrow(cdata) > 0) cbind(x, y, cdata) else cbind(x, y)
      d <- stats::na.omit(d)

      MN[xi, yi] <- nrow(d)

      est <- NA_real_
      pval <- NA_real_

      # need enough rows + non-constant variance
      if (nrow(d) >= min_n && stats::var(d[, 1]) > eps && stats::var(d[, 2]) > eps) {
        if (nrow(cdata) > 0) {
          tmp <- tryCatch(
            ppcor::pcor.test(d[, 1], d[, 2], d[, seq(3, ncol(d)), drop = FALSE], method = method),
            error = function(e) NULL
          )
          if (!is.null(tmp)) {
            est <- unname(tmp$estimate)
            pval <- tmp$p.value
          }
        } else {
          tmp <- tryCatch(stats::cor.test(d[, 1], d[, 2], method = method), error = function(e) NULL)
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

  if (removediag) {
    diag(MC) <- NaN
    diag(MP) <- NaN
    diag(MN) <- NaN
  }

  MC[MN < 2] <- NaN
  MP[MN < 2] <- NaN

  padj <- stats::p.adjust(as.vector(MP), method = "fdr")
  M_FDR_p <- matrix(padj, nrow = nrow(MP), ncol = ncol(MP), dimnames = dimnames(MP))
  if (removediag) diag(M_FDR_p) <- NaN

  # ---- long plot data (guard: if matrices are 0-dim, return empty) ----
  if (nrow(MC) == 0 || ncol(MC) == 0) {
    empty_plot_data <- data.frame(
      XVar = character(0), YVar = character(0),
      R = numeric(0), P = numeric(0), P_adj = numeric(0),
      nPairs = integer(0),
      stars = character(0), stars_FDR = character(0),
      XLabel = character(0), YLabel = character(0),
      PlotText = character(0),
      stringsAsFactors = FALSE
    )

    P0 <- ggplot2::ggplot(empty_plot_data, ggplot2::aes(x = YLabel, y = XLabel)) +
      ggplot2::geom_blank() +
      ggplot2::theme_bw() +
      ggplot2::labs(subtitle = "No numeric variables available for correlation") +
      ggplot2::xlab("") + ggplot2::ylab("")

    z <- matrix(numeric(0), nrow = 0, ncol = 0)
    M0 <- list(r = z, p = z, npairs = z, plot = P0)

    return(list(
      Unadjusted = M0,
      FDRCorrected = M0,
      method = method,
      Relabel = Relabel,
      Covariates = covars
    ))
  }

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

  plot.data$stars <- cut(plot.data$P,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***", "**", "*", ""))
  plot.data$stars_FDR <- cut(plot.data$P_adj,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))

  # ---- labels ----
  if (isTRUE(Relabel)) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[as.character(plot.data$XVar)], def.value = plot.data$XVar) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(plot.data$YVar)], def.value = plot.data$YVar) |>
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

  fmt_num <- function(x, digits = 3) ifelse(is.na(x) | is.nan(x), "NA", formatC(x, format = "f", digits = digits))
  fmt_p <- function(x) ifelse(is.na(x) | is.nan(x), "NA", format.pval(x, digits = 3, eps = 1e-3))

  plot.data$PlotText <- paste0(
    "<b>Y:</b> ", plot.data$YVar,
    "<br><b>X:</b> ", plot.data$XVar,
    "<br><b>r:</b> ", fmt_num(plot.data$R, digits = 3),
    "<br><b>p:</b> ", fmt_p(plot.data$P), " ", plot.data$stars,
    "<br><b>FDR p:</b> ", fmt_p(plot.data$P_adj), " ", plot.data$stars_FDR,
    "<br><b>nPairs:</b> ", ifelse(is.na(plot.data$nPairs) | is.nan(plot.data$nPairs), "NA", plot.data$nPairs)
  )

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

  M <- list(r = MC, p = MP, npairs = MN, plot = P)
  M_FDR <- list(r = MC, p = M_FDR_p, npairs = MN, plot = P_FDR)

  list(Unadjusted = M, FDRCorrected = M_FDR, method = method, Relabel = Relabel, Covariates = covars)
}
