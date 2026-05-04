#' Plot correlations heatmap
#'
#' Calculates correlations (or partial correlations) and plots a heatmap.
#' Fully label-aware, robust to non-syntactic names, and safe for real-world data.
#'
#' @param Data A data frame.
#' @param xVars Character vector of x variables.
#' @param yVars Character vector of y variables.
#' @param covars Optional covariates.
#' @param method Correlation method.
#' @param Relabel Use labels if available.
#' @param Ordinal Convert ordinal variables.
#' @param min_n Minimum N required.
#' @param eps Variance tolerance.
#' @return A list with correlation matrices and plots.
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

  # ---------------------------
  # Validate inputs
  # ---------------------------
  if (!is.data.frame(Data)) stop("Data must be a data.frame.")

  method <- tolower(method)
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("Invalid method.")
  }

  Data <- ReplaceMissingLabels(Data)

  # ---------------------------
  # Determine variables
  # ---------------------------
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

  covars_in <- covars
  covars <- intersect(as.character(covars %||% character(0)), names(Data))
  covars_missing <- setdiff(covars_in, covars)

  # ---------------------------
  # Early exit
  # ---------------------------
  if (length(xVars) == 0 || length(yVars) == 0) {

    empty_plot <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::theme_bw() +
      ggplot2::labs(subtitle = "No valid variables for correlation")

    z <- matrix(numeric(0), 0, 0)

    return(list(
      Unadjusted = list(r = z, p = z, npairs = z, plot = empty_plot),
      FDRCorrected = list(r = z, p = z, npairs = z, plot = empty_plot),
      method = method,
      Relabel = Relabel,
      Covariates = covars,
      CovariatesMissing = covars_missing
    ))
  }

  # ---------------------------
  # Prepare data (safe numeric coercion)
  # ---------------------------
  prep_numeric <- function(df, vars) {
    out <- df[, vars, drop = FALSE]
    names_out <- names(out)

    out <- as.data.frame(lapply(out, function(z) {
      if (is.numeric(z)) return(z)
      suppressWarnings(as.numeric(as.character(z)))
    }), check.names = FALSE)

    names(out) <- names_out
    out <- as.data.frame(lapply(out, function(z) ifelse(is.finite(z), z, NA_real_)),
                         check.names = FALSE)
    out
  }

  xdata <- prep_numeric(Data, xVars)
  ydata <- prep_numeric(Data, yVars)

  # ---------------------------
  # Initialize outputs
  # ---------------------------
  MC <- matrix(NA_real_, length(xVars), length(yVars),
               dimnames = list(xVars, yVars))
  MP <- MC
  MN <- MC

  # ---------------------------
  # Core computation
  # ---------------------------
  for (xi in seq_along(xVars)) {
    for (yi in seq_along(yVars)) {

      vars_needed <- unique(c(xVars[xi], yVars[yi], covars))
      df_tmp <- stats::na.omit(Data[, vars_needed, drop = FALSE])

      MN[xi, yi] <- nrow(df_tmp)

      est <- NA_real_
      pval <- NA_real_

      if (nrow(df_tmp) >= min_n) {

        x <- df_tmp[[xVars[xi]]]
        y <- df_tmp[[yVars[yi]]]

        if (stats::var(x) > eps && stats::var(y) > eps) {

          if (length(covars) > 0) {

            mf <- stats::model.frame(~ . - 1,
                                     data = df_tmp[, covars, drop = FALSE],
                                     na.action = stats::na.pass)

            Z <- stats::model.matrix(~ . - 1, data = mf)
            cdata <- as.data.frame(Z, check.names = FALSE)

            keep <- vapply(cdata, function(z) {
              v <- stats::var(z)
              is.finite(v) && v > eps
            }, logical(1))

            cdata <- cdata[, keep, drop = FALSE]

            if (ncol(cdata) > 0 && nrow(df_tmp) >= (ncol(cdata) + 3)) {

              df_model <- cbind(x = x, y = y, cdata)
              cov_terms <- paste0("`", colnames(cdata), "`", collapse = " + ")

              r1 <- tryCatch(
                residuals(lm(as.formula(paste("x ~", cov_terms)), data = df_model)),
                error = function(e) NULL
              )

              r2 <- tryCatch(
                residuals(lm(as.formula(paste("y ~", cov_terms)), data = df_model)),
                error = function(e) NULL
              )

              if (!is.null(r1) && !is.null(r2)) {
                tmp <- tryCatch(cor.test(r1, r2, method = method), error = function(e) NULL)
                if (!is.null(tmp)) {
                  est <- unname(tmp$estimate)
                  pval <- tmp$p.value
                }
              }
            }

            if (is.na(est)) {
              tmp <- tryCatch(cor.test(x, y, method = method), error = function(e) NULL)
              if (!is.null(tmp)) {
                est <- unname(tmp$estimate)
                pval <- tmp$p.value
              }
            }

          } else {
            tmp <- tryCatch(cor.test(x, y, method = method), error = function(e) NULL)
            if (!is.null(tmp)) {
              est <- unname(tmp$estimate)
              pval <- tmp$p.value
            }
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
  }

  MC[MN < 2] <- NaN
  MP[MN < 2] <- NaN

  # ---------------------------
  # FDR correction
  # ---------------------------
  pvec <- as.vector(MP)
  padj <- rep(NA_real_, length(pvec))
  ok <- is.finite(pvec)
  padj[ok] <- stats::p.adjust(pvec[ok], method = "fdr")

  M_FDR <- matrix(padj, nrow = nrow(MP), ncol = ncol(MP),
                  dimnames = dimnames(MP))

  # ---------------------------
  # Build plot data
  # ---------------------------
  plot.data <- expand.grid(XVar = xVars, YVar = yVars, stringsAsFactors = FALSE)

  plot.data$R <- as.vector(MC)
  plot.data$P <- as.vector(MP)
  plot.data$P_adj <- as.vector(M_FDR)
  plot.data$nPairs <- as.vector(MN)

  plot.data$stars <- cut(plot.data$P,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***","**","*",""))

  plot.data$stars_FDR <- cut(plot.data$P_adj,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***","**","*",""))

  # ---------------------------
  # Labels
  # ---------------------------
  if (Relabel) {
    xlabels <- sjlabelled::get_label(Data[, xVars, drop = FALSE], def.value = xVars)
    ylabels <- sjlabelled::get_label(Data[, yVars, drop = FALSE], def.value = yVars)

    plot.data$XLabel <- xlabels[plot.data$XVar]
    plot.data$YLabel <- ylabels[plot.data$YVar]

    plot.data$XLabel[is.na(plot.data$XLabel)] <- plot.data$XVar[is.na(plot.data$XLabel)]
    plot.data$YLabel[is.na(plot.data$YLabel)] <- plot.data$YVar[is.na(plot.data$YLabel)]

  } else {
    plot.data$XLabel <- plot.data$XVar
    plot.data$YLabel <- plot.data$YVar
  }

  plot.data$XLabel <- factor(plot.data$XLabel, levels = rev(unique(plot.data$XLabel)))
  plot.data$YLabel <- factor(plot.data$YLabel, levels = unique(plot.data$YLabel))

  # ---------------------------
  # Plot
  # ---------------------------
  P <- ggplot2::ggplot(plot.data,
                       ggplot2::aes(x = YLabel, y = XLabel, fill = R)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = stars)) +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), na.value = "grey90") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title = ggplot2::element_blank())

  P_FDR <- ggplot2::ggplot(plot.data,
                           ggplot2::aes(x = YLabel, y = XLabel, fill = R)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = stars_FDR)) +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), na.value = "grey90") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title = ggplot2::element_blank())

  # ---------------------------
  # Return
  # ---------------------------
  list(
    Unadjusted = list(r = MC, p = MP, npairs = MN, plot = P),
    FDRCorrected = list(r = MC, p = M_FDR, npairs = MN, plot = P_FDR),
    method = method,
    Relabel = Relabel,
    Covariates = covars,
    CovariatesMissing = covars_missing
  )
}
