#' Plot correlations heatmap
#'
#' Calculate correlations (or partial correlations if covariates are supplied) between variables
#' and plot them as a heatmap with significance stars.
#'
#' This version is robust to non-syntactic column names (spaces, hyphens, Greek letters like β),
#' by preserving names during data.frame coercions and using tidyselect::all_of() in pivot_longer().
#'
#' @param Data A data frame containing variables to correlate.
#' @param xVars Character vector of x-axis variable names. If NULL, numeric variables are auto-selected.
#' @param yVars Character vector of y-axis variable names. If NULL, uses xVars and removes diagonal.
#' @param covars Character vector of covariate names. If provided, computes partial correlations via ppcor.
#' @param method Correlation method: "pearson", "spearman", or "kendall".
#' @param Relabel Logical. If TRUE, uses sjlabelled variable labels when present.
#' @param Ordinal Logical. If TRUE, ordinal variables are converted to numeric using ConvertOrdinalToNumeric().
#' @param min_n Minimum complete cases required for a correlation (or partial correlation).
#' @param eps Small tolerance to treat near-constant variables as constant.
#'
#' @return A list with:
#' \itemize{
#'   \item Unadjusted: list(r, p, npairs, plot)
#'   \item FDRCorrected: list(r, p, npairs, plot)
#'   \item method, Relabel, Covariates, CovariatesMissing
#' }
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
  if (!is.data.frame(Data)) {
    stop("Data must be a data.frame.")
  }

  method <- tolower(method)
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("method must be one of: 'pearson', 'spearman', 'kendall'.")
  }

  # Ensure labels exist where expected (package-specific behavior)
  Data <- ReplaceMissingLabels(Data)

  # ---------------------------
  # Determine vars
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

  # Keep only vars that actually exist in Data
  xVars <- intersect(as.character(xVars), names(Data))
  yVars <- intersect(as.character(yVars), names(Data))

  # Covariates: track missing covars (report only, no crash)
  covars_in <- covars
  if (!is.null(covars) && length(covars) > 0) {
    covars <- intersect(as.character(covars), names(Data))
  } else {
    covars <- character(0)
  }
  covars_missing <- if (!is.null(covars_in) && length(covars_in) > 0) {
    setdiff(as.character(covars_in), covars)
  } else {
    character(0)
  }

  # ---------------------------
  # Early exit if nothing to correlate
  # ---------------------------
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
      Covariates = covars,
      CovariatesMissing = covars_missing
    ))
  }

  # ---------------------------
  # Prepare data
  # ---------------------------

  # Optional ordinal conversion (package-specific function)
  if (isTRUE(Ordinal)) {
    Variables <- unique(c(xVars, yVars))
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # Build xdata and ydata and COERCE TO NUMERIC while PRESERVING NAMES.
  # This is critical for weird names like "IFN-β" or "IL-10" not becoming "IFN.β" or "IL.10".
  xdata <- Data[, xVars, drop = FALSE]
  ydata <- Data[, yVars, drop = FALSE]

  x_names <- names(xdata)
  y_names <- names(ydata)

  xdata <- as.data.frame(
    lapply(xdata, function(z) {
      if (is.numeric(z) || is.integer(z)) return(as.numeric(z))
      suppressWarnings(as.numeric(as.character(z)))
    }),
    check.names = FALSE
  )
  names(xdata) <- x_names

  ydata <- as.data.frame(
    lapply(ydata, function(z) {
      if (is.numeric(z) || is.integer(z)) return(as.numeric(z))
      suppressWarnings(as.numeric(as.character(z)))
    }),
    check.names = FALSE
  )
  names(ydata) <- y_names

  # Scrub non-finite values globally (render-safe)
  xdata <- as.data.frame(lapply(xdata, function(z) ifelse(is.finite(z), z, NA_real_)), check.names = FALSE)
  names(xdata) <- x_names

  ydata <- as.data.frame(lapply(ydata, function(z) ifelse(is.finite(z), z, NA_real_)), check.names = FALSE)
  names(ydata) <- y_names

  # ---------------------------
  # Prepare covariates
  # ---------------------------
  if (length(covars) > 0) {
    cdata_raw <- Data[, covars, drop = FALSE]

    # Use model.frame with na.pass to prevent row dropping, then model.matrix for dummy coding.
    mf <- stats::model.frame(~ . - 1, data = cdata_raw, na.action = stats::na.pass)
    Z <- stats::model.matrix(~ . - 1, data = mf)

    # Preserve names (model.matrix usually makes syntactic names, but we only need numeric stability)
    cdata <- as.data.frame(Z, check.names = FALSE)

    # Coerce to numeric and keep NAs
    cdata <- as.data.frame(lapply(cdata, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)
    cdata <- as.data.frame(lapply(cdata, function(z) ifelse(is.finite(z), z, NA_real_)), check.names = FALSE)

    # Drop columns that are all NA
    keep <- vapply(cdata, function(z) any(!is.na(z)), logical(1))
    cdata <- cdata[, keep, drop = FALSE]

    # Drop near-constant columns (prevents NaNs and singularities)
    if (ncol(cdata) > 0) {
      v <- vapply(cdata, function(z) stats::var(z, na.rm = TRUE), numeric(1))
      keep2 <- is.finite(v) & (v > eps)
      cdata <- cdata[, keep2, drop = FALSE]
    }

    # Optional scaling for numeric stability (partial correlation does not require it, but helps conditioning)
    if (ncol(cdata) > 0) {
      sds <- vapply(cdata, stats::sd, numeric(1), na.rm = TRUE)
      sds[!is.finite(sds) | sds <= eps] <- 1
      cdata <- as.data.frame(scale(cdata, center = TRUE, scale = sds), check.names = FALSE)
    }
  } else {
    cdata <- data.frame(check.names = FALSE)
  }

  # Ensure row alignment (should hold with na.pass)
  if (ncol(cdata) > 0 && nrow(cdata) != nrow(Data)) {
    stop("Covariate matrix row count does not match Data row count. This should not happen with na.pass.")
  }

  # ---------------------------
  # Build outputs (correlations + p-values)
  # ---------------------------
  MC <- matrix(
    NA_real_,
    nrow = ncol(xdata), ncol = ncol(ydata),
    dimnames = list(colnames(xdata), colnames(ydata))
  )
  MP <- MC
  MN <- MC

  has_cov <- (ncol(cdata) > 0)

  for (xi in seq_len(ncol(xdata))) {
    for (yi in seq_len(ncol(ydata))) {

      x <- xdata[[xi]]
      y <- ydata[[yi]]

      d <- if (has_cov) cbind(x, y, cdata) else cbind(x, y)
      d <- stats::na.omit(d)

      MN[xi, yi] <- nrow(d)

      est <- NA_real_
      pval <- NA_real_

      # For partial correlation, you need enough rows relative to number of covariates
      min_needed <- if (has_cov) max(min_n, ncol(cdata) + 3) else min_n

      if (nrow(d) >= min_needed) {

        vx <- stats::var(d[, 1])
        vy <- stats::var(d[, 2])

        ok_var <- is.finite(vx) && is.finite(vy) && (vx > eps) && (vy > eps)

        if (ok_var) {
          if (has_cov) {
            tmp <- tryCatch(
              ppcor::pcor.test(
                d[, 1], d[, 2],
                d[, seq(3, ncol(d)), drop = FALSE],
                method = method
              ),
              error = function(e) NULL
            )
            if (!is.null(tmp)) {
              est <- unname(tmp$estimate)
              pval <- tmp$p.value
            }
          } else {
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
      }

      MC[xi, yi] <- est
      MP[xi, yi] <- pval
    }
  }

  # Remove diagonal if xVars == yVars
  if (removediag) {
    diag(MC) <- NaN
    diag(MP) <- NaN
    diag(MN) <- NaN
  }

  # Mark unusable cells (very low n) as NaN to avoid downstream issues
  MC[MN < 2] <- NaN
  MP[MN < 2] <- NaN

  # ---------------------------
  # FDR adjustment
  # ---------------------------
  pvec <- as.vector(MP)
  padj <- rep(NA_real_, length(pvec))
  okp <- is.finite(pvec)
  if (any(okp)) {
    padj[okp] <- stats::p.adjust(pvec[okp], method = "fdr")
  }

  M_FDR_p <- matrix(padj, nrow = nrow(MP), ncol = ncol(MP), dimnames = dimnames(MP))
  if (removediag) diag(M_FDR_p) <- NaN

  # ---------------------------
  # Build plot data (long form) safely for weird names
  # ---------------------------

  # as.data.frame(matrix) will otherwise name-repair when check.names is TRUE
  df_MC <- tibble::rownames_to_column(as.data.frame(MC, check.names = FALSE), var = "XVar")
  df_MP <- tibble::rownames_to_column(as.data.frame(MP, check.names = FALSE), var = "XVar")
  df_MF <- tibble::rownames_to_column(as.data.frame(M_FDR_p, check.names = FALSE), var = "XVar")
  df_MN <- tibble::rownames_to_column(as.data.frame(MN, check.names = FALSE), var = "XVar")

  plot.data_R <- tidyr::pivot_longer(
    df_MC,
    cols = tidyselect::all_of(colnames(MC)),
    names_to = "YVar", values_to = "R"
  )

  plot.data_P <- tidyr::pivot_longer(
    df_MP,
    cols = tidyselect::all_of(colnames(MP)),
    names_to = "YVar", values_to = "P"
  )

  plot.data_P_adj <- tidyr::pivot_longer(
    df_MF,
    cols = tidyselect::all_of(colnames(M_FDR_p)),
    names_to = "YVar", values_to = "P_adj"
  )

  plot.data_npairs <- tidyr::pivot_longer(
    df_MN,
    cols = tidyselect::all_of(colnames(MN)),
    names_to = "YVar", values_to = "nPairs"
  )

  plot.data <- suppressMessages(
    plot.data_R %>%
      dplyr::left_join(plot.data_P, by = c("XVar", "YVar")) %>%
      dplyr::left_join(plot.data_P_adj, by = c("XVar", "YVar")) %>%
      dplyr::left_join(plot.data_npairs, by = c("XVar", "YVar"))
  )

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

  # ---------------------------
  # Apply labels (label-aware, robust to mismatches)
  # ---------------------------
  if (isTRUE(Relabel)) {
    Data <- ReplaceMissingLabels(Data)

    # Intersect protects against any accidental name-repair upstream
    x_unique <- intersect(unique(as.character(plot.data$XVar)), names(Data))
    y_unique <- intersect(unique(as.character(plot.data$YVar)), names(Data))

    xlabels <- sjlabelled::get_label(Data[, x_unique, drop = FALSE], def.value = x_unique)
    ylabels <- sjlabelled::get_label(Data[, y_unique, drop = FALSE], def.value = y_unique)

    plot.data$XLabel <- unname(xlabels[as.character(plot.data$XVar)])
    plot.data$YLabel <- unname(ylabels[as.character(plot.data$YVar)])

    # Any that did not match get_label fallback to variable name
    plot.data$XLabel[is.na(plot.data$XLabel)] <- plot.data$XVar[is.na(plot.data$XLabel)]
    plot.data$YLabel[is.na(plot.data$YLabel)] <- plot.data$YVar[is.na(plot.data$YLabel)]
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
    "<br><b>nPairs:</b> ", ifelse(is.na(plot.data$nPairs) | is.nan(plot.data$nPairs), "NA", plot.data$nPairs),
    if (length(covars) > 0) paste0("<br><b>Covars:</b> ", paste(covars, collapse = ", ")) else "",
    if (length(covars_missing) > 0) paste0("<br><b>Missing covars:</b> ", paste(covars_missing, collapse = ", ")) else ""
  )

  # ---------------------------
  # Plot
  # ---------------------------
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

  # ---------------------------
  # Return result
  # ---------------------------
  M <- list(r = MC, p = MP, npairs = MN, plot = P)
  M_FDR <- list(r = MC, p = M_FDR_p, npairs = MN, plot = P_FDR)

  list(
    Unadjusted = M,
    FDRCorrected = M_FDR,
    method = method,
    Relabel = Relabel,
    Covariates = covars,
    CovariatesMissing = covars_missing
  )
}
