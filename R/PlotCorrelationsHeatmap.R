
#' Plot correlations heatmap
#'
#' Computes correlations or partial correlations and plots a heatmap.
#' Handles:
#' - continuous + categorical covariates
#' - labelled data
#' - non-syntactic names
#' - sparse real-world datasets
#' - ordinal variables
#' - partial correlations via residualization
#'
#' @param Data data.frame
#' @param xVars character vector
#' @param yVars character vector
#' @param covars optional covariates
#' @param method pearson/spearman/kendall
#' @param Relabel use labels
#' @param Ordinal include ordinal vars
#' @param min_n minimum complete rows
#' @param eps variance tolerance
#'
#' @return list
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

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  # =========================================================
  # Validation
  # =========================================================

  if (!is.data.frame(Data)) {
    stop("Data must be a data.frame.")
  }

  method <- tolower(method)

  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("method must be pearson, spearman, or kendall")
  }

  Data <- ReplaceMissingLabels(Data)

  # =========================================================
  # Determine variables
  # =========================================================

  if (is.null(xVars)) {

    xVars <- getNumVars(
      Data,
      Ordinal = isTRUE(Ordinal)
    )
  }

  if (is.null(yVars)) {

    yVars <- xVars
    removediag <- TRUE

  } else {

    removediag <- FALSE
  }

  xVars <- intersect(
    as.character(xVars),
    names(Data)
  )

  yVars <- intersect(
    as.character(yVars),
    names(Data)
  )

  covars_in <- covars

  covars <- intersect(
    as.character(covars %||% character(0)),
    names(Data)
  )

  covars_missing <- setdiff(
    covars_in %||% character(0),
    covars
  )

  all_vars <- unique(
    c(xVars, yVars, covars)
  )

  # =========================================================
  # Early exit
  # =========================================================

  if (
    length(xVars) == 0 ||
    length(yVars) == 0
  ) {

    empty_plot <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        subtitle = "No valid variables"
      )

    z <- matrix(
      numeric(0),
      0,
      0
    )

    return(list(
      Unadjusted = list(
        r = z,
        p = z,
        npairs = z,
        plot = empty_plot
      ),
      FDRCorrected = list(
        r = z,
        p = z,
        npairs = z,
        plot = empty_plot
      ),
      method = method,
      Relabel = Relabel,
      Covariates = covars,
      CovariatesMissing = covars_missing
    ))
  }

  # =========================================================
  # Prepare numeric vars ONLY
  # =========================================================

  corr_vars <- unique(
    c(xVars, yVars)
  )

  DataCorr <- PrepNumericData(
    Data,
    corr_vars
  )

  DataCorr <- as.data.frame(
    DataCorr,
    check.names = FALSE
  )

  # =========================================================
  # Preserve covariates separately
  # =========================================================

  CovarData <- Data[, covars, drop = FALSE]

  # safer categorical handling
  for (nm in names(CovarData)) {

    x <- CovarData[[nm]]

    if (
      is.character(x) ||
      is.logical(x)
    ) {
      CovarData[[nm]] <- factor(x)
    }

    # preserve labelled categorical variables
    if (
      inherits(x, "haven_labelled") &&
      !is.numeric(x)
    ) {
      CovarData[[nm]] <- factor(as.character(x))
    }
  }

  # =========================================================
  # Initialize matrices
  # =========================================================

  MC <- matrix(
    NA_real_,
    nrow = length(xVars),
    ncol = length(yVars),
    dimnames = list(xVars, yVars)
  )

  MP <- MC
  MN <- MC

  # =========================================================
  # Main loop
  # =========================================================

  for (xi in seq_along(xVars)) {

    for (yi in seq_along(yVars)) {

      xname <- xVars[xi]
      yname <- yVars[yi]

      vars_needed <- unique(
        c(xname, yname, covars)
      )

      # -----------------------------------
      # Build temporary dataset
      # -----------------------------------

      df_tmp <- data.frame(
        x = DataCorr[[xname]],
        y = DataCorr[[yname]],
        CovarData[, covars, drop = FALSE],
        check.names = FALSE
      )

      df_tmp <- stats::na.omit(df_tmp)

      n_complete <- nrow(df_tmp)

      MN[xi, yi] <- n_complete

      est <- NA_real_
      pval <- NA_real_

      if (n_complete < min_n) {

        MC[xi, yi] <- NA_real_
        MP[xi, yi] <- NA_real_

        next
      }

      # -----------------------------------
      # Variance checks
      # -----------------------------------

      vx <- suppressWarnings(
        stats::var(df_tmp$x)
      )

      vy <- suppressWarnings(
        stats::var(df_tmp$y)
      )

      if (
        !is.finite(vx) ||
        !is.finite(vy) ||
        is.na(vx) ||
        is.na(vy) ||
        vx <= eps ||
        vy <= eps
      ) {

        MC[xi, yi] <- NA_real_
        MP[xi, yi] <- NA_real_

        next
      }

      # =====================================================
      # PARTIAL CORRELATION
      # =====================================================

      if (length(covars) > 0) {

        # -----------------------------------
        # Build model matrix safely
        # -----------------------------------

        cov_df <- df_tmp[, covars, drop = FALSE]

        mm <- tryCatch(

          stats::model.matrix(
            ~ .,
            data = cov_df
          ),

          error = function(e) NULL
        )

        # fallback if model matrix fails
        if (is.null(mm)) {

          tmp <- tryCatch(

            suppressWarnings(
              stats::cor.test(
                df_tmp$x,
                df_tmp$y,
                method = method
              )
            ),

            error = function(e) NULL
          )

          if (!is.null(tmp)) {

            est <- unname(tmp$estimate)
            pval <- tmp$p.value
          }

        } else {

          # remove intercept
          mm <- mm[, -1, drop = FALSE]

          # remove zero variance columns
          keep <- apply(
            mm,
            2,
            function(z) {

              v <- suppressWarnings(
                stats::var(z)
              )

              is.finite(v) &&
                !is.na(v) &&
                v > eps
            }
          )

          if (length(keep) > 0) {

            mm <- mm[, keep, drop = FALSE]
          }

          # -----------------------------------
          # If ALL covariates collapsed
          # -----------------------------------

          if (ncol(mm) == 0) {

            tmp <- tryCatch(

              suppressWarnings(
                stats::cor.test(
                  df_tmp$x,
                  df_tmp$y,
                  method = method
                )
              ),

              error = function(e) NULL
            )

            if (!is.null(tmp)) {

              est <- unname(tmp$estimate)
              pval <- tmp$p.value
            }

          } else {

            # -----------------------------------
            # Prevent overparameterization
            # -----------------------------------

            if (
              n_complete <= (ncol(mm) + 2)
            ) {

              # fallback to standard correlation
              tmp <- tryCatch(

                suppressWarnings(
                  stats::cor.test(
                    df_tmp$x,
                    df_tmp$y,
                    method = method
                  )
                ),

                error = function(e) NULL
              )

              if (!is.null(tmp)) {

                est <- unname(tmp$estimate)
                pval <- tmp$p.value
              }

            } else {

              # -----------------------------------
              # Residualization
              # -----------------------------------

              fitx <- tryCatch(

                stats::lm.fit(
                  x = cbind(1, mm),
                  y = df_tmp$x
                ),

                error = function(e) NULL
              )

              fity <- tryCatch(

                stats::lm.fit(
                  x = cbind(1, mm),
                  y = df_tmp$y
                ),

                error = function(e) NULL
              )

              if (
                !is.null(fitx) &&
                !is.null(fity)
              ) {

                rx <- fitx$residuals
                ry <- fity$residuals

                tmp <- tryCatch(

                  suppressWarnings(
                    stats::cor.test(
                      rx,
                      ry,
                      method = method
                    )
                  ),

                  error = function(e) NULL
                )

                if (!is.null(tmp)) {

                  est <- unname(tmp$estimate)
                  pval <- tmp$p.value
                }
              }
            }
          }
        }

      } else {

        # =====================================================
        # STANDARD CORRELATION
        # =====================================================

        tmp <- tryCatch(

          suppressWarnings(
            stats::cor.test(
              df_tmp$x,
              df_tmp$y,
              method = method
            )
          ),

          error = function(e) NULL
        )

        if (!is.null(tmp)) {

          est <- unname(tmp$estimate)
          pval <- tmp$p.value
        }
      }

      MC[xi, yi] <- est
      MP[xi, yi] <- pval
    }
  }

  # =========================================================
  # Symmetric cleanup
  # =========================================================

  if (removediag) {

    diag(MC) <- NaN
    diag(MP) <- NaN
  }

  MC[MN < min_n] <- NaN
  MP[MN < min_n] <- NaN

  # =========================================================
  # FDR correction
  # =========================================================

  pvec <- as.vector(MP)

  padj <- rep(
    NA_real_,
    length(pvec)
  )

  ok <- is.finite(pvec)

  if (any(ok)) {

    padj[ok] <- stats::p.adjust(
      pvec[ok],
      method = "fdr"
    )
  }

  M_FDR <- matrix(
    padj,
    nrow = nrow(MP),
    ncol = ncol(MP),
    dimnames = dimnames(MP)
  )

  # =========================================================
  # Plot dataframe
  # =========================================================

  plot.data <- expand.grid(
    XVar = xVars,
    YVar = yVars,
    stringsAsFactors = FALSE
  )

  plot.data$R <- as.vector(MC)
  plot.data$P <- as.vector(MP)
  plot.data$P_adj <- as.vector(M_FDR)
  plot.data$nPairs <- as.vector(MN)

  plot.data$stars <- as.character(

    cut(
      plot.data$P,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    )
  )

  plot.data$stars_FDR <- as.character(

    cut(
      plot.data$P_adj,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    )
  )

  plot.data$stars[
    is.na(plot.data$stars)
  ] <- ""

  plot.data$stars_FDR[
    is.na(plot.data$stars_FDR)
  ] <- ""

  # =========================================================
  # Labels
  # =========================================================

  if (Relabel) {

    xlabels <- stats::setNames(

      sjlabelled::get_label(
        Data[, xVars, drop = FALSE],
        def.value = xVars
      ),

      xVars
    )

    ylabels <- stats::setNames(

      sjlabelled::get_label(
        Data[, yVars, drop = FALSE],
        def.value = yVars
      ),

      yVars
    )

    plot.data$XLabel <- xlabels[
      plot.data$XVar
    ]

    plot.data$YLabel <- ylabels[
      plot.data$YVar
    ]

  } else {

    plot.data$XLabel <- plot.data$XVar
    plot.data$YLabel <- plot.data$YVar
  }

  plot.data$XLabel <- factor(
    plot.data$XLabel,
    levels = rev(unique(plot.data$XLabel))
  )

  plot.data$YLabel <- factor(
    plot.data$YLabel,
    levels = unique(plot.data$YLabel)
  )

  # =========================================================
  # Plot helper
  # =========================================================

  BuildPlot <- function(
    plot.data,
    starvar
  ) {

    ggplot2::ggplot(
      plot.data,
      ggplot2::aes(
        x = YLabel,
        y = XLabel,
        fill = R
      )
    ) +

      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.2
      ) +

      ggplot2::geom_text(
        ggplot2::aes(
          label = .data[[starvar]]
        ),
        size = 4
      ) +

      ggplot2::scale_fill_gradient2(
        high = "#2166AC",
        mid = "white",
        low = "#B2182B",
        midpoint = 0,
        limits = c(-1, 1),
        na.value = "grey90",
        name = "r"
      ) +

      ggplot2::theme_bw() +

      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          hjust = 1
        ),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  P <- BuildPlot(
    plot.data,
    "stars"
  )

  P_FDR <- BuildPlot(
    plot.data,
    "stars_FDR"
  )

  # =========================================================
  # Return
  # =========================================================

  list(

    Unadjusted = list(
      r = MC,
      p = MP,
      npairs = MN,
      plot = P
    ),

    FDRCorrected = list(
      r = MC,
      p = M_FDR,
      npairs = MN,
      plot = P_FDR
    ),

    method = method,
    Relabel = Relabel,
    Covariates = covars,
    CovariatesMissing = covars_missing
  )
}

