#' Plot Chi-Square Tests for Categorical Associations (optionally stratified by covariates)
#'
#' Conducts Chi-square tests between sets of categorical variables and visualizes the results.
#' NOTE: Chi-square tests do not natively "adjust" for covariates. If `covars` are provided,
#' this function can (optionally) run tests *within strata* (each combination of covariate levels),
#' and combine p-values across strata (Fisher's method) for a single summary p-value per pair.
#' If you need true covariate adjustment, use regression-based models (logistic/multinomial).
#'
#' @param Data A data.frame containing the dataset.
#' @param xVars Character vector of x-axis categorical variables.
#' @param yVars Character vector of y-axis categorical variables. If NULL, uses xVars.
#' @param covars Optional character vector of covariate variables used for stratification (not adjustment).
#' @param Relabel Logical; whether to use variable labels (sjlabelled) in the plot.
#' @param Ordinal Logical; included for backward compatibility (currently unused here).
#' @param Stratify Logical; if TRUE and covars provided, run chi-square within covariate strata and
#'   combine p-values with Fisher's method. If FALSE, covars are ignored.
#' @param MinExpected Minimum expected cell count threshold for chi-square validity warning.
#'   (used to annotate; test still runs).
#'
#' @return A list with:
#' \item{p}{ggplot for unadjusted p-values}
#' \item{pvaltable}{wide table of unadjusted p-values}
#' \item{p_FDR}{ggplot for FDR-adjusted p-values}
#' \item{pvaltable_FDR}{wide table of FDR-adjusted p-values}
#' \item{details}{long table with diagnostics (n, warnings, strata info)}
#' @export
PlotChiSqCovar <- function(
    Data,
    xVars,
    yVars = NULL,
    covars = NULL,
    Relabel = TRUE,
    Ordinal = TRUE,
    Stratify = FALSE,
    MinExpected = 5
) {

  # ---- helpers ----
  .stop_missing <- function(vars, Data) {
    miss <- setdiff(vars, names(Data))
    if (length(miss) > 0) {
      stop("PlotChiSqCovar: Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
  }

  .as_factor <- function(x) {
    # keep NA as NA; preserve labelled values if present
    if (is.factor(x)) return(x)
    as.factor(x)
  }

  .fisher_combine <- function(pvals) {
    pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
    if (length(pvals) == 0) return(NA_real_)
    stat <- -2 * sum(log(pvals))
    df <- 2 * length(pvals)
    stats::pchisq(stat, df = df, lower.tail = FALSE)
  }

  .chisq_one <- function(df, x, y) {
    # df has x and y as factors; NA already removed
    if (nrow(df) == 0) {
      return(list(p = NA_real_, n = 0L, warn = "No complete cases", exp_min = NA_real_))
    }

    tab <- table(df[[x]], df[[y]])
    # Must have at least 2 levels in each to test
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(list(p = NA_real_, n = nrow(df), warn = "Insufficient levels", exp_min = NA_real_))
    }

    # suppress warnings but capture expected counts for diagnostics
    tst <- suppressWarnings(stats::chisq.test(tab))
    exp_min <- tryCatch(min(tst$expected), error = function(e) NA_real_)

    warn <- NA_character_
    if (!is.na(exp_min) && exp_min < MinExpected) {
      warn <- paste0("Low expected counts (min=", signif(exp_min, 3), ")")
    }

    list(p = unname(tst$p.value), n = nrow(df), warn = warn, exp_min = exp_min)
  }

  # ---- input defaults / checks ----
  if (is.null(yVars)) yVars <- xVars

  xVars <- unique(as.character(xVars))
  yVars <- unique(as.character(yVars))
  covars <- if (!is.null(covars) && length(covars) > 0) unique(as.character(covars)) else character(0)

  .stop_missing(unique(c(xVars, yVars, covars)), Data)

  # subset + coerce categorical vars ----
  vars_keep <- unique(c(xVars, yVars, covars))
  DataSubset <- Data[, vars_keep, drop = FALSE]

  # Coerce x/y to factor
  for (v in unique(c(xVars, yVars))) {
    DataSubset[[v]] <- .as_factor(DataSubset[[v]])
  }
  # Coerce covars to factor for stratification (if requested)
  if (length(covars) > 0) {
    for (cv in covars) {
      # keep numeric covars numeric unless stratifying
      if (isTRUE(Stratify)) {
        DataSubset[[cv]] <- .as_factor(DataSubset[[cv]])
      }
    }
  }

  # ---- compute tests for all pairs ----
  pairs <- tidyr::expand_grid(
    XVar = xVars,
    YVar = yVars
  )

  # remove diagonal (XVar == YVar) if overlap
  pairs <- pairs[pairs$XVar != pairs$YVar, , drop = FALSE]

  # main loop
  results <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    xv <- pairs$XVar[[i]]
    yv <- pairs$YVar[[i]]

    # base complete cases for x/y (+ covars if stratifying)
    needed <- unique(c(xv, yv, if (isTRUE(Stratify)) covars else character(0)))
    df0 <- DataSubset[, needed, drop = FALSE]
    df0 <- df0[stats::complete.cases(df0), , drop = FALSE]

    if (nrow(df0) == 0) {
      results[[i]] <- data.frame(
        XVar = xv, YVar = yv,
        pval = NA_real_, n = 0L,
        Test = "Chi Squared",
        warn = "No complete cases",
        strata = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    if (isTRUE(Stratify) && length(covars) > 0) {
      # split by covariate strata, compute pvals, combine with Fisher
      df0$strata <- interaction(df0[, covars, drop = FALSE], drop = TRUE, sep = " | ")

      p_by <- tapply(
        X = seq_len(nrow(df0)),
        INDEX = df0$strata,
        FUN = function(idx) {
          out <- .chisq_one(df0[idx, , drop = FALSE], xv, yv)
          c(p = out$p, n = out$n, exp_min = out$exp_min)
        }
      )

      # unpack
      pvals <- vapply(p_by, function(z) z[["p"]], numeric(1))
      nsum  <- sum(vapply(p_by, function(z) z[["n"]], numeric(1)), na.rm = TRUE)

      p_comb <- .fisher_combine(pvals)

      # diagnostic: any low expected?
      expmins <- vapply(p_by, function(z) z[["exp_min"]], numeric(1))
      warn <- NA_character_
      if (any(is.finite(expmins) & expmins < MinExpected, na.rm = TRUE)) {
        warn <- paste0("Low expected counts in >=1 stratum (min=", signif(min(expmins, na.rm = TRUE), 3), ")")
      }

      results[[i]] <- data.frame(
        XVar = xv, YVar = yv,
        pval = p_comb,
        n = as.integer(nsum),
        Test = "Chi Squared (stratified, Fisher combined)",
        warn = warn,
        strata = paste(names(p_by), collapse = "; "),
        stringsAsFactors = FALSE
      )
    } else {
      out <- .chisq_one(df0, xv, yv)
      results[[i]] <- data.frame(
        XVar = xv, YVar = yv,
        pval = out$p,
        n = as.integer(out$n),
        Test = "Chi Squared",
        warn = out$warn,
        strata = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  stat.test <- dplyr::bind_rows(results)

  # adjust + significance ----
  stat.test <- stat.test |>
    dplyr::mutate(
      pval.adj = stats::p.adjust(pval, method = "fdr"),
      logp = -log10(pval),
      logp_FDR = -log10(pval.adj)
    )

  # Significance bins (match your existing style)
  stat.test$pval.signif <- cut(
    stat.test$pval,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )
  stat.test$pval.adj.signif <- cut(
    stat.test$pval.adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )

  # factor order for consistent legend sizing
  sig_levels <- c("ns", "*", "**", "***", "****")
  # keep compatibility with your shape scale (5 levels); map "***" to "***" etc
  # (We don't generate "****" here, but leave it for drop=FALSE)
  stat.test$`p<.05` <- factor(stat.test$pval.signif, levels = c("ns", "*", "**", "***", "****"))
  stat.test$pval.adj.signif <- factor(stat.test$pval.adj.signif, levels = c("ns", "*", "**", "***", "****"))

  # cap logp for point size limits like your original
  stat.test$logp_lim <- dplyr::if_else(is.finite(stat.test$logp) & stat.test$logp > 5, 5, stat.test$logp)
  stat.test$logp_FDR_lim <- dplyr::if_else(is.finite(stat.test$logp_FDR) & stat.test$logp_FDR > 5, 5, stat.test$logp_FDR)

  # labels ----
  if (isTRUE(Relabel)) {
    Data <- ReplaceMissingLabels(Data)

    xlabels <- sjlabelled::get_label(Data[as.character(stat.test$XVar)], def.value = stat.test$XVar) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(stat.test$YVar)], def.value = stat.test$YVar) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$XVar
    stat.test$YLabel <- stat.test$YVar
  }

  # Plot hover text
  PlotText <- paste0(
    "<br><b>XVar Label:</b> ", stat.test$XLabel,
    "<br><b>XVar:</b> ", stat.test$XVar,
    "<br><b>YVar Label:</b> ", stat.test$YLabel,
    "<br><b>YVar:</b> ", stat.test$YVar,
    "<br><b>p:</b> ", stat.test$pval, " ", stat.test$pval.signif,
    "<br><b>FDR p:</b> ", stat.test$pval.adj, " ", stat.test$pval.adj.signif,
    "<br><b>n:</b> ", stat.test$n,
    ifelse(!is.na(stat.test$warn), paste0("<br><b>Note:</b> ", stat.test$warn), ""),
    ifelse(isTRUE(Stratify) && !is.na(stat.test$strata), paste0("<br><b>Strata:</b> ", stat.test$strata), "")
  )

  # For axis ordering, keep original xVars/yVars order; if Relabel, use labels
  if (isTRUE(Relabel)) {
    x_levels <- sjlabelled::get_label(Data[as.character(xVars)], def.value = xVars)
    y_levels <- sjlabelled::get_label(Data[as.character(yVars)], def.value = yVars)
  } else {
    x_levels <- xVars
    y_levels <- yVars
  }

  stat.test$XLabel <- factor(stat.test$XLabel, levels = y_levels) # x-axis variables (yVars) in your original plot
  stat.test$YLabel <- factor(stat.test$YLabel, levels = x_levels) # y-axis variables (xVars)


  p_g <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = `p<.05`, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = logp_lim, colour = pval)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(
      values = c(7, 16, 17, 15, 18),
      breaks = c("ns", "*", "**", "***", "****"),
      drop = FALSE
    ) +
    ggplot2::scale_color_gradientn(
      trans = "log",
      colours = rev(RColorBrewer::brewer.pal(9, "Oranges")[-1]),
      limits = c(1e-7, 1),
      oob = scales::squish,
      breaks = c(1, 0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
    ) +
    ggplot2::scale_size_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction")

  p_g_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = pval.adj.signif, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = logp_FDR_lim, colour = pval.adj)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(
      values = c(7, 16, 17, 15, 18),
      breaks = c("ns", "*", "**", "***", "****"),
      drop = FALSE
    ) +
    ggplot2::scale_color_gradientn(
      trans = "log",
      colours = rev(RColorBrewer::brewer.pal(9, "Oranges")),
      limits = c(1e-7, 1),
      oob = scales::squish,
      breaks = c(1, 0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
    ) +
    ggplot2::scale_size_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction")


  pvaltable <- stat.test |>
    dplyr::select(YVar, XVar, pval) |>
    tidyr::pivot_wider(names_from = YVar, values_from = pval) |>
    dplyr::filter(stats::complete.cases(XVar))

  pvaltable_FDR <- stat.test |>
    dplyr::select(YVar, XVar, pval.adj) |>
    tidyr::pivot_wider(names_from = YVar, values_from = pval.adj) |>
    dplyr::filter(stats::complete.cases(XVar))

  list(
    p = p_g,
    pvaltable = pvaltable,
    p_FDR = p_g_FDR,
    pvaltable_FDR = pvaltable_FDR,
    details = stat.test
  )
}
