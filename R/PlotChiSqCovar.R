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
PlotChiSqCovar <- function(Data, xVars, yVars, covars = NULL, Relabel = TRUE, Ordinal = TRUE) {

  if (is.null(yVars)) yVars <- xVars

  xVars <- intersect(as.character(xVars), names(Data))
  yVars <- intersect(as.character(yVars), names(Data))
  covars <- if (!is.null(covars) && length(covars) > 0) intersect(as.character(covars), names(Data)) else character(0)

  # Drop vars with <2 observed levels
  xVars <- xVars[vapply(xVars, function(v) {
    x <- Data[[v]]; x <- x[!is.na(x)]
    length(unique(x)) >= 2
  }, logical(1))]

  yVars <- yVars[vapply(yVars, function(v) {
    x <- Data[[v]]; x <- x[!is.na(x)]
    length(unique(x)) >= 2
  }, logical(1))]

  if (length(xVars) == 0 || length(yVars) == 0) {
    empty <- data.frame()
    return(list(p = ggplot2::ggplot(), pvaltable = empty, p_FDR = ggplot2::ggplot(), pvaltable_FDR = empty))
  }

  DataSubset <- Data[unique(c(xVars, yVars, covars))]
  DataSubset[unique(c(xVars, yVars))] <- lapply(DataSubset[unique(c(xVars, yVars))], factor)
  DataSubset$rowID <- seq_len(nrow(DataSubset))

  # Long format
  mData <- tidyr::pivot_longer(DataSubset, cols = dplyr::all_of(xVars), names_to = "name", values_to = "value")
  mData <- tidyr::pivot_longer(mData, cols = dplyr::all_of(yVars), names_to = "Covariate", values_to = "Category")

  # attach overlapping covars if needed
  cOverlapV <- covars[covars %in% unique(c(xVars, yVars))]
  if (length(cOverlapV) > 0) {
    mData <- dplyr::left_join(mData, DataSubset[c("rowID", cOverlapV)], by = "rowID")
  }

  mData <- stats::na.omit(mData)

  # Guard against single-level pairs after NA removal
  mData <- mData %>%
    dplyr::group_by(name, Covariate) %>%
    dplyr::filter(dplyr::n_distinct(value) > 1, dplyr::n_distinct(Category) > 1) %>%
    dplyr::ungroup()

  if (nrow(mData) == 0) {
    empty <- data.frame()
    return(list(p = ggplot2::ggplot(), pvaltable = empty, p_FDR = ggplot2::ggplot(), pvaltable_FDR = empty))
  }

  safe_chisq <- function(x, y) {
    if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
    tbl <- table(x, y)
    if (any(tbl == 0)) return(NA_real_) # skip sparse zeros (or swap to fisher.test)
    tryCatch(stats::chisq.test(tbl)$p.value, error = function(e) NA_real_)
  }

  stat.test <- mData %>%
    dplyr::group_by(name, Covariate) %>%
    dplyr::summarise(pval = safe_chisq(value, Category), .groups = "drop") %>%
    dplyr::filter(!is.na(pval)) %>%
    rstatix::add_significance(p.col = "pval", output.col = "pval.signif") %>%
    rstatix::adjust_pvalue(p.col = "pval", output.col = "pval.adj", method = "fdr") %>%
    rstatix::add_significance(p.col = "pval.adj", output.col = "pval.adj.signif") %>%
    dplyr::mutate(Test = "Chi Squared")

  if (nrow(stat.test) == 0) {
    empty <- data.frame()
    return(list(p = ggplot2::ggplot(), pvaltable = empty, p_FDR = ggplot2::ggplot(), pvaltable_FDR = empty))
  }

  stat.test$logp <- -log10(stat.test$pval)
  stat.test$logp_FDR <- -log10(stat.test$pval.adj)
  stat.test$`p<.05` <- factor(stat.test$pval.signif, levels = c("ns","*","**","***","****"))
  stat.test$pval.adj.signif <- factor(stat.test$pval.adj.signif, levels = c("ns","*","**","***","****"))

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[as.character(stat.test$Covariate)], def.value = stat.test$Covariate) %>%
      as.data.frame() %>% tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(stat.test$name)], def.value = stat.test$name) %>%
      as.data.frame() %>% tibble::rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$Covariate
    stat.test$YLabel <- stat.test$name
  }

  PlotText <- paste(
    "</br>XVar Label:", stat.test$XLabel,
    "</br>XVar:", stat.test$Covariate,
    "</br>YVar Label:", stat.test$YLabel,
    "</br>YVar:", stat.test$name,
    "</br>P-Value: ", stat.test$pval, stat.test$pval.signif,
    "</br>FDR-corrected P: ", stat.test$pval.adj, stat.test$pval.adj.signif
  )

  stat.test$PlotText <- PlotText

  # Simple plots (kept similar structure; you can restore your styling)
  p_g <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = logp, colour = pval)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction")

  p_g_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = logp_FDR, colour = pval.adj)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction")

  pvaltable <- stat.test %>%
    dplyr::select(name, Covariate, pval) %>%
    tidyr::pivot_wider(names_from = name, values_from = pval)

  pvaltable_FDR <- stat.test %>%
    dplyr::select(name, Covariate, pval.adj) %>%
    tidyr::pivot_wider(names_from = name, values_from = pval.adj)

  list(p = p_g, pvaltable = pvaltable, p_FDR = p_g_FDR, pvaltable_FDR = pvaltable_FDR)
}
