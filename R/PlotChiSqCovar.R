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
#' @param min_n Minimum number of complete observations required for a tested association.
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
    yVars,
    covars = NULL,
    Relabel = TRUE,
    Ordinal = TRUE,
    min_n = 4
) {

  if (is.null(yVars)) yVars <- xVars

  xVars <- unique(intersect(as.character(xVars), names(Data)))
  yVars <- unique(intersect(as.character(yVars), names(Data)))
  covars <- if (!is.null(covars) && length(covars) > 0) unique(intersect(as.character(covars), names(Data))) else character(0)

  empty_return <- function() {
    empty <- data.frame()
    p0 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_blank() + ggplot2::theme_void()
    list(p = p0, pvaltable = empty, p_FDR = p0, pvaltable_FDR = empty)
  }

  if (length(xVars) == 0 || length(yVars) == 0) return(empty_return())

  keep_cat <- function(v) {
    z <- Data[[v]]
    z <- z[!is.na(z)]
    length(unique(z)) >= 2
  }

  xVars <- xVars[vapply(xVars, keep_cat, logical(1))]
  yVars <- yVars[vapply(yVars, keep_cat, logical(1))]

  if (length(xVars) == 0 || length(yVars) == 0) return(empty_return())

  df0 <- Data %>%
    dplyr::select(dplyr::all_of(unique(c(xVars, yVars, covars)))) %>%
    dplyr::mutate(.rowID = dplyr::row_number()) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(unique(c(xVars, yVars))), ~ as.factor(.x)))

  cov_df <- if (length(covars) > 0) {
    df0 %>% dplyr::select(.rowID, dplyr::all_of(covars))
  } else {
    df0 %>% dplyr::select(.rowID)
  }

  mX <- df0 %>%
    tidyr::pivot_longer(cols = dplyr::all_of(xVars), names_to = "XVar", values_to = "XVal") %>%
    dplyr::select(.rowID, XVar, XVal)

  mY <- df0 %>%
    tidyr::pivot_longer(cols = dplyr::all_of(yVars), names_to = "YVar", values_to = "YVal") %>%
    dplyr::select(.rowID, YVar, YVal)

  mData <- mX %>%
    dplyr::inner_join(mY, by = ".rowID") %>%
    dplyr::left_join(cov_df, by = ".rowID") %>%
    tidyr::drop_na(XVal, YVal)

  if (length(covars) > 0) {
    mData <- mData %>% tidyr::drop_na(dplyr::all_of(covars))
  }

  mData <- mData %>%
    dplyr::group_by(XVar, YVar) %>%
    dplyr::filter(
      dplyr::n() >= min_n,
      dplyr::n_distinct(XVal) > 1,
      dplyr::n_distinct(YVal) > 1
    ) %>%
    dplyr::ungroup()

  if (nrow(mData) == 0) return(empty_return())

  stat.test <- mData %>%
    dplyr::group_by(XVar, YVar) %>%
    dplyr::summarise(
      test = list(tryCatch(stats::chisq.test(XVal, YVal), error = function(e) NULL)),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pval = if (!is.null(test)) test$p.value else NA_real_,
      cramers_v = if (!is.null(test)) {
        tbl <- table(Data[[XVar]], Data[[YVar]])
        n <- sum(tbl)
        if (n > 0 && all(dim(tbl) > 1)) {
          sqrt(test$statistic / (n * (min(dim(tbl)) - 1)))
        } else {
          NA_real_
        }
      } else {
        NA_real_
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(pval)) %>%
    rstatix::add_significance(p.col = "pval", output.col = "pval.signif") %>%
    rstatix::adjust_pvalue(p.col = "pval", output.col = "pval.adj", method = "fdr") %>%
    rstatix::add_significance(p.col = "pval.adj", output.col = "pval.adj.signif") %>%
    dplyr::mutate(
      Test = "Chi Squared",
      logp = -log10(pval),
      logp_FDR = -log10(pval.adj),
      `p<.05` = factor(pval.signif, levels = c("ns","*","**","***","****")),
      pval.adj.signif = factor(pval.adj.signif, levels = c("ns","*","**","***","****"))
    )

  if (nrow(stat.test) == 0) return(empty_return())

  if (Relabel) {

    Data <- ReplaceMissingLabels(Data)

    xlab <- sjlabelled::get_label(Data[as.character(stat.test$XVar)], def.value = stat.test$XVar) %>%
      as.data.frame() %>%
      tibble::rownames_to_column()
    colnames(xlab) <- c("Variable", "label")

    ylab <- sjlabelled::get_label(Data[as.character(stat.test$YVar)], def.value = stat.test$YVar) %>%
      as.data.frame() %>%
      tibble::rownames_to_column()
    colnames(ylab) <- c("Variable", "label")

    stat.test$XLabel <- xlab$label
    stat.test$YLabel <- ylab$label

  } else {

    stat.test$XLabel <- stat.test$XVar
    stat.test$YLabel <- stat.test$YVar
  }

  PlotText <- paste(
    "</br>X:", stat.test$XVar,
    "</br>Y:", stat.test$YVar,
    "</br>p:", stat.test$pval, stat.test$pval.signif,
    "</br>FDR p:", stat.test$pval.adj, stat.test$pval.adj.signif,
    "</br>Cramer's V:", round(stat.test$cramers_v, 3)
  )

  p_g <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = `p<.05`, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = pmin(logp, 5), colour = cramers_v)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction")

  p_g_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = pval.adj.signif, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = pmin(logp_FDR, 5), colour = cramers_v)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction")

  pvaltable <- stat.test %>%
    dplyr::select(YVar, XVar, pval) %>%
    tidyr::pivot_wider(names_from = YVar, values_from = pval)

  pvaltable_FDR <- stat.test %>%
    dplyr::select(YVar, XVar, pval.adj) %>%
    tidyr::pivot_wider(names_from = YVar, values_from = pval.adj)

  list(p = p_g, pvaltable = pvaltable, p_FDR = p_g_FDR, pvaltable_FDR = pvaltable_FDR)
}
