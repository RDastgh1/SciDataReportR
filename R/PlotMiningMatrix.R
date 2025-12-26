#' PlotMiningMatrix
#'
#' This function generates a matrix of statistical relationships between specified outcome and predictor variables
#' in a dataset. It includes visualizations for correlations, ANOVA results, and FDR corrections.
#'
#' @param Data A data frame containing the dataset to analyze.
#' @param OutcomeVars A vector of outcome variables to be analyzed.
#' @param PredictorVars A vector of predictor variables to be analyzed.
#' @param Covariates An optional vector of covariates to adjust the analysis (default is NULL).
#' @param Relabel Logical flag indicating whether to relabel the variables in the output (default is TRUE).
#' @param Parametric Logical flag indicating whether to use parametric methods (default is TRUE).
#' If FALSE, non-parametric methods will be used.
#'
#' @return A list containing the following elements:
#' \item{Unadjusted}{A list with the unadjusted p-value table and corresponding plot.}
#' \item{FDRCorrected}{A list with the FDR-adjusted p-value table and corresponding plot.}
#' \item{method}{The method used for correlation ("pearson" for parametric, "spearman" for non-parametric).}
#' \item{Relabel}{The value of the Relabel parameter.}
#' \item{Covariates}{The covariates used in the analysis, if any.}
#'
#' @import dplyr ggplot2 rstatix sjlabelled paletteer
#' @importFrom stats p.adjust
#' @importFrom magrittr %>%
#' @export
PlotMiningMatrix <- function(
    Data,
    OutcomeVars,
    PredictorVars,
    Covariates = NULL,
    Relabel = TRUE,
    Parametric = TRUE
) {

  method <- ifelse(Parametric, "pearson", "spearman")

  OutcomeVars   <- unique(intersect(as.character(OutcomeVars), names(Data)))
  PredictorVars <- unique(intersect(as.character(PredictorVars), names(Data)))
  Covariates <- if (!is.null(Covariates) && length(Covariates) > 0) {
    unique(intersect(as.character(Covariates), names(Data)))
  } else NULL

  empty_return <- function() {
    empty <- data.frame()
    p0 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_blank() + ggplot2::theme_void()
    list(
      Unadjusted    = list(PvalTable = empty, plot = p0),
      FDRCorrected  = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    )
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) return(empty_return())

  ok_df <- function(x) is.data.frame(x) && nrow(x) > 0

  num_Outcomes   <- getNumVars(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)))
  cat_Outcomes   <- getCatVars(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)))
  num_Predictors <- getNumVars(Data %>% dplyr::select(dplyr::all_of(PredictorVars)))
  cat_Predictors <- getCatVars(Data %>% dplyr::select(dplyr::all_of(PredictorVars)))

  P_Correlations <- NULL
  if (length(num_Outcomes) > 0 && length(num_Predictors) > 0) {
    tmp <- PlotCorrelationsHeatmap(
      Data,
      xVars = num_Predictors,
      yVars = num_Outcomes,
      covars = Covariates,
      method = method,
      Relabel = Relabel
    )

    cor_df <- tmp$Unadjusted$plot$data
    if (ok_df(cor_df)) {
      P_Correlations <- cor_df %>%
        dplyr::select(-dplyr::any_of("P_adj")) %>%
        dplyr::mutate(Test = method, p = P) %>%
        dplyr::rename(XLabel = YLabel, YLabel = XLabel)
    }
  }

  P_Anova1 <- NULL
  if (length(cat_Predictors) > 0 && length(num_Outcomes) > 0) {
    tmp <- PlotAnovaRelationshipsMatrix(
      Data,
      CatVars     = cat_Predictors,
      ContVars    = num_Outcomes,
      Covariates  = Covariates,
      Relabel     = Relabel,
      Parametric  = Parametric
    )
    tab <- tmp$Unadjusted$PvalTable
    if (ok_df(tab)) {
      P_Anova1 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) as.numeric(.data$DFd) else rep(NA_real_, dplyr::n()),
          XVar   = ContinuousVariable,
          YVar   = CategoricalVariable,
          p      = .data$p
        ) %>%
        dplyr::rename(XLabel = YLabel, YLabel = XLabel)
    }
  }

  P_Anova2 <- NULL
  if (length(cat_Outcomes) > 0 && length(num_Predictors) > 0) {
    tmp <- PlotAnovaRelationshipsMatrix(
      Data,
      CatVars     = cat_Outcomes,
      ContVars    = num_Predictors,
      Covariates  = Covariates,
      Relabel     = Relabel,
      Parametric  = Parametric
    )
    tab <- tmp$Unadjusted$PvalTable
    if (ok_df(tab)) {
      P_Anova2 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) as.numeric(.data$DFd) else rep(NA_real_, dplyr::n()),
          XVar   = CategoricalVariable,
          YVar   = ContinuousVariable,
          p      = .data$p
        )
    }
  }

  P_Chi <- NULL
  if (length(cat_Outcomes) > 0 && length(cat_Predictors) > 0) {
    tmp <- PlotChiSqCovar(
      Data,
      xVars  = cat_Predictors,
      yVars  = cat_Outcomes,
      covars = Covariates,
      Relabel = Relabel
    )
    chi_df <- tmp$p$data
    if (ok_df(chi_df)) {
      P_Chi <- chi_df %>%
        dplyr::mutate(
          p = dplyr::coalesce(.data$pval, .data$p, NA_real_),
          Test = "Chi Squared"
        )
      if ("pval.adj" %in% names(chi_df)) P_Chi <- P_Chi %>% dplyr::mutate(p_adj = .data$pval.adj)
    }
  }

  non_null <- list(P_Correlations, P_Anova1, P_Anova2, P_Chi) %>% purrr::keep(ok_df)
  if (length(non_null) == 0) return(empty_return())

  # safer join keys to avoid accidental cartesian explosion
  join_keys <- c("XLabel","YLabel","XVar","YVar")
  join_two <- function(a, b) dplyr::full_join(a, b, by = intersect(join_keys, intersect(names(a), names(b))))
  P_Combined <- Reduce(join_two, non_null)

  out_lab  <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)),   def.value = OutcomeVars)
  pred_lab <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(PredictorVars)), def.value = PredictorVars)

  out_map  <- setNames(OutcomeVars, out_lab)
  pred_map <- setNames(PredictorVars, pred_lab)
  any_map  <- c(out_map, pred_map)

  if (!"XVar" %in% names(P_Combined)) P_Combined$XVar <- NA_character_
  if (!"YVar" %in% names(P_Combined)) P_Combined$YVar <- NA_character_

  P_Combined <- P_Combined %>%
    dplyr::filter(!is.na(p)) %>%
    dplyr::mutate(
      XVar = ifelse(is.na(XVar) & !is.na(XLabel) & as.character(XLabel) %in% names(any_map),
                    any_map[as.character(XLabel)], XVar),
      YVar = ifelse(is.na(YVar) & !is.na(YLabel) & as.character(YLabel) %in% names(any_map),
                    any_map[as.character(YLabel)], YVar),
      p_adj = stats::p.adjust(p, method = "fdr")
    ) %>%
    rstatix::add_significance(p.col = "p",     output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr") %>%
    dplyr::mutate(
      stars     = ifelse(stars     == "****", "***", as.character(stars)),
      stars_fdr = ifelse(stars_fdr == "****", "***", as.character(stars_fdr))
    )

  if (Relabel) {
    xorder <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)),   def.value = OutcomeVars)
    yorder <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(PredictorVars)), def.value = PredictorVars)
  } else {
    xorder <- OutcomeVars
    yorder <- PredictorVars
  }

  sig_levels <- c("ns","*","**","***","****")

  P_Combined <- P_Combined %>%
    dplyr::mutate(
      XLabel    = factor(XLabel, levels = unique(xorder)),
      YLabel    = factor(YLabel, levels = unique(yorder)),
      stars     = factor(stars,     levels = sig_levels),
      stars_fdr = factor(stars_fdr, levels = sig_levels)
    )

  p <- ggplot2::ggplot(P_Combined, ggplot2::aes(y = YLabel, x = XLabel)) +
    ggplot2::geom_point(ggplot2::aes(size = stars, colour = stars), shape = 18) +
    ggplot2::scale_color_manual(values = paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5, direction = -1)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.title = ggplot2::element_blank())

  p_FDR <- ggplot2::ggplot(P_Combined, ggplot2::aes(y = YLabel, x = XLabel)) +
    ggplot2::geom_point(ggplot2::aes(size = stars_fdr, colour = stars_fdr), shape = 18) +
    ggplot2::scale_color_manual(values = paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5, direction = -1)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.title = ggplot2::element_blank())

  list(
    Unadjusted   = list(PvalTable = P_Combined, plot = p),
    FDRCorrected = list(PvalTable = P_Combined, plot = p_FDR),
    method = method, Relabel = Relabel, Covariates = Covariates
  )
}
