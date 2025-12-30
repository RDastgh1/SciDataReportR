#' PlotMiningMatrix
#'
#' Generates a matrix of statistical relationships between specified outcome and predictor variables.
#' Produces both unadjusted and FDR-corrected significance plots.
#'
#' @param Data A data frame containing the dataset to analyze.
#' @param OutcomeVars A vector of outcome variables to be analyzed.
#' @param PredictorVars A vector of predictor variables to be analyzed.
#' @param Covariates Optional vector of covariates to adjust the analysis (default NULL).
#' @param Relabel Logical; relabel variables in output (default TRUE).
#' @param Parametric Logical; if TRUE uses Pearson / parametric ANOVA, else Spearman / nonparametric options where supported.
#'
#' @return A list with Unadjusted and FDRCorrected tables/plots plus metadata.
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
      Unadjusted = list(PvalTable = empty, plot = p0),
      FDRCorrected = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    )
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) return(empty_return())

  ok_df <- function(x) is.data.frame(x) && nrow(x) > 0

  safe_call <- function(expr) {
    tryCatch(expr, error = function(e) NULL)
  }

  # Safely create/normalize a `p` column without ever referencing missing columns.
  ensure_p_col <- function(df) {
    if (!ok_df(df)) return(df)

    if (!"p" %in% names(df)) {
      candidates <- c("p", "P", "pval", "p.value", "p_value", "p.value.", "pvalue")
      found <- intersect(candidates, names(df))
      if (length(found) > 0) {
        df$p <- df[[found[1]]]
      } else {
        df$p <- NA_real_
      }
    }
    df
  }

  # split numeric vs categorical using SciDataReportR helpers
  num_Outcomes   <- SciDataReportR::getNumVars(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)))
  cat_Outcomes   <- SciDataReportR::getCatVars(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)))
  num_Predictors <- SciDataReportR::getNumVars(Data %>% dplyr::select(dplyr::all_of(PredictorVars)))
  cat_Predictors <- SciDataReportR::getCatVars(Data %>% dplyr::select(dplyr::all_of(PredictorVars)))

  P_Correlations <- NULL
  if (length(num_Outcomes) > 0 && length(num_Predictors) > 0) {

    tmp <- safe_call(
      SciDataReportR::PlotCorrelationsHeatmap(
        Data,
        xVars   = num_Predictors,
        yVars   = num_Outcomes,
        covars  = Covariates,
        method  = method,
        Relabel = Relabel
      )
    )

    if (!is.null(tmp) && !is.null(tmp$Unadjusted$plot) && ok_df(tmp$Unadjusted$plot$data)) {
      cor_df <- tmp$Unadjusted$plot$data

      # Remove any existing adjusted column so downstream is consistent
      cor_df <- cor_df %>% dplyr::select(-dplyr::any_of("P_adj"))

      cor_df <- ensure_p_col(cor_df)
      if ("p" %in% names(cor_df)) {
        P_Correlations <- cor_df %>%
          dplyr::mutate(Test = method) %>%
          # Some internal plots use X/Y label conventions differently; keep your prior swap behavior
          dplyr::rename(XLabel = YLabel, YLabel = XLabel)

        if (!"XVar" %in% names(P_Correlations)) P_Correlations$XVar <- NA_character_
        if (!"YVar" %in% names(P_Correlations)) P_Correlations$YVar <- NA_character_
      }
    }
  }

  P_Anova1 <- NULL
  if (length(cat_Predictors) > 0 && length(num_Outcomes) > 0) {
    tmp <- safe_call(
      SciDataReportR::PlotAnovaRelationshipsMatrix(
        Data, CatVars = cat_Predictors, ContVars = num_Outcomes,
        Covariates = Covariates, Relabel = Relabel, Parametric = Parametric
      )
    )

    if (!is.null(tmp) && ok_df(tmp$Unadjusted$PvalTable)) {
      tab <- tmp$Unadjusted$PvalTable
      tab <- ensure_p_col(tab)

      P_Anova1 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) .data$DFd else NA_real_,
          XVar = .data$ContinuousVariable,
          YVar = .data$CategoricalVariable,
          Test = "ANOVA"
        ) %>%
        dplyr::rename(XLabel = YLabel, YLabel = XLabel)

      P_Anova1 <- ensure_p_col(P_Anova1)
    }
  }

  P_Anova2 <- NULL
  if (length(cat_Outcomes) > 0 && length(num_Predictors) > 0) {
    tmp <- safe_call(
      SciDataReportR::PlotAnovaRelationshipsMatrix(
        Data, CatVars = cat_Outcomes, ContVars = num_Predictors,
        Covariates = Covariates, Relabel = Relabel, Parametric = Parametric
      )
    )

    if (!is.null(tmp) && ok_df(tmp$Unadjusted$PvalTable)) {
      tab <- tmp$Unadjusted$PvalTable
      tab <- ensure_p_col(tab)

      P_Anova2 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) .data$DFd else NA_real_,
          XVar = .data$CategoricalVariable,
          YVar = .data$ContinuousVariable,
          Test = "ANOVA"
        )

      P_Anova2 <- ensure_p_col(P_Anova2)
    }
  }

  P_Chi <- NULL
  if (length(cat_Outcomes) > 0 && length(cat_Predictors) > 0) {
    tmp <- safe_call(
      SciDataReportR::PlotChiSqCovar(
        Data, xVars = cat_Predictors, yVars = cat_Outcomes,
        covars = Covariates, Relabel = Relabel
      )
    )

    if (!is.null(tmp) && !is.null(tmp$p) && ok_df(tmp$p$data)) {
      chi_df <- tmp$p$data

      # Normalize p if chi_df uses pval
      if (!"p" %in% names(chi_df)) {
        if ("pval" %in% names(chi_df)) chi_df$p <- chi_df[["pval"]]
        if ("P" %in% names(chi_df))    chi_df$p <- chi_df[["P"]]
      }

      P_Chi <- chi_df %>%
        dplyr::mutate(Test = "Chi Squared")

      P_Chi <- ensure_p_col(P_Chi)
      if (!"XVar" %in% names(P_Chi)) P_Chi$XVar <- NA_character_
      if (!"YVar" %in% names(P_Chi)) P_Chi$YVar <- NA_character_
    }
  }

  non_null <- list(P_Correlations, P_Anova1, P_Anova2, P_Chi) %>% purrr::keep(ok_df)
  if (length(non_null) == 0) return(empty_return())

  P_Combined <- dplyr::bind_rows(non_null) %>%
    dplyr::filter(!is.na(.data$XLabel), !is.na(.data$YLabel))

  P_Combined <- ensure_p_col(P_Combined)
  if (!"p" %in% names(P_Combined)) return(empty_return())

  P_Combined <- P_Combined %>%
    dplyr::filter(!is.na(.data$p)) %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "fdr")) %>%
    rstatix::add_significance(p.col = "p", output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr") %>%
    dplyr::mutate(
      stars     = ifelse(stars     == "****", "***", as.character(stars)),
      stars_fdr = ifelse(stars_fdr == "****", "***", as.character(stars_fdr))
    )

  # axis ordering
  if (Relabel) {
    xorder <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(OutcomeVars)), def.value = OutcomeVars)
    yorder <- sjlabelled::get_label(Data %>% dplyr::select(dplyr::all_of(PredictorVars)), def.value = PredictorVars)
  } else {
    xorder <- OutcomeVars
    yorder <- PredictorVars
  }
  xorder <- unique(xorder); yorder <- unique(yorder)

  sig_levels <- c("ns","*","**","***","****")

  P_Combined <- P_Combined %>%
    dplyr::mutate(
      XLabel = factor(.data$XLabel, levels = xorder),
      YLabel = factor(.data$YLabel, levels = yorder),
      stars = factor(.data$stars, levels = sig_levels),
      stars_fdr = factor(.data$stars_fdr, levels = sig_levels),
      logp = -log10(.data$p),
      logpfdr = -log10(.data$p_adj)
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
    Unadjusted = list(PvalTable = P_Combined, plot = p),
    FDRCorrected = list(PvalTable = P_Combined, plot = p_FDR),
    method = method, Relabel = Relabel, Covariates = Covariates
  )
}
