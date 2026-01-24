#' PlotMiningMatrix
#'
#' Generate a matrix of statistical relationships between outcome and predictor variables.
#' Returns unadjusted and FDR-adjusted results plus corresponding dot-matrix plots.
#'
#' @param Data A data frame.
#' @param OutcomeVars Outcome variables (typically numeric, e.g., biomarker-by-batch columns).
#' @param PredictorVars Predictor variables (covariates to scan).
#' @param Covariates Optional adjustment covariates.
#' @param Relabel Whether to use variable labels instead of names.
#' @param Parametric If TRUE uses Pearson; if FALSE uses Spearman and non-parametric ANOVA where relevant.
#'
#' @return List with Unadjusted (PvalTable, plot), FDRCorrected (PvalTable, plot),
#'   method, Relabel, Covariates.
#'
#' @importFrom dplyr select mutate filter bind_rows coalesce any_of all_of
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#' @importFrom stats p.adjust
#' @importFrom sjlabelled get_label
#' @importFrom rstatix add_significance
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_manual
#' @importFrom ggplot2 labs xlab ylab theme_bw theme element_text element_blank geom_blank theme_void
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

  out_synt  <- make.names(OutcomeVars)
  pred_synt <- make.names(PredictorVars)
  map_out   <- stats::setNames(OutcomeVars, out_synt)
  map_pred  <- stats::setNames(PredictorVars, pred_synt)

  map_back <- function(x, map) {
    x <- as.character(x)
    idx <- !is.na(x) & x %in% names(map)
    x[idx] <- unname(map[x[idx]])
    x
  }

  col_or_na <- function(df, nm, type = c("numeric", "character")) {
    type <- match.arg(type)
    if (!nm %in% names(df)) {
      if (type == "numeric") return(rep(NA_real_, nrow(df)))
      return(rep(NA_character_, nrow(df)))
    }
    if (type == "numeric") return(suppressWarnings(as.numeric(df[[nm]])))
    as.character(df[[nm]])
  }

  coalesce_p <- function(df) {
    out <- rep(NA_real_, nrow(df))
    for (nm in c("p", "P", "pval", "p.value", "p_value")) {
      if (nm %in% names(df)) {
        out <- dplyr::coalesce(out, suppressWarnings(as.numeric(df[[nm]])))
      }
    }
    out
  }

  ok_df <- function(x) is.data.frame(x) && nrow(x) > 0
  safe_call <- function(expr) tryCatch(expr, error = function(e) NULL)

  normalize_p_cols <- function(df) {
    if (!ok_df(df)) return(df)
    if (!"p" %in% names(df)) {
      if ("P" %in% names(df)) df$p <- df$P
      else if ("pval" %in% names(df)) df$p <- df$pval
      else if ("p.value" %in% names(df)) df$p <- df[["p.value"]]
      else if ("p_value" %in% names(df)) df$p <- df[["p_value"]]
    }
    df
  }

  empty_return <- function() {
    empty <- data.frame()
    p0 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_blank() + ggplot2::theme_void()
    list(
      Unadjusted   = list(PvalTable = empty, plot = p0),
      FDRCorrected = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    )
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) return(empty_return())

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

      orig_XVar   <- col_or_na(cor_df, "XVar", "character")
      orig_YVar   <- col_or_na(cor_df, "YVar", "character")
      orig_XLabel <- dplyr::coalesce(col_or_na(cor_df, "XLabel", "character"), orig_XVar)
      orig_YLabel <- dplyr::coalesce(col_or_na(cor_df, "YLabel", "character"), orig_YVar)

      P_Correlations <- cor_df
      P_Correlations$p    <- coalesce_p(cor_df)
      P_Correlations$Test <- method

      P_Correlations$XVar   <- orig_YVar
      P_Correlations$YVar   <- orig_XVar
      P_Correlations$XLabel <- orig_YLabel
      P_Correlations$YLabel <- orig_XLabel

      P_Correlations <- P_Correlations %>% dplyr::select(-dplyr::any_of("P_adj"))
      P_Correlations <- normalize_p_cols(P_Correlations)
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
      tab <- normalize_p_cols(tmp$Unadjusted$PvalTable)

      P_Anova1 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) .data$DFd else NA_real_,
          XVar   = trimws(as.character(.data$ContinuousVariable)),
          YVar   = trimws(as.character(.data$CategoricalVariable)),
          XLabel = dplyr::coalesce(trimws(as.character(.data$YLabel)), trimws(as.character(.data$ContinuousVariable))),
          YLabel = dplyr::coalesce(trimws(as.character(.data$XLabel)), trimws(as.character(.data$CategoricalVariable))),
          Test   = "ANOVA"
        )
      P_Anova1 <- normalize_p_cols(P_Anova1)
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
      tab <- normalize_p_cols(tmp$Unadjusted$PvalTable)

      P_Anova2 <- tab %>%
        dplyr::select(-dplyr::any_of(c("p.adj","p.adj.signif","logp_FDR"))) %>%
        dplyr::mutate(
          nPairs = if ("DFd" %in% names(.)) .data$DFd else NA_real_,
          XVar   = trimws(as.character(.data$CategoricalVariable)),
          YVar   = trimws(as.character(.data$ContinuousVariable)),
          XLabel = dplyr::coalesce(trimws(as.character(.data$XLabel)), trimws(as.character(.data$CategoricalVariable))),
          YLabel = dplyr::coalesce(trimws(as.character(.data$YLabel)), trimws(as.character(.data$ContinuousVariable))),
          Test   = "ANOVA"
        )
      P_Anova2 <- normalize_p_cols(P_Anova2)
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
      if ("pval" %in% names(chi_df)) chi_df$p <- chi_df$pval
      if ("pval.adj" %in% names(chi_df)) chi_df$p_adj <- chi_df$pval.adj

      P_Chi <- chi_df %>%
        dplyr::mutate(
          Test   = "Chi Squared",
          XVar   = trimws(as.character(.data$XVar)),
          YVar   = trimws(as.character(.data$YVar)),
          XLabel = dplyr::coalesce(trimws(as.character(.data$XLabel)), trimws(as.character(.data$XVar))),
          YLabel = dplyr::coalesce(trimws(as.character(.data$YLabel)), trimws(as.character(.data$YVar)))
        )
      P_Chi <- normalize_p_cols(P_Chi)
    }
  }

  non_null <- purrr::keep(list(P_Correlations, P_Anova1, P_Anova2, P_Chi), ok_df)
  if (length(non_null) == 0) return(empty_return())

  P_Combined <- dplyr::bind_rows(non_null) %>% normalize_p_cols()

  for (nm in c("XVar","YVar","XLabel","YLabel")) {
    if (!nm %in% names(P_Combined)) P_Combined[[nm]] <- NA_character_
  }

  P_Combined <- P_Combined %>%
    dplyr::mutate(
      XVar   = trimws(as.character(.data$XVar)),
      YVar   = trimws(as.character(.data$YVar)),
      XLabel = trimws(as.character(.data$XLabel)),
      YLabel = trimws(as.character(.data$YLabel))
    )

  if (!Relabel) {
    P_Combined <- P_Combined %>%
      dplyr::mutate(
        XVar   = map_back(.data$XVar, map_out),
        YVar   = map_back(.data$YVar, map_pred),
        XLabel = .data$XVar,
        YLabel = .data$YVar
      )
  } else {
    P_Combined <- P_Combined %>%
      dplyr::mutate(
        XLabel = dplyr::coalesce(.data$XLabel, .data$XVar),
        YLabel = dplyr::coalesce(.data$YLabel, .data$YVar)
      )
  }

  P_Combined <- P_Combined %>%
    dplyr::filter(
      !is.na(.data$XVar), !is.na(.data$YVar),
      .data$XVar != "NA", .data$YVar != "NA",
      !is.na(.data$p)
    )

  if (!ok_df(P_Combined)) return(empty_return())

  P_Combined <- P_Combined %>%
    dplyr::filter(.data$XVar %in% OutcomeVars, .data$YVar %in% PredictorVars)

  if (!ok_df(P_Combined)) return(empty_return())

  if (Relabel) {
    xlab_map <- sjlabelled::get_label(Data[, OutcomeVars, drop = FALSE], def.value = OutcomeVars)
    names(xlab_map) <- OutcomeVars
    ylab_map <- sjlabelled::get_label(Data[, PredictorVars, drop = FALSE], def.value = PredictorVars)
    names(ylab_map) <- PredictorVars

    P_Combined <- P_Combined %>%
      dplyr::mutate(
        XLabel = unname(xlab_map[.data$XVar]),
        YLabel = unname(ylab_map[.data$YVar])
      )
    xorder <- unname(xlab_map[OutcomeVars])
    yorder <- unname(ylab_map[PredictorVars])
  } else {
    P_Combined <- P_Combined %>%
      dplyr::mutate(
        XLabel = .data$XVar,
        YLabel = .data$YVar
      )
    xorder <- OutcomeVars
    yorder <- PredictorVars
  }

  P_Combined <- P_Combined %>%
    dplyr::filter(!is.na(.data$XLabel), !is.na(.data$YLabel))

  if (!ok_df(P_Combined)) return(empty_return())

  P_Combined <- P_Combined %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "fdr")) %>%
    rstatix::add_significance(p.col = "p", output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr") %>%
    dplyr::mutate(
      stars     = ifelse(stars     == "****", "***", as.character(stars)),
      stars_fdr = ifelse(stars_fdr == "****", "***", as.character(stars_fdr))
    )

  sig_levels <- c("ns","*","**","***","****")

  P_Combined <- P_Combined %>%
    dplyr::mutate(
      XLabel = factor(.data$XLabel, levels = unique(as.character(xorder))),
      YLabel = factor(.data$YLabel, levels = unique(as.character(yorder))),
      stars = factor(.data$stars, levels = sig_levels),
      stars_fdr = factor(.data$stars_fdr, levels = sig_levels)
    )

  size_num_map <- c("ns" = 1, "*" = 2, "**" = 3, "***" = 4, "****" = 5)
  P_Combined$stars_num     <- unname(size_num_map[as.character(P_Combined$stars)])
  P_Combined$stars_fdr_num <- unname(size_num_map[as.character(P_Combined$stars_fdr)])

  p <- ggplot2::ggplot(P_Combined, ggplot2::aes(y = YLabel, x = XLabel)) +
    ggplot2::geom_point(ggplot2::aes(size = stars_num, colour = stars), shape = 18) +
    ggplot2::scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
    ggplot2::scale_color_manual(
      values = paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5, direction = -1),
      drop = FALSE
    ) +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_blank()
    )

  p_FDR <- ggplot2::ggplot(P_Combined, ggplot2::aes(y = YLabel, x = XLabel)) +
    ggplot2::geom_point(ggplot2::aes(size = stars_fdr_num, colour = stars_fdr), shape = 18) +
    ggplot2::scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
    ggplot2::scale_color_manual(
      values = paletteer::paletteer_c("grDevices::Purple-Yellow", n = 5, direction = -1),
      drop = FALSE
    ) +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_blank()
    )

  list(
    Unadjusted   = list(PvalTable = P_Combined, plot = p),
    FDRCorrected = list(PvalTable = P_Combined, plot = p_FDR),
    method = method, Relabel = Relabel, Covariates = Covariates
  )
}
