#' Create directional heatmaps across continuous & binary variables
#'
#' Combines:
#'  - Continuous~Continuous (Pearson/Spearman)
#'  - Binary~Binary (Phi; 1 == PositiveLevel)
#'  - Binary~Continuous (r_pb; 1 == PositiveLevel)
#' into a single square heatmap with raw and FDR-star overlays.
#'
#' Constant variables (no variation in the current data) are
#' automatically excluded before computing any tiles.
#'
#' @param Data A dataframe.
#' @param xVars Character vector of variables for x-axis (subset of Data cols).
#'              If NULL, uses all detected continuous + binary vars.
#' @param yVars Character vector for y-axis (defaults to xVars if NULL).
#' @param Relabel Logical; use sjlabelled variable labels if present.
#' @param Ordinal Logical; reserved for future use.
#' @return list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)
#' @export
PlotDirectionalHeatmaps <- function (Data, xVars = NULL, yVars = NULL, Relabel = TRUE, Ordinal = TRUE) {

  # ---- helpers -------------------------------------------------------------
  uniq2 <- function(x) length(unique(stats::na.omit(x)))
  has_two_levels <- function(v) uniq2(Data[[v]]) == 2
  has_variance   <- function(v) {
    x <- Data[[v]]
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    length(x) >= 2 && stats::var(x) > 0
  }

  # ---- 1) Determine candidate vars ----------------------------------------
  if (is.null(xVars)) {
    Cont_all <- getNumVars(Data)
    Cat_all  <- getBinaryVars(Data)
    xVars <- colnames(Data)[colnames(Data) %in% c(Cat_all, Cont_all)]
  }
  if (is.null(yVars)) yVars <- xVars

  cand <- unique(c(xVars, yVars))

  # initial typing from your helpers, but we'll enforce non-constant below
  Cat0  <- intersect(getBinaryVars(dplyr::select(Data, dplyr::all_of(cand))), cand)
  Cont0 <- intersect(getNumVars(   dplyr::select(Data, dplyr::all_of(cand))), cand)

  # ---- 2) Drop constants BEFORE any subcalls ------------------------------
  CatVars  <- Cat0[vapply(Cat0, has_two_levels, logical(1))]
  ContVars <- Cont0[vapply(Cont0, has_variance,   logical(1))]

  dropped <- setdiff(cand, c(CatVars, ContVars))
  if (length(dropped)) {
    reasons <- vapply(dropped, function(v) {
      if (v %in% Cont0 && !has_variance(v))         return("no variance")
      if (v %in% Cat0 && !has_two_levels(v))        return("not two levels")
      if (!(v %in% c(Cat0, Cont0)))                 return("not binary/continuous")
      "excluded"
    }, character(1))
    Excluded <- tibble::tibble(Variable = dropped, Reason = reasons)
  } else {
    Excluded <- tibble::tibble(Variable = character(), Reason = character())
  }

  # refresh axes after exclusions
  xVars <- xVars[xVars %in% c(CatVars, ContVars)]
  yVars <- yVars[yVars %in% c(CatVars, ContVars)]

  # Build a single binary mapping (stable 1 == PositiveLevel)
  BinaryMapping <- NULL
  if (length(CatVars) > 0) {
    BinaryMapping <- createBinaryMapping(Data, CatVars)
  }

  # ---- 3) Collect tiles from each sub-plotter -----------------------------
  df_Combined <- tibble::tibble(XVar = character(), YVar = character(),
                                correlation = numeric(), test = character())

  # Continuous ~ Continuous
  if (length(ContVars) > 0) {
    O_ContCont <- PlotCorrelationsHeatmap(Data, ContVars, Relabel = Relabel)
    df_ContCont <- O_ContCont$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = R, test = O_ContCont$method)
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_ContCont))
  }

  # Binary ~ Binary (Phi)
  if (length(CatVars) > 1) {  # need at least two
    O_CatCat <- PlotPhiHeatmap(Data, CatVars, Relabel = Relabel, binary_map = BinaryMapping)
    df_CatCat <- O_CatCat$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = Phi, test = "Phi")
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_CatCat))
  }

  # Binary ~ Continuous (point-biserial)
  if (length(CatVars) > 0 & length(ContVars) > 0) {
    O_CatCont <- PlotPointCorrelationsHeatmap(
      Data, CatVars, ContVars, Relabel = Relabel, Ordinal = Ordinal, binary_map = BinaryMapping
    )

    df_CatCont <- O_CatCont$Unadjusted$plot$data %>%
      dplyr::mutate(stars = `p<.05`, stars_FDR = p.adj.signif) %>%
      dplyr::select(-`p<.05`, -p.adj.signif) %>%
      dplyr::mutate(xVar = CategoricalVariable, yVar = ContinuousVariable)

    # make square by duplicating with swapped axes/labels
    df_CatCont1 <- df_CatCont %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable) %>%
      dplyr::mutate(XVar = xVar, YVar = yVar) %>%
      dplyr::select(-xVar, -yVar)

    df_CatCont2 <- df_CatCont %>%
      dplyr::mutate(XVar = yVar, YVar = xVar) %>%
      dplyr::select(-xVar, -yVar) %>%
      dplyr::mutate(XL = XLabel, YL = YLabel) %>%
      dplyr::mutate(YLabel = XL, XLabel = YL) %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable, -XL, -YL)

    df_CatContSquare <- rbind(df_CatCont1, df_CatCont2)
    df_CatContSquare$test <- "Point Correlation"
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_CatContSquare))
  }

  # ---- 4) Filter to requested axes, lock label order ----------------------
  df_Combined_plot <- df_Combined %>%
    dplyr::filter(XVar %in% xVars, YVar %in% yVars)

  # handle fully-empty safely
  if (nrow(df_Combined_plot) == 0L) {
    empty_plot <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::ggtitle("No valid variable pairs after excluding constants")
    Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)
    FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)
    return(list(Unadjusted = Unadjusted, FDRCorrected = FDRCorrected,
                Relabel = Relabel, BinaryMapping = BinaryMapping, Excluded = Excluded))
  }

  ordered_xlabels <- sapply(xVars, function(var) df_Combined_plot$XLabel[df_Combined_plot$XVar == var][1])
  ordered_xlabels <- unique(ordered_xlabels[xVars])
  ordered_ylabels <- sapply(yVars, function(var) df_Combined_plot$YLabel[df_Combined_plot$YVar == var][1])
  ordered_ylabels <- unique(ordered_ylabels)

  df_Combined_plot$XLabel <- factor(df_Combined_plot$XLabel, levels = ordered_xlabels)
  df_Combined_plot$YLabel <- factor(df_Combined_plot$YLabel, levels = ordered_ylabels)

  # ---- 5) Build layered heatmap ------------------------------------------
  p <- ggplot2::ggplot()

  if (nrow(dplyr::filter(df_Combined_plot, test == "pearson"))) {
    p <- p +
      ggplot2::geom_tile(data = dplyr::filter(df_Combined_plot, test == "pearson"),
                         ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = "r") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "spearman"))) {
    p <- p +
      ggplot2::geom_tile(data = dplyr::filter(df_Combined_plot, test == "spearman"),
                         ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = "\u03C1") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "Phi"))) {
    p <- p +
      ggplot2::geom_tile(data = dplyr::filter(df_Combined_plot, test == "Phi"),
                         ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1),
                                    name = "\u03A6",
                                    low  = scales::muted("purple"),
                                    high = scales::muted("green")) +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "Point Correlation"))) {
    p <- p +
      ggplot2::geom_tile(data = dplyr::filter(df_Combined_plot, test == "Point Correlation"),
                         ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1),
                                    name = expression(r[pb]),
                                    low  = scales::muted("#FFA500"),
                                    high = scales::muted("#008080")) +
      ggnewscale::new_scale_fill()
  }

  # reset limits to preserve ordering
  p <- p +
    ggplot2::scale_x_discrete(limits = levels(df_Combined_plot$XLabel)) +
    ggplot2::scale_y_discrete(limits = levels(df_Combined_plot$YLabel))

  # annotate with raw and FDR stars
  p_raw <- p +
    ggplot2::geom_text(data = df_Combined_plot,
                       ggplot2::aes(x = XLabel, y = YLabel, label = stars),
                       color = "black") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  p_FDR <- p +
    ggplot2::geom_text(data = df_Combined_plot,
                       ggplot2::aes(x = XLabel, y = YLabel, label = stars_FDR),
                       color = "black") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # ---- 6) Return ----------------------------------------------------------
  Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_raw)
  FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_FDR)

  list(Unadjusted   = Unadjusted,
       FDRCorrected = FDRCorrected,
       Relabel      = Relabel,
       BinaryMapping= BinaryMapping,
       Excluded     = Excluded)
}
