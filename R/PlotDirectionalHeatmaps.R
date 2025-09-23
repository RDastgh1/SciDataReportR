#' Create directional heatmaps across continuous & binary variables
#'
#' Combines:
#'  - Continuous~Continuous (Pearson/Spearman)
#'  - Binary~Binary (Phi; 1 == PositiveLevel)
#'  - Binary~Continuous (r_pb; 1 == PositiveLevel)
#' into a single square heatmap with raw and FDR-star overlays.
#'
#' @param Data A dataframe.
#' @param xVars Character vector of variables for x-axis (subset of Data cols).
#'              If NULL, uses all detected continuous + binary vars.
#' @param yVars Character vector for y-axis (defaults to xVars if NULL).
#' @param Relabel Logical; use sjlabelled variable labels if present.
#' @param Ordinal Logical; reserved for future use.
#' @return list(Unadjusted, FDRCorrected, Relabel, BinaryMapping)
#' @export
PlotDirectionalHeatmaps <- function (Data, xVars = NULL, yVars = NULL, Relabel = TRUE, Ordinal = TRUE) {

  # 1) Determine candidate vars
  if (is.null(xVars)) {
    ContVars_all <- getNumVars(Data)
    CatVars_all  <- getBinaryVars(Data)
    xVars <- colnames(Data)[colnames(Data) %in% c(CatVars_all, ContVars_all)]
  }
  if (is.null(yVars)) {
    yVars <- xVars
  }

  # restrict to vars that are continuous or binary among x/y
  CatVars <- getBinaryVars(dplyr::select(Data, dplyr::all_of(unique(c(xVars, yVars)))))
  ContVars <- getNumVars(dplyr::select(Data, dplyr::all_of(unique(c(xVars, yVars)))))

  xVars <- xVars[xVars %in% c(CatVars, ContVars)]
  yVars <- yVars[yVars %in% c(CatVars, ContVars)]

  # Build a single binary mapping once and reuse it (stable 1 == PositiveLevel)
  BinaryMapping <- NULL
  if (length(CatVars) > 0) {
    BinaryMapping <- createBinaryMapping(Data, CatVars)
  }

  # 2) Collect tiles from each sub-plotter
  df_Combined <- tibble::tibble(XVar = character(), YVar = character(), correlation = numeric(), test = character())

  # Continuous ~ Continuous
  if (length(ContVars) > 0) {
    O_ContCont <- PlotCorrelationsHeatmap(Data, ContVars, Relabel = Relabel)
    df_ContCont <- O_ContCont$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = R, test = O_ContCont$method)
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_ContCont))
  }

  # Binary ~ Binary (Phi), pass mapping so sign is meaningful
  if (length(CatVars) > 0) {
    O_CatCat <- PlotPhiHeatmap(Data, CatVars, Relabel = Relabel, binary_map = BinaryMapping)
    df_CatCat <- O_CatCat$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = Phi, test = "Phi")
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_CatCat))
  }

  # Binary ~ Continuous (point-biserial), pass mapping so sign is meaningful
  if (length(CatVars) > 0 & length(ContVars) > 0) {
    O_CatCont <- PlotPointCorrelationsHeatmap(Data, CatVars, ContVars, Relabel = Relabel, Ordinal = Ordinal, binary_map = BinaryMapping)

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
    df_CatContSquare$test <- "Point Correlation"  # keep legacy label expected downstream
    df_Combined <- suppressMessages(dplyr::full_join(df_Combined, df_CatContSquare))
  }

  # 3) Filter to requested x/y, and lock label order to xVars/yVars sequence
  df_Combined_plot <- df_Combined %>%
    dplyr::filter(XVar %in% xVars, YVar %in% yVars)

  ordered_xlabels <- sapply(xVars, function(var) df_Combined_plot$XLabel[df_Combined_plot$XVar == var][1])
  ordered_xlabels <- unique(ordered_xlabels[xVars])
  ordered_ylabels <- sapply(yVars, function(var) df_Combined_plot$YLabel[df_Combined_plot$YVar == var][1])
  ordered_ylabels <- unique(ordered_ylabels)

  df_Combined_plot$XLabel <- factor(df_Combined_plot$XLabel, levels = ordered_xlabels)
  df_Combined_plot$YLabel <- factor(df_Combined_plot$YLabel, levels = ordered_ylabels)

  # 4) Build the layered heatmap, preserving separate scales per test
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
                                    low = scales::muted("purple"),
                                    high = scales::muted("green")) +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "Point Correlation"))) {
    p <- p +
      ggplot2::geom_tile(data = dplyr::filter(df_Combined_plot, test == "Point Correlation"),
                         ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1),
                                    name = expression(r[pb]),
                                    low = scales::muted("#FFA500"),
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

  # 5) Return bundle (keep BinaryMapping for transparency/debugging)
  Unadjusted <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_raw)
  FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_FDR)

  return(list(Unadjusted = Unadjusted,
              FDRCorrected = FDRCorrected,
              Relabel = Relabel,
              BinaryMapping = BinaryMapping))
}
