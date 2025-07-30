
#'This function creates directional heatmaps for the associations between continuous and binary variables.
#'It identifies continuous and binary variables in the dataset and computes correlations using appropriate methods.
#' @param Data A dataframe containing the data to analyze.
#'@param xVars A vector of variable names to be used on the x-axis. If NULL, all continuous and binary variables are included.
#'@param yVars A vector of variable names to be used on the y-axis. If NULL, defaults to the same as xVars.
#'@param Relabel Logical. If TRUE, variables are relabeled for better readability in the plots.
#'@param Ordinal Logical. If TRUE, treats binary variables as ordinal for certain calculations (future functionality).
#'@return A list containing two plots (Unadjusted and FDR Corrected), along with corresponding data.
#' @export
PlotDirectionalHeatmaps <- function (Data, xVars = NULL, yVars = NULL,  Relabel = TRUE, Ordinal = TRUE) {

  if (is.null(xVars)) {
    ContVars <- getNumVars(Data)
    CatVars <- getBinaryVars(Data)
    xVars <- colnames(Data)[colnames(Data) %in% c(CatVars, ContVars)]
  }
  #xVars <- rev(xVars)
  if (is.null(yVars)) {
    yVars <- xVars
  }

  CatVars <- getBinaryVars(Data %>% select(all_of(unique(c(xVars, yVars)))))
  ContVars <- getNumVars(Data %>% select(all_of(unique(c(xVars, yVars)))))

  xVars <- xVars[xVars %in% c(CatVars, ContVars)]
  yVars <- yVars[yVars %in% c(CatVars, ContVars)]

  df_Combined <- tibble(XVar = character(), YVar = character(), correlation = numeric(), test = character())

  if (length(ContVars) > 0) {
    O_ContCont <- PlotCorrelationsHeatmap(Data, ContVars, Relabel = Relabel)
    df_ContCont <- O_ContCont$Unadjusted$plot$data %>% mutate(correlation = R)
    df_ContCont$test <- O_ContCont$method
    df_Combined <- suppressMessages(full_join(df_Combined, df_ContCont))
  }

  if (length(CatVars) > 0) {
    O_CatCat <- PlotPhiHeatmap(Data, CatVars, Relabel = Relabel)
    df_CatCat <- O_CatCat$Unadjusted$plot$data %>% mutate(correlation = Phi)
    df_CatCat$test <- "Phi"
    df_Combined <- suppressMessages(full_join(df_Combined, df_CatCat))
  }

  if (length(CatVars) > 0 & length(ContVars) > 0) {
    O_CatCont <- PlotPointCorrelationsHeatmap(Data, CatVars, ContVars, Relabel = Relabel)
    df_CatCont <- O_CatCont$Unadjusted$plot$data %>%
      mutate(stars = `p<.05`, stars_FDR = p.adj.signif) %>%
      select(-`p<.05`, -p.adj.signif)

    df_CatCont <- df_CatCont %>% mutate(xVar = CategoricalVariable, yVar = ContinuousVariable)

    df_CatCont1 <- df_CatCont %>%
      select(-CategoricalVariable, -ContinuousVariable) %>%
      mutate(XVar = xVar, YVar = yVar) %>% select(-xVar, -yVar)
    df_CatCont2 <- df_CatCont %>% mutate(XVar = yVar, YVar = xVar) %>%
      select(-xVar, -yVar) %>%
      mutate(XL = XLabel, YL = YLabel) %>%
      mutate(YLabel = XL, XLabel = YL) %>% select(-CategoricalVariable, -ContinuousVariable, -XL, -YL)

    df_CatContSquare <- rbind(df_CatCont1, df_CatCont2)
    df_CatContSquare$test <- "Point Correlation"
    df_Combined <- suppressMessages(full_join(df_Combined, df_CatContSquare))
  }

  df_Combined_plot <- df_Combined %>% filter(XVar %in% xVars) %>% filter(YVar %in% yVars)
  #df_Combined_plot$XLabel <- factor(df_Combined_plot$XLabel, levels = unique(df_Combined_plot$XLabel[df_Combined_plot$XVar %in% xVars]))
  #df_Combined_plot$YLabel <- factor(df_Combined_plot$YLabel, levels = unique(df_Combined_plot$YLabel[df_Combined_plot$YVar %in% yVars]))
  ordered_xlabels <- sapply(xVars, function(var) {
    # take the first label that corresponds to each xVar
    df_Combined_plot$XLabel[df_Combined_plot$XVar == var][1]
  })
  ordered_xlabels <- ordered_xlabels[xVars]
  ordered_xlabels <- unique(ordered_xlabels)
  ordered_ylabels <- sapply(yVars, function(v) {
    df_Combined_plot$YLabel[df_Combined_plot$YVar == v][1]
  })
  ordered_ylabels <- unique(ordered_ylabels)
  #ordered_ylabels<- ordered_ylabels[yVars]
  # 3) Now use that vector (with unique()) to set your factor levels:
  df_Combined_plot$XLabel <- factor(
    df_Combined_plot$XLabel,
    levels = ordered_xlabels
  )
  df_Combined_plot$YLabel <- factor(
    df_Combined_plot$YLabel,
    levels = ordered_ylabels
  )
  p <- ggplot()

  if (nrow(df_Combined_plot %>% filter(test == "pearson"))) {
    p <- p + geom_tile(data = df_Combined_plot %>% filter(test == "pearson"),
                       aes(x = XLabel, y = YLabel, fill = correlation)) +
      scale_fill_gradient2(limits = c(-1, 1),
                           name = "r") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(df_Combined_plot %>% filter(test == "spearman"))) {
    p <- p + geom_tile(data = df_Combined_plot %>% filter(test == "spearman"),
                       aes(x = XLabel, y = YLabel, fill = correlation)) +
      scale_fill_gradient2(limits = c(-1, 1),
                           name = "\u03C1") +  # Unicode for "ρ"
      ggnewscale::new_scale_fill()
  }

  if (nrow(df_Combined_plot %>% filter(test == "Phi"))) {
    p <- p + geom_tile(data = df_Combined_plot %>% filter(test == "Phi"),
                       aes(x = XLabel, y = YLabel, fill = correlation)) +
      scale_fill_gradient2(limits = c(-1, 1),
                           name = "\u03A6",  # Unicode for "Φ"
                           low = scales::muted("purple"),
                           high = scales::muted("green")) +
      ggnewscale::new_scale_fill()

  }

  if (nrow(df_Combined_plot %>% filter(test == "Point Correlation"))) {
    p <- p + geom_tile(data = df_Combined_plot %>% filter(test == "Point Correlation"),
                       aes(x = XLabel, y = YLabel, fill = correlation)) +
      scale_fill_gradient2(limits = c(-1, 1),
                           name = expression(r[pb]),
                           low = scales::muted("#FFA500"), high = scales::muted("#008080")) +
      ggnewscale::new_scale_fill()
  }

  # reset limits
  p<- p + scale_x_discrete(limits = levels(df_Combined_plot$XLabel)) +
    scale_y_discrete(limits = levels(df_Combined_plot$YLabel))

  p_raw <- p + geom_text(data = df_Combined_plot, aes(x = XLabel, y = YLabel, label = stars), color = "black")
  p_raw <- p_raw +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  p_FDR <- p + geom_text(data = df_Combined_plot, aes(x = XLabel, y = YLabel, label = stars_FDR), color = "black")
  p_FDR <- p_FDR +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  Unadjusted <- list()
  Unadjusted$Relabel <- Relabel
  Unadjusted$data <- df_Combined_plot
  Unadjusted$plot <- p_raw

  FDRCorrected <- Unadjusted
  FDRCorrected$plot <- p_FDR
  BinaryMapping <- createBinaryMapping(Data, CatVars)
  return(list(Unadjusted = Unadjusted, FDRCorrected = FDRCorrected, Relabel = Relabel, BinaryMapping = BinaryMapping))
}
