#' Create PCA Table and Visualization
#'
#' Perform principal component analysis (PCA) on specified variables and create visualizations.
#'
#' @param Data The dataset containing the variables for PCA.
#' @param VarsToReduce A character vector specifying the variables to include in the PCA.
#' @param minThresh The minimum threshold for cumulative proportion of variance (default is 0.85).
#' @param scale Logical, indicating whether to scale the data (default is TRUE).
#' @param VariableCategories, categorical, annotates the lollipop plot with color
#' @param center Logical, indicating whether to center the data (default is TRUE).
#' @return A list containing PCA results and visualizations.
#' @import paletteer
#' @importFrom psych principal fa.sort
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_line geom_point geom_col geom_hline geom_segment geom_point scale_color_manual facet_wrap theme element_text
#' @importFrom xtable xtable
#' @export
CreatePCATable<- function (Data, VarsToReduce, VariableCategories = NULL, minThresh = 0.85, scale = TRUE,
          center = TRUE, Relabel = T, Ordinal = FALSE, numComponents = NULL)
{
  classcolors <- c(paletteer::paletteer_d("calecopal::superbloom2"),
                   paletteer::paletteer_d("calecopal::vermillion"), paletteer::paletteer_d("fishualize::Antennarius_commerson"),
                   paletteer::paletteer_d("fishualize::Bodianus_rufus"))

  if (Relabel == TRUE) {
    Data <- ReplaceMissingLabels(Data)
    LabelCodebook <- data.frame(Variable = colnames(Data),
                                Labels = sjlabelled::get_label(Data))
  }else {
    LabelCodebook <- data.frame(Variable = colnames(Data),
                                Labels = colnames(Data))
  }
  DataSubset <- Data[VarsToReduce]
  if (sum(is.na(DataSubset)) > 0) {
    set.seed(123456)
    colnames(DataSubset) <- make.names(VarsToReduce)
    DataSubset <- missRanger::missRanger(DataSubset, num.trees = 100)
    colnames(DataSubset) <- VarsToReduce
  }

  n <- length(VarsToReduce)

  fit1 <- psych::principal(scale(DataSubset, center = center,
                                 scale = center), nfactors = n, rotate = "none")
  if(is.null(numComponents)){
  nc <- min(which(fit1$Vaccounted[3, ] > minThresh))}else{
    nc <- numComponents
  }
  Vaccounted <- as.data.frame(t(fit1$Vaccounted))
  Vaccounted$Component <- factor(rownames(Vaccounted), levels = rownames(Vaccounted))
  p_scree <- ggplot(Vaccounted, aes(x = Component)) + geom_line(aes(y = `Cumulative Var`,
                                                                    group = 1)) + geom_point(aes(y = `Cumulative Proportion`)) +
    geom_col(aes(y = `Proportion Var`)) + geom_hline(aes(yintercept = minThresh),
                                                     linetype = "dashed") + theme_bw()
  fit <- psych::principal(DataSubset, nfactors = nc, rotate = "varimax")
  fit <- psych::fa.sort(fit)
  LoadingTable <- as.data.frame(xtable(unclass(fit$loadings)))
  scores <- fit$scores
  CombinedData <- cbind(Data, scores)
  LoadingTable$Variable <- rownames(LoadingTable)
  LoadingTable <- left_join(LoadingTable, LabelCodebook, by = "Variable", multiple = "first")
  meltedLoading <- LoadingTable %>% tidyr::pivot_longer(cols = matches("RC"))
  meltedLoading$name <- factor(meltedLoading$name, levels = colnames(LoadingTable)[1:length(LoadingTable) -
                                                                                     1])
  meltedLoading <- meltedLoading %>% group_by(Variable) %>%
    arrange(value) %>% mutate(color = abs(value) > 0.4)
  meltedLoading$TopContributor <- meltedLoading$color

  if (is.null(VariableCategories)) {
    meltedLoading$color <- "black"
  }else {
    meltedLoading$color <- VariableCategories[match(meltedLoading$Variable,
                                                    VarsToReduce)]
  }
  meltedLoading$Labels <- factor(meltedLoading$Labels, levels = rev(LoadingTable$Labels))
  if(is.null(VariableCategories)){
    p <- ggplot(meltedLoading, aes(x = Labels, y = value, alpha = TopContributor,
                                   color = color)) + coord_flip()
    p <- p + geom_segment(size = 2, aes(x = Labels, xend = Labels,
                                        y = 0, yend = value))
    p <- p + facet_wrap(vars(name), nrow = 1)
    p <- p + theme(legend.position = "none", strip.text.y = element_text(angle = 45),
                   axis.title.x = element_blank(), axis.title.y = element_blank())
    p <- p + geom_point() + scale_color_manual(values = c("black"))} else{
      p <- ggplot(meltedLoading, aes(x = Labels, y = value, alpha = TopContributor,
                                     color = color)) + coord_flip()
      p <- p + geom_segment(size = 2, aes(x = Labels, xend = Labels,
                                          y = 0, yend = value))
      p <- p + facet_wrap(vars(name), nrow = 1)
      p <- p + theme(legend.position = "none", strip.text.y = element_text(angle = 45),
                     axis.title.x = element_blank(), axis.title.y = element_blank())
      p <- p + geom_point() + scale_color_manual(values = classcolors)
    }
  return(list(p_scree = p_scree, pcaresults = fit, LoadingTable = LoadingTable,
              Scores = scores, CombinedData = CombinedData, Lollipop = p))
}
