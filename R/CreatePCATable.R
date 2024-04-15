#' Create PCA Table and Visualization
#'
#' Perform principal component analysis (PCA) on specified variables and create visualizations.
#'
#' @param Data The dataset containing the variables for PCA.
#' @param VarsToReduce A character vector specifying the variables to include in the PCA.
#' @param minThresh The minimum threshold for cumulative proportion of variance (default is 0.85).
#' @param scale Logical, indicating whether to scale the data (default is TRUE).
#' @param center Logical, indicating whether to center the data (default is TRUE).
#' @return A list containing PCA results and visualizations.
#' @importFrom psych principal fa.sort
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_line geom_point geom_col geom_hline geom_segment geom_point scale_color_manual facet_wrap theme element_text
#' @importFrom xtable xtable
#' @export
CreatePCATable <- function(Data,  VarsToReduce, minThresh = 0.85, scale = TRUE, center = TRUE) {

  # Subset data to include only specified variables
  DataSubset <- Data[VarsToReduce]

  # Check for missing values and handle them if present
  if (sum(is.na(DataSubset)) > 0) {
    set.seed(123456)
    colnames(DataSubset) <- make.names(VarsToReduce)
    DataSubset <- missRanger::missRanger(DataSubset, num.trees = 100)
    colnames(DataSubset) <- VarsToReduce
  }

  # Determine the number of components based on cumulative proportion of variance
  if (length(VarsToReduce) < 20) {
    n <- length(VarsToReduce)
  } else {
    n <- 20
  }

  # Perform PCA
  fit1 <- psych::principal(scale(DataSubset, center = center, scale = center), nfactors = n, rotate = "none")
  nc <- min(which(fit1$Vaccounted[5, ] > minThresh))

  # Create scree plot
  Vaccounted <- as.data.frame(t(fit1$Vaccounted))
  Vaccounted$Component <- factor(rownames(Vaccounted), levels = rownames(Vaccounted))
  p_scree <- ggplot(Vaccounted, aes(x = Component)) +
    geom_line(aes(y = `Cumulative Proportion`, group = 1)) +
    geom_point(aes(y = `Cumulative Proportion`)) +
    geom_col(aes(y = `Proportion Var`)) +
    geom_hline(aes(yintercept = minThresh), linetype = "dashed")

  # Perform PCA with determined number of components and rotate
  fit <- psych::principal(DataSubset, nfactors = nc, rotate = "varimax")
  fit <- fa.sort(fit)
  LoadingTable <- as.data.frame(xtable(unclass(fit$loadings)))
  scores <- fit$scores
  CombinedData <- cbind(Data, scores)

  # Create lollipop figure
  LoadingTable$Variable <- factor(rownames(LoadingTable), levels = rev(rownames(LoadingTable)))
  meltedLoading <- LoadingTable %>% tidyr::pivot_longer(cols = matches("RC"))
  meltedLoading$name <- factor(meltedLoading$name, levels = colnames(LoadingTable)[1:length(LoadingTable) - 1])
  meltedLoading <- meltedLoading %>% group_by(Variable) %>% arrange(value) %>% mutate(color = abs(value) > .4)
  meltedLoading$TopContributor <- meltedLoading$color
  meltedLoading$color <- "black"

  p <- ggplot(meltedLoading, aes(x = Variable, y = value, alpha = TopContributor, color = color)) + coord_flip()
  p <- p + geom_segment(size = 2, aes(x = Variable, xend = Variable, y = 0, yend = value))
  p <- p + facet_wrap(vars(name), nrow = 1)
  p <- p + theme(legend.position = "none", strip.text.y = element_text(angle = 45), axis.title.x = element_blank(), axis.title.y = element_blank())
  p <- p + geom_point() + scale_color_manual(values = c("black"))

  # Return a list of results and visualizations
  return(list(p_scree = p_scree, pcaresults = fit, LoadingTable = LoadingTable, Scores = scores, CombinedData = CombinedData, Lollipop = p))
}
