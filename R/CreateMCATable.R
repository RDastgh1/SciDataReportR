#' CreateMCATable
#'
#' This function performs Multiple Correspondence Analysis (MCA) on a set of categorical variables, imputes missing data if needed, and generates a set of visualizations and tables to interpret the results.
#'
#' @param Data A dataframe containing the data to be analyzed.
#' @param VarsToReduce A vector of column names in `Data` to be included in the MCA.
#' @param VariableCategories An optional vector to assign specific categories to the variables in `VarsToReduce`. These will be used to color the loadings plot.
#' @param minThresh A numeric value representing the minimum cumulative variance threshold to determine the number of components. Default is 75%.
#' @param scale Logical, indicating whether to scale the variables. Default is TRUE.
#' @param center Logical, indicating whether to center the variables. Default is TRUE.
#' @param Relabel Logical, if TRUE, the function will replace missing labels in the data using an external helper function `ReplaceMissingLabels`. Default is TRUE.
#' @param Ordinal Logical, if TRUE, the function will treat variables as ordinal for MCA. Default is FALSE.
#' @param numComponents An optional integer specifying the number of components to retain. If NULL, the number of components will be determined based on `minThresh`.
#'
#' @return A list with the following elements:
#' \item{p_scree}{A `ggplot` object representing the scree plot, showing the cumulative and percentage of variance explained by each component.}
#' \item{pcaresults}{The MCA results object, which includes component scores and contributions.}
#' \item{LoadingTable}{A data frame with the variable loadings for each component.}
#' \item{Scores}{A data frame with the MCA scores for each individual in the data.}
#' \item{CombinedData}{The original data combined with the MCA scores.}
#' \item{Lollipop}{A `ggplot` object showing a lollipop plot of variable loadings across components.}
#'
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_col geom_segment coord_flip theme_bw facet_wrap geom_hline scale_x_discrete
#' @importFrom paletteer paletteer_d
#' @importFrom missRanger missRanger
#' @importFrom FactoMineR MCA
#' @importFrom sjlabelled get_label
#' @importFrom psych fa.sort
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select group_by arrange mutate
#' @importFrom tibble rownames_to_column
#' @importFrom scales removeString
#'
#' @export
#'
CreateMCATable <- function(Data, VarsToReduce, VariableCategories = NULL,
                           minThresh = 75, scale = TRUE, center = TRUE,
                           Relabel = TRUE, Ordinal = FALSE,
                           numComponents = NULL) {



  # Define custom colors for different classes
  classcolors <- c(paletteer::paletteer_d("calecopal::superbloom2"),
                   paletteer::paletteer_d("calecopal::vermillion"),
                   paletteer::paletteer_d("fishualize::Antennarius_commerson"),
                   paletteer::paletteer_d("fishualize::Bodianus_rufus"))

  # Optionally relabel data using an external helper function ReplaceMissingLabels
  if (Relabel == TRUE) {
    Data <- ReplaceMissingLabels(Data)
    LabelCodebook <- data.frame(Variable = colnames(Data),
                                Labels = sjlabelled::get_label(Data))
  } else {
    LabelCodebook <- data.frame(Variable = colnames(Data),
                                Labels = colnames(Data))
  }

  # Subset the data to the selected variables for reduction
  DataSubset <- Data[VarsToReduce]
  DataSubset[] <- lapply(DataSubset, factor)  # Ensure all variables are factors

  # Handle missing data using missRanger for imputation if necessary
  if (sum(is.na(DataSubset)) > 0) {
    set.seed(123456)
    colnames(DataSubset) <- make.names(VarsToReduce)  # Ensure valid column names
    DataSubset <- missRanger::missRanger(DataSubset, num.trees = 100)
    colnames(DataSubset) <- VarsToReduce  # Restore original column names
  }

  # Perform MCA with all components, select appropriate number of components
  n <- length(VarsToReduce)
  mca_result <- MCA(DataSubset %>% select(all_of(VarsToReduce)),
                    ncp = n, graph = FALSE)

  # Determine number of components if not provided
  if (is.null(numComponents)) {
    nc <- min(which(mca_result$eig[, 3] >= minThresh))  # Select based on variance
  } else {
    nc <- numComponents
  }

  # Create a scree plot showing the percentage of variance explained by each component
  Vaccounted <- as.data.frame(mca_result$eig)
  Vaccounted$Dimension <- factor(rownames(Vaccounted), levels = rownames(Vaccounted))

  p_scree <- ggplot(Vaccounted, aes(x = Dimension)) +
    geom_line(aes(y = `cumulative percentage of variance`, group = 1)) +
    geom_point(aes(y = `cumulative percentage of variance`)) +
    geom_col(aes(y = `percentage of variance`)) +
    geom_hline(aes(yintercept = minThresh), linetype = "dashed") +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(angle = 90))

  # Perform MCA with selected number of components
  mca_result <- MCA(DataSubset %>% select(all_of(VarsToReduce)),
                    ncp = nc, graph = FALSE)

  # Extract contributions and individual participant values from the MCA result
  MCAvar <- mca_result$var
  MCAind <- mca_result$ind
  var_contributions <- MCAvar$contrib
  ParticipantsValues <- MCAind$coord

  # Create a loading table and sort by contribution
  LoadingTable <- mca_result$var$coord %>% as.data.frame()
  LoadingTable <- psych::fa.sort(LoadingTable)

  # Reshape loading table for plotting
  mLoading <- LoadingTable %>% rownames_to_column() %>%
    select(-order) %>%
    pivot_longer(cols = colnames(LoadingTable %>% select(-order)),
                 names_to = "Dimension", values_to = "Contribution")

  # Highlight top contributors based on contribution threshold (>=1)
  mLoading <- mLoading %>% group_by(rowname) %>%
    arrange(Contribution) %>%
    mutate(TopContributor = abs(Contribution) >= 1)

  # Assign colors based on variable categories if provided
  if (!is.null(VariableCategories)) {
    mLoading$color <- VariableCategories[match(mLoading$rowname, VarsToReduce)]
  } else {
    mLoading$color <- "black"  # Default to black if no categories provided
  }

  # Adjust factor levels for proper display in the plot
  mLoading$rowname <- factor(mLoading$rowname, levels = rev(rownames(LoadingTable)))
  mLoading$Dimension <- factor(mLoading$Dimension, levels = colnames(LoadingTable) %>% removeString("order"))

  # Create a lollipop plot for loadings
  p <- ggplot(mLoading, aes(x = rowname, y = Contribution,
                            alpha = TopContributor, color = color)) +
    geom_segment(aes(xend = rowname, y = 0, yend = Contribution)) +
    geom_point() +
    facet_wrap(~Dimension, nrow = 1) +
    coord_flip() +
    geom_hline(yintercept = 0) +
    theme_bw()

  # Combine original data with participant scores from MCA
  CombinedData <- cbind(Data, MCAind)

  # Return results as a list
  return(list(p_scree = p_scree,
              pcaresults = mca_result,
              LoadingTable = LoadingTable,
              Scores = MCAind,
              CombinedData = CombinedData,
              Lollipop = p))
}
