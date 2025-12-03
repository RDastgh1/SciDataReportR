#' Create PCA table and visualization
#'
#' Perform principal component analysis (PCA) on specified variables and create visualizations.
#'
#' @param Data The dataset containing the variables for PCA.
#' @param VarsToReduce A character vector specifying the variables to include in the PCA.
#' @param VariableCategories Optional categorical vector, used to color the lollipop plot.
#'   Must be the same length and order as VarsToReduce if provided.
#' @param minThresh The minimum threshold for cumulative proportion of variance (default 0.85).
#' @param scale Logical, indicating whether to scale the data (divide by SD). Default TRUE.
#' @param center Logical, indicating whether to center the data (subtract mean). Default TRUE.
#' @param Ordinal Logical, indicating whether ordinal variables should be handled.
#'   Currently not used inside this function.
#' @param numComponents Number of principal components to compute. If NULL, chosen by minThresh.
#'
#' @return A list containing PCA results and visualizations:
#'   \item{p_scree}{Scree and variance plot (ggplot object).}
#'   \item{pcaresults}{psych::principal result (after psych::fa.sort).}
#'   \item{LoadingTable}{Data frame of loadings with variable names and labels.}
#'   \item{Scores}{Matrix or data frame of component scores.}
#'   \item{CombinedData}{Original Data with component scores appended.}
#'   \item{Lollipop}{Lollipop style loading plot (ggplot object).}
#'   \item{ScaleParams}{List with elements means, sds, center, scale.}
#'   \item{VarsUsed}{Character vector of variables actually used in PCA.}
#'   \item{Center}{Logical flag indicating if centering was used.}
#'   \item{Scale}{Logical flag indicating if scaling was used.}
#' @export
CreatePCATable <- function(Data,
                           VarsToReduce,
                           VariableCategories = NULL,
                           minThresh = 0.85,
                           scale = TRUE,
                           center = TRUE,
                           Ordinal = FALSE,
                           numComponents = NULL) {

  classcolors <- c(
    paletteer::paletteer_d("calecopal::superbloom2"),
    paletteer::paletteer_d("calecopal::vermillion"),
    paletteer::paletteer_d("fishualize::Antennarius_commerson"),
    paletteer::paletteer_d("fishualize::Bodianus_rufus")
  )

  # 1. Label codebook (always use labels if present)
  Data <- ReplaceMissingLabels(Data)
  lbls <- sjlabelled::get_label(Data)
  # fallback: where labels are NULL or "", use variable names
  lbls[is.null(lbls)] <- ""
  lbls[lbls == ""] <- colnames(Data)[lbls == ""]

  LabelCodebook <- data.frame(
    Variable = colnames(Data),
    Labels   = lbls,
    stringsAsFactors = FALSE
  )

  # 2. Subset to PCA variables
  DataSubset <- Data[VarsToReduce]

  is_num <- vapply(DataSubset, is.numeric, logical(1))
  if (!all(is_num)) {
    warning(
      "Dropping non numeric variables from PCA: ",
      paste(VarsToReduce[!is_num], collapse = ", ")
    )
    VarsToReduce <- VarsToReduce[is_num]
    DataSubset   <- DataSubset[VarsToReduce]
  }

  if (length(VarsToReduce) == 0) {
    stop("No numeric VarsToReduce available for PCA.")
  }

  # 3. Impute missing data with missRanger (same as your original)
  if (sum(is.na(DataSubset)) > 0) {
    set.seed(123456)
    colnames(DataSubset) <- make.names(VarsToReduce)
    DataSubset <- missRanger::missRanger(DataSubset, num.trees = 100)
    colnames(DataSubset) <- VarsToReduce
  }

  # 4. Scale and center once, and store parameters
  if (center || scale) {
    DataSubset_scaled <- scale(DataSubset, center = center, scale = scale)

    scale_means <- attr(DataSubset_scaled, "scaled:center")
    scale_sds   <- attr(DataSubset_scaled, "scaled:scale")

    if (is.null(scale_means)) {
      scale_means <- rep(0, length(VarsToReduce))
      names(scale_means) <- VarsToReduce
    }
    if (is.null(scale_sds)) {
      scale_sds <- rep(1, length(VarsToReduce))
      names(scale_sds) <- VarsToReduce
    }
  } else {
    DataSubset_scaled <- as.matrix(DataSubset)
    scale_means <- rep(0, length(VarsToReduce))
    scale_sds   <- rep(1, length(VarsToReduce))
    names(scale_means) <- VarsToReduce
    names(scale_sds)   <- VarsToReduce
  }

  n <- length(VarsToReduce)

  # 5. Scree PCA on scaled data
  fit1 <- psych::principal(
    r        = DataSubset_scaled,
    nfactors = n,
    rotate   = "none",
    scores   = FALSE
  )

  if (is.null(numComponents)) {
    nc <- min(which(fit1$Vaccounted[3, ] > minThresh))
  } else {
    nc <- numComponents
  }

  Vaccounted <- as.data.frame(t(fit1$Vaccounted))
  Vaccounted$Component <- factor(
    rownames(Vaccounted),
    levels = rownames(Vaccounted)
  )

  p_scree <- ggplot(Vaccounted, aes(x = Component)) +
    geom_line(aes(y = `Cumulative Var`, group = 1)) +
    geom_point(aes(y = `Cumulative Proportion`)) +
    geom_col(aes(y = `Proportion Var`)) +
    geom_hline(aes(yintercept = minThresh), linetype = "dashed") +
    theme_bw()

  # 6. Final PCA (varimax) on same scaled data
  fit <- psych::principal(
    r        = DataSubset_scaled,
    nfactors = nc,
    rotate   = "varimax",
    scores   = TRUE
  )
  fit <- psych::fa.sort(fit)

  # 7. Loading table
  LoadingTable <- as.data.frame(xtable::xtable(unclass(fit$loadings)))
  LoadingTable$Variable <- rownames(LoadingTable)

  LoadingTable <- dplyr::left_join(
    LoadingTable,
    LabelCodebook,
    by = "Variable",
    multiple = "first"
  )

  # 8. Scores and combined data
  scores <- fit$scores
  CombinedData <- cbind(Data, scores)

  # 9. Lollipop plot data
  rc_cols <- names(LoadingTable)[grepl("^RC", names(LoadingTable))]

  meltedLoading <- LoadingTable %>%
    tidyr::pivot_longer(cols = dplyr::all_of(rc_cols))

  meltedLoading$name <- factor(
    meltedLoading$name,
    levels = rc_cols
  )

  meltedLoading <- meltedLoading %>%
    dplyr::group_by(Variable) %>%
    dplyr::arrange(value, .by_group = TRUE) %>%
    dplyr::mutate(color = abs(value) > 0.4)

  meltedLoading$TopContributor <- meltedLoading$color

  if (is.null(VariableCategories)) {
    meltedLoading$color <- "black"
  } else {
    meltedLoading$color <- VariableCategories[
      match(meltedLoading$Variable, VarsToReduce)
    ]
  }

  meltedLoading$Labels <- factor(
    meltedLoading$Labels,
    levels = rev(LoadingTable$Labels)
  )

  if (is.null(VariableCategories)) {
    p <- ggplot(
      meltedLoading,
      aes(
        x     = Labels,
        y     = value,
        alpha = TopContributor,
        color = color
      )
    ) +
      coord_flip() +
      geom_segment(
        size = 2,
        aes(x = Labels, xend = Labels, y = 0, yend = value)
      ) +
      facet_wrap(vars(name), nrow = 1) +
      theme(
        legend.position = "none",
        strip.text.y    = element_text(angle = 45),
        axis.title.x    = element_blank(),
        axis.title.y    = element_blank()
      ) +
      geom_point() +
      scale_color_manual(values = c("black"))
  } else {
    p <- ggplot(
      meltedLoading,
      aes(
        x     = Labels,
        y     = value,
        alpha = TopContributor,
        color = color
      )
    ) +
      coord_flip() +
      geom_segment(
        size = 2,
        aes(x = Labels, xend = Labels, y = 0, yend = value)
      ) +
      facet_wrap(vars(name), nrow = 1) +
      theme(
        legend.position = "none",
        strip.text.y    = element_text(angle = 45),
        axis.title.x    = element_blank(),
        axis.title.y    = element_blank()
      ) +
      geom_point() +
      scale_color_manual(values = classcolors)
  }

  ScaleParams <- list(
    means  = scale_means,
    sds    = scale_sds,
    center = center,
    scale  = scale
  )

  return(list(
    p_scree      = p_scree,
    pcaresults   = fit,
    LoadingTable = LoadingTable,
    Scores       = scores,
    CombinedData = CombinedData,
    Lollipop     = p,
    ScaleParams  = ScaleParams,
    VarsUsed     = VarsToReduce,
    Center       = center,
    Scale        = scale
  ))
}
