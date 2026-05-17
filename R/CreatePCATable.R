#' Create PCA table and visualization
#'
#' Perform principal component analysis (PCA) on specified variables and create
#' visualizations. The default `"classic"` mode preserves the original
#' `psych::principal()` workflow for compatibility with existing
#' SciDataReportR analyses. The optional `"omics"` mode is designed for
#' high-dimensional data such as proteomics, metabolomics, and other settings
#' where the number of variables is much larger than the number of participants.
#'
#' In `"omics"` mode, the function avoids computing one component per input
#' variable. Instead, it computes a capped number of components using a faster
#' PCA backend, applies optional variance filtering, and then applies rotation
#' and sorting so loadings remain interpretable. This is useful when PCA is being
#' used for exploratory structure, visualization, latent biological signal, or
#' clustering rather than exhaustive factor decomposition.
#'
#' @param Data The dataset containing the variables for PCA.
#' @param VarsToReduce A character vector specifying the variables to include in
#'   the PCA.
#' @param VariableCategories Optional categorical vector, used to color the
#'   lollipop plot. Must be the same length and order as `VarsToReduce` if
#'   provided.
#' @param minThresh The minimum threshold for cumulative proportion of variance.
#'   Default is `0.85`. In `"omics"` mode, this threshold is evaluated only
#'   across the capped scree components.
#' @param scale Logical, indicating whether to scale the data by standard
#'   deviation. Default is `TRUE`.
#' @param center Logical, indicating whether to center the data by subtracting
#'   the mean. Default is `TRUE`.
#' @param Ordinal Logical, indicating whether ordinal variables should be
#'   handled. Currently not used inside this function.
#' @param numComponents Number of principal components to compute. If `NULL`,
#'   chosen by `minThresh` in `"classic"` mode and by capped scree results in
#'   `"omics"` mode.
#' @param Mode Character. Either `"classic"` or `"omics"`. `"classic"` preserves
#'   the original `psych::principal()` behavior. `"omics"` uses capped PCA logic
#'   for high-dimensional data.
#' @param backend Character. PCA backend to use. Options are `"psych"`,
#'   `"prcomp"`, and `"irlba"`. Default is `"psych"` to preserve backward
#'   compatibility. In `"omics"` mode, `"irlba"` is recommended for speed.
#' @param rotate Character. Rotation method. Options are `"varimax"` and
#'   `"none"`. Default is `"varimax"` to preserve interpretability and existing
#'   loading-table behavior.
#' @param maxComponents Maximum number of final components to retain in
#'   `"omics"` mode when `numComponents` is `NULL`. Default is `20`.
#' @param maxScreeComponents Maximum number of components used to estimate and
#'   plot the scree curve in `"omics"` mode. Default is `20`.
#' @param VarianceFilter Optional numeric value for variance filtering before
#'   PCA in `"omics"` mode. If `VarianceFilterMethod = "top_n"`, this keeps the
#'   top `VarianceFilter` most variable variables. If
#'   `VarianceFilterMethod = "variance_quantile"`, this keeps variables with
#'   variance at or above the specified quantile.
#' @param VarianceFilterMethod Character. Either `"top_n"` or
#'   `"variance_quantile"`.
#' @param MissingnessWarningThreshold Numeric threshold for warning about high
#'   variable-level missingness. Default is `0.20`.
#' @param ParticipantMissingnessWarningThreshold Numeric threshold for warning
#'   about high participant-level missingness. Default is `0.20`.
#' @param imputeMethod Character. Missing-data imputation method. Options are
#'   `"missRanger"` and `"median"`. Default is `"missRanger"` to preserve the
#'   original behavior.
#' @param SuppressWarnings Logical. If `TRUE`, suppresses PCA-specific warning
#'   messages from this function. Default is `FALSE`.
#'
#' @return A list containing PCA results and visualizations:
#'   \item{p_scree}{Scree and variance plot as a ggplot object.}
#'   \item{pcaresults}{PCA result object. In `"classic"` mode this is the
#'   `psych::principal()` result after `psych::fa.sort()`. In `"omics"` mode
#'   this is a harmonized list with rotated and sorted loadings, scores, and
#'   variance information.}
#'   \item{LoadingTable}{Data frame of rotated loadings with variable names and
#'   labels.}
#'   \item{Scores}{Matrix or data frame of component scores.}
#'   \item{CombinedData}{Original data with component scores appended.}
#'   \item{Lollipop}{Lollipop-style loading plot as a ggplot object.}
#'   \item{ScaleParams}{List with means, standard deviations, center, and scale.}
#'   \item{VarsUsed}{Character vector of variables actually used in PCA.}
#'   \item{VarianceTable}{Table of variance explained by component.}
#'   \item{Mode}{PCA mode used.}
#'   \item{Backend}{PCA backend used.}
#'   \item{Center}{Logical flag indicating if centering was used.}
#'   \item{Scale}{Logical flag indicating if scaling was used.}
#'
#' @examples
#' PCA <- CreatePCATable(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   numComponents = 3
#' )
#'
#' PCA_Omics <- CreatePCATable(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   Mode = "omics",
#'   backend = "prcomp",
#'   maxComponents = 5,
#'   maxScreeComponents = 5,
#'   imputeMethod = "median"
#' )
#'
#' @export
CreatePCATable <- function(Data,
                           VarsToReduce,
                           VariableCategories = NULL,
                           minThresh = 0.85,
                           scale = TRUE,
                           center = TRUE,
                           Ordinal = FALSE,
                           numComponents = NULL,
                           Mode = c("classic", "omics"),
                           backend = c("psych", "prcomp", "irlba"),
                           rotate = c("varimax", "none"),
                           maxComponents = 20,
                           maxScreeComponents = 20,
                           VarianceFilter = NULL,
                           VarianceFilterMethod = c("top_n", "variance_quantile"),
                           MissingnessWarningThreshold = 0.20,
                           ParticipantMissingnessWarningThreshold = 0.20,
                           imputeMethod = c("missRanger", "median"),
                           SuppressWarnings = FALSE) {

  Mode <- match.arg(Mode)
  backend <- match.arg(backend)
  rotate <- match.arg(rotate)
  VarianceFilterMethod <- match.arg(VarianceFilterMethod)
  imputeMethod <- match.arg(imputeMethod)

  if (Mode == "omics" && backend == "psych" && !SuppressWarnings) {
    warning(
      "Mode = 'omics' works best with backend = 'irlba' or backend = 'prcomp'. ",
      "backend = 'psych' is allowed for compatibility, but may be slow with high-dimensional omics data."
    )
  }

  classcolors <- c(
    paletteer::paletteer_d("calecopal::superbloom2"),
    paletteer::paletteer_d("calecopal::vermillion"),
    paletteer::paletteer_d("fishualize::Antennarius_commerson"),
    paletteer::paletteer_d("fishualize::Bodianus_rufus")
  )

  # 1. Label codebook (always use labels if present)
  Data <- ReplaceMissingLabels(Data)
  lbls <- sjlabelled::get_label(Data)

  lbls[is.null(lbls)] <- ""
  lbls[lbls == ""] <- colnames(Data)[lbls == ""]

  LabelCodebook <- data.frame(
    Variable = colnames(Data),
    Labels   = lbls,
    stringsAsFactors = FALSE
  )

  # 2. Subset to PCA variables
  missing_vars <- setdiff(VarsToReduce, colnames(Data))

  if (length(missing_vars) > 0) {
    stop(
      "Some VarsToReduce are missing from Data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

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

  if (Mode == "classic" &&
      length(VarsToReduce) > nrow(DataSubset) &&
      !SuppressWarnings) {
    warning(
      "Detected high-dimensional PCA setting: ",
      length(VarsToReduce), " variables and ", nrow(DataSubset), " rows. ",
      "Classic mode preserves the original psych::principal() behavior, but ",
      "Mode = 'omics' may be much faster and more stable."
    )
  }

  # 3. Warn about missingness
  variable_missingness <- colMeans(is.na(DataSubset))
  participant_missingness <- rowMeans(is.na(DataSubset))

  high_missing_vars <- names(variable_missingness)[
    variable_missingness > MissingnessWarningThreshold
  ]

  high_missing_participants <- which(
    participant_missingness > ParticipantMissingnessWarningThreshold
  )

  if (length(high_missing_vars) > 0 && !SuppressWarnings) {
    warning(
      length(high_missing_vars),
      " PCA variables have missingness above ",
      MissingnessWarningThreshold,
      ". Consider filtering before PCA."
    )
  }

  if (length(high_missing_participants) > 0 && !SuppressWarnings) {
    warning(
      length(high_missing_participants),
      " rows have missingness above ",
      ParticipantMissingnessWarningThreshold,
      ". Consider filtering participants before PCA."
    )
  }

  # 4. Optional variance filtering for omics-scale PCA
  VarianceTable <- NULL

  if (Mode == "omics" && !is.null(VarianceFilter)) {
    variable_variance <- vapply(
      DataSubset,
      stats::var,
      numeric(1),
      na.rm = TRUE
    )

    VarianceTable <- data.frame(
      Variable = names(variable_variance),
      Variance = as.numeric(variable_variance),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(dplyr::desc(Variance))

    if (VarianceFilterMethod == "top_n") {
      VarsKeep <- VarianceTable %>%
        dplyr::slice_head(n = VarianceFilter) %>%
        dplyr::pull(Variable)
    }

    if (VarianceFilterMethod == "variance_quantile") {
      variance_cutoff <- stats::quantile(
        variable_variance,
        probs = VarianceFilter,
        na.rm = TRUE
      )

      VarsKeep <- VarianceTable %>%
        dplyr::filter(Variance >= variance_cutoff) %>%
        dplyr::pull(Variable)
    }

    DataSubset <- DataSubset[VarsKeep]
    VarsToReduce <- VarsKeep

    if (!is.null(VariableCategories)) {
      VariableCategories <- VariableCategories[
        match(VarsToReduce, names(Data[VarsToReduce]))
      ]
    }
  }

  # 5. Impute missing data
  if (sum(is.na(DataSubset)) > 0) {
    if (imputeMethod == "missRanger") {
      set.seed(123456)
      colnames(DataSubset) <- make.names(VarsToReduce)
      DataSubset <- missRanger::missRanger(DataSubset, num.trees = 100)
      colnames(DataSubset) <- VarsToReduce
    }

    if (imputeMethod == "median") {
      DataSubset <- DataSubset %>%
        dplyr::mutate(
          dplyr::across(
            dplyr::everything(),
            ~ ifelse(is.na(.x), stats::median(.x, na.rm = TRUE), .x)
          )
        )
    }
  }

  # 6. Scale and center once, and store parameters
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

  zero_sd_vars <- names(scale_sds)[scale_sds == 0 | is.na(scale_sds)]

  if (length(zero_sd_vars) > 0) {
    stop(
      "Some PCA variables have zero or undefined standard deviation after preprocessing: ",
      paste(zero_sd_vars, collapse = ", ")
    )
  }

  n_vars <- length(VarsToReduce)

  # 7. Classic PCA path
  if (Mode == "classic") {
    fit1 <- psych::principal(
      r        = DataSubset_scaled,
      nfactors = n_vars,
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

    fit <- psych::principal(
      r        = DataSubset_scaled,
      nfactors = nc,
      rotate   = rotate,
      scores   = TRUE
    )

    fit <- psych::fa.sort(fit)

    LoadingTable <- as.data.frame(xtable::xtable(unclass(fit$loadings)))
    LoadingTable$Variable <- rownames(LoadingTable)

    scores <- fit$scores
  }

  # 8. Omics PCA path
  if (Mode == "omics") {
    max_possible_components <- min(
      nrow(DataSubset_scaled) - 1,
      ncol(DataSubset_scaled)
    )

    scree_n <- min(
      maxScreeComponents,
      max_possible_components
    )

    if (scree_n < 1) {
      stop("Not enough rows or variables to compute PCA.")
    }

    if (backend == "irlba") {
      if (!requireNamespace("irlba", quietly = TRUE)) {
        stop(
          "The irlba package is required when backend = 'irlba'. ",
          "Install it with install.packages('irlba') or use backend = 'prcomp'."
        )
      }

      fit_scree <- irlba::prcomp_irlba(
        DataSubset_scaled,
        n = scree_n,
        center = FALSE,
        scale. = FALSE
      )
    }

    if (backend == "prcomp") {
      fit_scree <- stats::prcomp(
        DataSubset_scaled,
        center = FALSE,
        scale. = FALSE,
        rank. = scree_n
      )
    }

    if (backend == "psych") {
      fit_scree <- psych::principal(
        r        = DataSubset_scaled,
        nfactors = scree_n,
        rotate   = "none",
        scores   = FALSE
      )
    }

    if (backend %in% c("irlba", "prcomp")) {
      eigen_values <- fit_scree$sdev^2
      proportion_var <- eigen_values / sum(eigen_values)
      cumulative_var <- cumsum(proportion_var)

      Vaccounted <- data.frame(
        `SS loadings` = eigen_values,
        `Proportion Var` = proportion_var,
        `Cumulative Var` = cumulative_var,
        `Proportion Explained` = proportion_var,
        `Cumulative Proportion` = cumulative_var,
        check.names = FALSE
      )

      rownames(Vaccounted) <- paste0("PC", seq_len(nrow(Vaccounted)))
    }

    if (backend == "psych") {
      Vaccounted <- as.data.frame(t(fit_scree$Vaccounted))
    }

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

    if (is.null(numComponents)) {
      nc_from_thresh <- min(which(Vaccounted$`Cumulative Proportion` > minThresh))

      if (is.infinite(nc_from_thresh) || is.na(nc_from_thresh)) {
        nc_from_thresh <- min(maxComponents, scree_n)
      }

      nc <- min(nc_from_thresh, maxComponents, scree_n)
    } else {
      nc <- min(numComponents, scree_n)
    }

    if (backend == "irlba") {
      fit_raw <- irlba::prcomp_irlba(
        DataSubset_scaled,
        n = nc,
        center = FALSE,
        scale. = FALSE
      )
    }

    if (backend == "prcomp") {
      fit_raw <- stats::prcomp(
        DataSubset_scaled,
        center = FALSE,
        scale. = FALSE,
        rank. = nc
      )
    }

    if (backend == "psych") {
      fit <- psych::principal(
        r        = DataSubset_scaled,
        nfactors = nc,
        rotate   = rotate,
        scores   = TRUE
      )

      fit <- psych::fa.sort(fit)

      LoadingTable <- as.data.frame(xtable::xtable(unclass(fit$loadings)))
      LoadingTable$Variable <- rownames(LoadingTable)

      scores <- fit$scores
    }

    if (backend %in% c("irlba", "prcomp")) {
      loadings_raw <- sweep(
        fit_raw$rotation,
        2,
        fit_raw$sdev,
        FUN = "*"
      )

      colnames(loadings_raw) <- paste0("RC", seq_len(ncol(loadings_raw)))
      rownames(loadings_raw) <- VarsToReduce

      scores_raw <- fit_raw$x
      colnames(scores_raw) <- paste0("RC", seq_len(ncol(scores_raw)))

      if (rotate == "varimax") {
        rotated <- stats::varimax(loadings_raw)

        loadings_final <- as.matrix(rotated$loadings)
        scores_final <- scores_raw %*% rotated$rotmat
      } else {
        loadings_final <- loadings_raw
        scores_final <- scores_raw
      }

      # Keep fa.sort-style interpretability while preserving score/loading matches
      sort_index <- apply(abs(loadings_final), 2, max)
      component_order <- order(sort_index, decreasing = TRUE)

      loadings_final <- loadings_final[, component_order, drop = FALSE]
      scores_final <- scores_final[, component_order, drop = FALSE]

      colnames(loadings_final) <- paste0("RC", seq_len(ncol(loadings_final)))
      colnames(scores_final) <- paste0("RC", seq_len(ncol(scores_final)))

      row_order <- apply(abs(loadings_final), 1, which.max)
      row_max <- apply(abs(loadings_final), 1, max)
      loading_row_order <- order(row_order, -row_max)

      loadings_final <- loadings_final[loading_row_order, , drop = FALSE]

      class(loadings_final) <- "loadings"

      LoadingTable <- as.data.frame(xtable::xtable(unclass(loadings_final)))
      LoadingTable$Variable <- rownames(LoadingTable)

      scores <- scores_final

      rotated_ss <- colSums(unclass(loadings_final)^2)
      rotated_prop <- rotated_ss / sum(rotated_ss)

      fit <- list(
        values = fit_raw$sdev^2,
        loadings = loadings_final,
        scores = scores,
        Vaccounted = rbind(
          `SS loadings` = rotated_ss,
          `Proportion Var` = rotated_prop,
          `Cumulative Var` = cumsum(rotated_prop),
          `Proportion Explained` = rotated_prop,
          `Cumulative Proportion` = cumsum(rotated_prop)
        ),
        rotation = rotate,
        backend = backend,
        mode = Mode,
        call = match.call()
      )

      class(fit) <- c("SciDataReportR_pca", "list")
    }
  }

  # 9. Loading table
  LoadingTable <- dplyr::left_join(
    LoadingTable,
    LabelCodebook,
    by = "Variable",
    multiple = "first"
  )

  # 10. Scores and combined data
  Scores <- as.data.frame(scores)
  CombinedData <- cbind(Data, Scores)

  # 11. Lollipop plot data
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
        linewidth = 2,
        aes(x = Labels, xend = Labels, y = 0, yend = value)
      )+
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
    Scores       = Scores,
    CombinedData = CombinedData,
    Lollipop     = p,
    ScaleParams  = ScaleParams,
    VarsUsed     = VarsToReduce,
    VarianceTable = VarianceTable,
    Mode         = Mode,
    Backend      = backend,
    Center       = center,
    Scale        = scale
  ))
}
