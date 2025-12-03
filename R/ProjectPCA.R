#' Project PCA scores onto new data
#'
#' Use an existing PCA solution (either a PCA object from CreatePCATable or
#' a loading table) to compute principal component scores on a new dataset.
#'
#' @param Data Data frame on which to project PCA scores.
#' @param VarsToReduce Optional character vector of variable names to use.
#'   If NULL, uses all variables that appear in both Data and the PCA solution.
#' @param PCAInput Either:
#'   - the full object returned by CreatePCATable (when InputType = "PCAObj"), or
#'   - a loading table like CreatePCATable()$LoadingTable (when InputType = "LoadingTable").
#' @param InputType One of "PCAObj" or "LoadingTable".
#' @param center Logical; only used when InputType is "LoadingTable".
#'   For "PCAObj", the centering choice is taken from PCAInput$ScaleParams$center
#'   and this argument is ignored (with a warning if it conflicts).
#' @param scale Logical; only used when InputType is "LoadingTable".
#'   For "PCAObj", the scaling choice is taken from PCAInput$ScaleParams$scale
#'   and this argument is ignored (with a warning if it conflicts).
#'
#' @return A list with:
#'   \item{Scores}{Data frame of projected PCA scores.}
#'   \item{CombinedData}{Original Data with scores appended as new columns.}
#'   \item{LoadingsUsed}{Matrix of loadings used for projection.}
#'   \item{PCAObj}{PCA object used (if InputType is "PCAObj"), otherwise NULL.}
#'   \item{VarsUsed}{Variables used from Data for projection.}
#'   \item{Center}{Logical flag indicating whether centering was applied for projection.}
#'   \item{Scale}{Logical flag indicating whether scaling was applied for projection.}
#' @export
ProjectPCA <- function(Data,
                       VarsToReduce = NULL,
                       PCAInput,
                       InputType = c("PCAObj", "LoadingTable"),
                       center = TRUE,
                       scale  = TRUE) {

  InputType <- match.arg(InputType)

  # 1. Resolve loadings/weights, variable names, and scale parameters ---------
  if (InputType == "PCAObj") {

    if (!is.list(PCAInput) || is.null(PCAInput$pcaresults)) {
      stop("For InputType = 'PCAObj', PCAInput must be the output of CreatePCATable().")
    }

    fit <- PCAInput$pcaresults

    if (is.null(fit$loadings)) {
      stop("PCA object does not contain loadings.")
    }

    loading_mat  <- as.matrix(unclass(fit$loadings))
    all_pca_vars <- rownames(loading_mat)

    # Use weights for scoring (matches psych::principal)
    if (is.null(fit$weights)) {
      stop("PCA object does not contain weights; cannot project consistently with psych::principal().")
    }
    weight_mat <- as.matrix(fit$weights)

    if (is.null(PCAInput$ScaleParams)) {
      stop("PCAInput is missing ScaleParams. Re run CreatePCATable with the updated version.")
    }

    scale_means   <- PCAInput$ScaleParams$means
    scale_sds     <- PCAInput$ScaleParams$sds
    train_center  <- isTRUE(PCAInput$ScaleParams$center)
    train_scale   <- isTRUE(PCAInput$ScaleParams$scale)

    # If user passes different center/scale, warn and override
    if (center != train_center || scale != train_scale) {
      warning(
        "center/scale arguments to ProjectPCA() are ignored for InputType = 'PCAObj'. ",
        "Using Center = ", train_center, ", Scale = ", train_scale, " from PCAInput$ScaleParams."
      )
    }

  } else {

    LT <- PCAInput
    if (!is.data.frame(LT)) {
      stop("For InputType = 'LoadingTable', PCAInput must be a data frame of loadings.")
    }

    num_cols <- vapply(LT, is.numeric, logical(1))
    if (!any(num_cols)) {
      stop("LoadingTable does not contain numeric loading columns.")
    }

    loading_mat <- as.matrix(LT[, num_cols, drop = FALSE])

    if (!"Variable" %in% names(LT)) {
      stop("LoadingTable must contain a 'Variable' column with variable names.")
    }
    rownames(loading_mat) <- LT$Variable
    all_pca_vars <- LT$Variable

    # No weights / no stored scaling; rely on arguments
    weight_mat   <- NULL
    scale_means  <- NULL
    scale_sds    <- NULL
    train_center <- center
    train_scale  <- scale
  }

  # 2. Determine which variables to use --------------------------------------

  if (is.null(VarsToReduce)) {
    VarsToReduce <- intersect(all_pca_vars, names(Data))
    if (length(VarsToReduce) == 0) {
      stop("No overlap between PCA variables and columns in Data.")
    }
  } else {
    missing_in_data <- setdiff(VarsToReduce, names(Data))
    if (length(missing_in_data) > 0) {
      stop(
        "The following VarsToReduce are not in Data: ",
        paste(missing_in_data, collapse = ", ")
      )
    }
    missing_in_pca <- setdiff(VarsToReduce, all_pca_vars)
    if (length(missing_in_pca) > 0) {
      stop(
        "The following VarsToReduce are not in the PCA solution: ",
        paste(missing_in_pca, collapse = ", ")
      )
    }
  }

  DataSubset <- Data[VarsToReduce]

  is_num <- vapply(DataSubset, is.numeric, logical(1))
  if (!all(is_num)) {
    warning(
      "Dropping non numeric variables from projection: ",
      paste(VarsToReduce[!is_num], collapse = ", ")
    )
    VarsToReduce <- VarsToReduce[is_num]
    DataSubset   <- DataSubset[VarsToReduce]
  }

  if (length(VarsToReduce) == 0) {
    stop("No numeric variables available for projection.")
  }

  if (!all(VarsToReduce %in% all_pca_vars)) {
    missing_in_loadings <- setdiff(VarsToReduce, all_pca_vars)
    stop(
      "The following VarsToReduce do not have loadings/weights: ",
      paste(missing_in_loadings, collapse = ", ")
    )
  }

  # Align loadings/weights rows to VarsToReduce
  if (InputType == "PCAObj") {
    loading_mat_use <- loading_mat[VarsToReduce, , drop = FALSE]
    weight_mat_use  <- weight_mat[VarsToReduce, , drop = FALSE]
  } else {
    loading_mat_use <- loading_mat[VarsToReduce, , drop = FALSE]
    weight_mat_use  <- NULL
  }

  # 3. Prepare X matrix and apply training scaling if available --------------

  X <- as.matrix(DataSubset)

  if (InputType == "PCAObj") {

    if (train_center || train_scale) {
      if (!all(VarsToReduce %in% names(scale_means)) ||
          !all(VarsToReduce %in% names(scale_sds))) {
        stop("ScaleParams do not cover all VarsToReduce. Check PCAInput$ScaleParams.")
      }
      mu  <- scale_means[VarsToReduce]
      sig <- scale_sds[VarsToReduce]

      if (train_center) {
        X <- sweep(X, 2, mu, FUN = "-")
      }
      if (train_scale) {
        sig[sig == 0 | is.na(sig)] <- 1
        X <- sweep(X, 2, sig, FUN = "/")
      }
    }

  } else {

    if (train_center || train_scale) {
      X <- scale(X, center = train_center, scale = train_scale)
    }
  }

  # 4. Compute scores --------------------------------------------------------

  if (InputType == "PCAObj") {
    # Use regression weights to match psych::principal scores
    Scores <- X %*% weight_mat_use
    colnames(Scores) <- colnames(weight_mat_use)
  } else {
    # Use loadings when only a loading table is available
    Scores <- X %*% loading_mat_use
    if (is.null(colnames(Scores))) {
      colnames(Scores) <- paste0("PC", seq_len(ncol(Scores)))
    }
  }

  Scores <- as.data.frame(Scores)
  CombinedData <- cbind(Data, Scores)

  out <- list(
    Scores       = Scores,
    CombinedData = CombinedData,
    LoadingsUsed = loading_mat_use,
    PCAObj       = if (InputType == "PCAObj") fit else NULL,
    VarsUsed     = VarsToReduce,
    Center       = train_center,
    Scale        = train_scale
  )
  class(out) <- c("ProjectPCAObj", class(out))
  out
}
