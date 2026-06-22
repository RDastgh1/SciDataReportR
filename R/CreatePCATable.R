#' Create a reusable PCA object and visualizations
#'
#' Perform principal component analysis (PCA) on specified variables and return
#' reusable PCA results, scores, loading tables, combined data, and plots.
#'
#' The default `"classic"` mode preserves the original `psych::principal()`
#' workflow for compatibility with existing SciDataReportR analyses. The
#' optional `"omics"` mode is designed for high-dimensional data such as
#' proteomics, metabolomics, flow cytometry, FACS, transcriptomics, and other
#' settings where the number of variables may be large relative to the number
#' of participants.
#'
#' Missing data are not imputed by default. With the default
#' `MissingDataStrategy = "complete_cases"`, PCA is fit using only rows that
#' are complete across the final PCA variables, and rows with missing PCA inputs
#' receive `NA` component scores in `Scores` and `CombinedData`. To impute
#' missing values, set `MissingDataStrategy = "impute"` explicitly and choose
#' an imputation method with `ImputeMethod`.
#'
#' The function uses raw variable names internally, but by default uses variable
#' labels in human-facing outputs when labels are available. Variables with
#' zero or undefined standard deviation after preprocessing are dropped
#' automatically with a warning so PCA can continue without manual preprocessing.
#'
#' @param Data A data frame containing the variables for PCA.
#' @param VarsToReduce Character vector of raw variable names to include in PCA.
#' @param VariableCategories Optional vector used to color variables in the
#'   lollipop loading plot. If supplied, it should be the same length and order
#'   as `VarsToReduce`.
#' @param Relabel Logical. If `TRUE`, variable labels are used in output tables
#'   and plots when available. If `FALSE`, raw variable names are used. Default
#'   is `TRUE`.
#' @param minThresh Numeric threshold for cumulative variance used to select the
#'   number of components when `numComponents = NULL`. Default is `0.85`.
#' @param scale Logical. If `TRUE`, variables are scaled by their standard
#'   deviation before PCA. Default is `TRUE`.
#' @param center Logical. If `TRUE`, variables are centered before PCA. Default
#'   is `TRUE`.
#' @param Ordinal Logical. Placeholder retained for backward compatibility.
#'   Currently not used inside this function.
#' @param numComponents Optional integer number of components to retain. If
#'   `NULL`, the number is selected using `minThresh`.
#' @param Mode Character. Either `"classic"` or `"omics"`. `"classic"` preserves
#'   the original `psych::principal()` behavior. `"omics"` uses capped PCA logic
#'   for higher-dimensional data.
#' @param backend Character. PCA backend to use. Options are `"psych"`,
#'   `"prcomp"`, and `"irlba"`. Default is `"psych"` for compatibility.
#' @param rotate Character. Rotation method. Options are `"varimax"` and
#'   `"none"`. Default is `"varimax"`.
#' @param maxComponents Integer. Maximum number of final components to retain in
#'   `"omics"` mode when `numComponents = NULL`. Default is `20`.
#' @param maxScreeComponents Integer. Maximum number of components used to
#'   estimate and plot the scree curve in `"omics"` mode. Default is `20`.
#' @param VarianceFilter Optional numeric value for variance filtering before
#'   PCA in `"omics"` mode. If `VarianceFilterMethod = "top_n"`, keeps the top
#'   `VarianceFilter` most variable variables. If
#'   `VarianceFilterMethod = "variance_quantile"`, keeps variables with
#'   variance at or above the specified quantile.
#' @param VarianceFilterMethod Character. Either `"top_n"` or
#'   `"variance_quantile"`. Default is `"top_n"`.
#' @param MissingnessWarningThreshold Numeric threshold for warning about
#'   variable-level missingness. Default is `0.20`.
#' @param ParticipantMissingnessWarningThreshold Numeric threshold for warning
#'   about participant-level missingness. Default is `0.20`.
#' @param MissingDataStrategy Character. Missing-data handling strategy. Options
#'   are `"complete_cases"`, `"impute"`, and `"stop"`. Default is
#'   `"complete_cases"`, which fits PCA only on rows complete across the final
#'   PCA variables and returns `NA` PCA scores for incomplete rows.
#' @param ImputeMethod Character. Missing-data imputation method used only when
#'   `MissingDataStrategy = "impute"`. Options are `"missRanger"` and
#'   `"median"`. Default is `"missRanger"`.
#' @param MaxMissingForImputation Numeric value between 0 and 1. When
#'   `MissingDataStrategy = "impute"`, rows with missingness greater than this
#'   value across the final PCA variables are excluded from PCA scoring and
#'   receive `NA` component scores. Default is `0.20`, meaning rows can be
#'   imputed if at least 80 percent of PCA variables are observed.
#' @param ImputeRowsWithAllMissing Logical. If `FALSE`, rows with 100 percent
#'   missingness across PCA variables are not imputed even if
#'   `MaxMissingForImputation = 1`. Default is `FALSE`.
#' @param imputeMethod Deprecated compatibility alias for `ImputeMethod`. If
#'   supplied, it sets `ImputeMethod` and uses `MissingDataStrategy = "impute"`
#'   unless `MissingDataStrategy` was explicitly set.
#' @param SuppressWarnings Logical. If `TRUE`, suppresses PCA-specific warning
#'   messages from this function. Default is `FALSE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{p_scree}{A ggplot scree and cumulative variance plot.}
#'   \item{pcaresults}{The PCA result object. In `"classic"` mode this is a
#'   `psych::principal()` result after `psych::fa.sort()`. In `"omics"` mode
#'   with `"prcomp"` or `"irlba"`, this is a harmonized list containing loadings,
#'   scores, variance information, backend, mode, and rotation.}
#'   \item{LoadingTable}{A data frame of component loadings with raw variable
#'   names, labels, and duplicate-safe plot labels.}
#'   \item{Scores}{A data frame of component scores aligned to the original row
#'   order of `Data`. Rows not used for PCA scoring receive `NA` scores.}
#'   \item{CombinedData}{The original input data with component scores appended.}
#'   \item{Lollipop}{A ggplot lollipop loading plot.}
#'   \item{ScaleParams}{A list containing centering and scaling parameters.}
#'   \item{VarsUsed}{Character vector of variables actually used in PCA after
#'   preprocessing.}
#'   \item{VarianceTable}{A variance table used for optional omics variance
#'   filtering, or `NULL`.}
#'   \item{Preprocessing}{A list documenting variables dropped during
#'   preprocessing, missingness summaries, rows used for PCA, rows excluded from
#'   PCA, and rows imputed.}
#'   \item{Mode}{The PCA mode used.}
#'   \item{Backend}{The PCA backend used.}
#'   \item{Center}{Logical flag indicating whether centering was used.}
#'   \item{Scale}{Logical flag indicating whether scaling was used.}
#' }
#'
#' @examples
#' PCA <- CreatePCAObject(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   numComponents = 3
#' )
#'
#' PCA_imputed <- CreatePCAObject(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   numComponents = 3,
#'   MissingDataStrategy = "impute",
#'   ImputeMethod = "median",
#'   MaxMissingForImputation = 0.20
#' )
#'
#' PCA_raw_names <- CreatePCAObject(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   Relabel = FALSE,
#'   numComponents = 3
#' )
#'
#' PCA_omics <- CreatePCAObject(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   Mode = "omics",
#'   backend = "prcomp",
#'   maxComponents = 5,
#'   maxScreeComponents = 5
#' )
#'
#' @export
CreatePCAObject <- function(Data,
                            VarsToReduce,
                            VariableCategories = NULL,
                            Relabel = TRUE,
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
                            MissingDataStrategy = c("complete_cases", "impute", "stop"),
                            ImputeMethod = c("missRanger", "median"),
                            MaxMissingForImputation = 0.20,
                            ImputeRowsWithAllMissing = FALSE,
                            imputeMethod = NULL,
                            SuppressWarnings = FALSE) {

  missing_data_strategy_missing <- missing(MissingDataStrategy)

  Mode <- match.arg(Mode)
  backend <- match.arg(backend)
  rotate <- match.arg(rotate)
  VarianceFilterMethod <- match.arg(VarianceFilterMethod)
  MissingDataStrategy <- match.arg(MissingDataStrategy)

  if (!is.null(imputeMethod)) {
    if (!is.character(imputeMethod) || length(imputeMethod) != 1) {
      stop("imputeMethod must be NULL or a single character value.")
    }

    ImputeMethod <- imputeMethod

    if (missing_data_strategy_missing) {
      MissingDataStrategy <- "impute"
    }

    if (!SuppressWarnings) {
      warning(
        "imputeMethod is deprecated. Please use ImputeMethod instead. ",
        "Because imputeMethod was supplied, MissingDataStrategy was set to 'impute' ",
        "unless MissingDataStrategy was explicitly provided."
      )
    }
  }

  ImputeMethod <- match.arg(ImputeMethod)

  if (Mode == "omics" && backend == "psych" && !SuppressWarnings) {
    warning(
      "Mode = 'omics' works best with backend = 'irlba' or backend = 'prcomp'. ",
      "backend = 'psych' is allowed for compatibility, but may be slow with high-dimensional omics data."
    )
  }

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }

  if (!is.character(VarsToReduce) || length(VarsToReduce) == 0) {
    stop("VarsToReduce must be a non-empty character vector of variable names.")
  }

  if (!is.logical(Relabel) || length(Relabel) != 1) {
    stop("Relabel must be a single logical value: TRUE or FALSE.")
  }

  if (!is.logical(SuppressWarnings) || length(SuppressWarnings) != 1) {
    stop("SuppressWarnings must be a single logical value: TRUE or FALSE.")
  }

  if (!is.logical(ImputeRowsWithAllMissing) || length(ImputeRowsWithAllMissing) != 1) {
    stop("ImputeRowsWithAllMissing must be a single logical value: TRUE or FALSE.")
  }

  if (!is.numeric(MaxMissingForImputation) ||
      length(MaxMissingForImputation) != 1 ||
      is.na(MaxMissingForImputation) ||
      MaxMissingForImputation < 0 ||
      MaxMissingForImputation > 1) {
    stop("MaxMissingForImputation must be a single numeric value between 0 and 1.")
  }

  if (!is.null(numComponents)) {
    if (!is.numeric(numComponents) || length(numComponents) != 1 || numComponents < 1) {
      stop("numComponents must be NULL or a single positive number.")
    }

    numComponents <- as.integer(numComponents)
  }

  if (!is.numeric(minThresh) || length(minThresh) != 1 || minThresh <= 0 || minThresh > 1) {
    stop("minThresh must be a single numeric value greater than 0 and less than or equal to 1.")
  }

  classcolors <- c(
    paletteer::paletteer_d("calecopal::superbloom2"),
    paletteer::paletteer_d("calecopal::vermillion"),
    paletteer::paletteer_d("fishualize::Antennarius_commerson"),
    paletteer::paletteer_d("fishualize::Bodianus_rufus")
  )

  # Validate inputs

  missing_vars <- setdiff(VarsToReduce, colnames(Data))

  if (length(missing_vars) > 0) {
    stop(
      "Some VarsToReduce are missing from Data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (!is.null(VariableCategories) && length(VariableCategories) != length(VarsToReduce)) {
    stop(
      "VariableCategories must be NULL or the same length as VarsToReduce. ",
      "Received ", length(VariableCategories), " categories for ",
      length(VarsToReduce), " variables."
    )
  }

  # Prepare labels

  Data <- ReplaceMissingLabels(Data)

  lbls <- sjlabelled::get_label(Data)

  if (is.null(lbls)) {
    lbls <- rep("", ncol(Data))
    names(lbls) <- colnames(Data)
  }

  lbls[is.na(lbls)] <- ""
  lbls[lbls == ""] <- colnames(Data)[lbls == ""]

  LabelCodebook <- data.frame(
    Variable = colnames(Data),
    Labels = as.character(lbls),
    stringsAsFactors = FALSE
  )

  if (!Relabel) {
    LabelCodebook$Labels <- LabelCodebook$Variable
  }

  # Prepare data

  VarsOriginal <- VarsToReduce
  VariableCategoryTable <- NULL

  if (!is.null(VariableCategories)) {
    VariableCategoryTable <- data.frame(
      Variable = VarsToReduce,
      Category = as.character(VariableCategories),
      stringsAsFactors = FALSE
    )
  }

  DataSubset <- Data[VarsToReduce]

  is_num <- vapply(DataSubset, is.numeric, logical(1))
  dropped_non_numeric <- VarsToReduce[!is_num]

  if (length(dropped_non_numeric) > 0) {
    if (!SuppressWarnings) {
      warning(
        "Dropping non-numeric PCA variables: ",
        paste(dropped_non_numeric, collapse = ", ")
      )
    }

    VarsToReduce <- VarsToReduce[is_num]
    DataSubset <- DataSubset[VarsToReduce]
  }

  if (length(VarsToReduce) < 2) {
    stop("At least 2 numeric VarsToReduce are required for PCA.")
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
      ": ",
      paste(high_missing_vars, collapse = ", "),
      ". Consider filtering before PCA."
    )
  }

  if (length(high_missing_participants) > 0 && !SuppressWarnings) {
    warning(
      length(high_missing_participants),
      " rows have missingness above ",
      ParticipantMissingnessWarningThreshold,
      ". MissingDataStrategy controls whether these rows are excluded, imputed, or trigger an error."
    )
  }

  # Optional variance filtering

  VarianceTable <- NULL
  dropped_by_variance_filter <- character(0)

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
      if (!is.numeric(VarianceFilter) ||
          length(VarianceFilter) != 1 ||
          VarianceFilter < 1) {
        stop("When VarianceFilterMethod = 'top_n', VarianceFilter must be a single positive number.")
      }

      VarsKeep <- VarianceTable %>%
        dplyr::slice_head(n = as.integer(VarianceFilter)) %>%
        dplyr::pull(Variable)
    }

    if (VarianceFilterMethod == "variance_quantile") {
      if (!is.numeric(VarianceFilter) ||
          length(VarianceFilter) != 1 ||
          VarianceFilter < 0 ||
          VarianceFilter > 1) {
        stop(
          "When VarianceFilterMethod = 'variance_quantile', ",
          "VarianceFilter must be a single numeric value between 0 and 1."
        )
      }

      variance_cutoff <- stats::quantile(
        variable_variance,
        probs = VarianceFilter,
        na.rm = TRUE
      )

      VarsKeep <- VarianceTable %>%
        dplyr::filter(Variance >= variance_cutoff) %>%
        dplyr::pull(Variable)
    }

    dropped_by_variance_filter <- setdiff(VarsToReduce, VarsKeep)

    DataSubset <- DataSubset[VarsKeep]
    VarsToReduce <- VarsKeep
  }

  if (length(VarsToReduce) < 2) {
    stop("Fewer than 2 PCA variables remain after variance filtering.")
  }

  participant_missingness_final_vars <- rowMeans(is.na(DataSubset))
  complete_case_rows <- which(stats::complete.cases(DataSubset))
  incomplete_case_rows <- which(!stats::complete.cases(DataSubset))

  rows_imputed <- integer(0)
  rows_excluded_missing <- integer(0)
  rows_too_missing_to_impute <- integer(0)
  rows_all_missing <- which(participant_missingness_final_vars == 1)

  # Handle missing data

  if (MissingDataStrategy == "stop" && anyNA(DataSubset)) {
    stop(
      "MissingDataStrategy = 'stop' but missing values were found in PCA variables. ",
      "Use MissingDataStrategy = 'complete_cases' to score only complete rows, ",
      "or MissingDataStrategy = 'impute' to impute eligible rows."
    )
  }

  if (MissingDataStrategy == "complete_cases") {
    rows_used_for_pca <- complete_case_rows
    rows_excluded_missing <- incomplete_case_rows

    if (length(rows_excluded_missing) > 0 && !SuppressWarnings) {
      warning(
        length(rows_excluded_missing),
        " rows have missing values in PCA variables and will receive NA PCA scores. ",
        "Set MissingDataStrategy = 'impute' to impute eligible rows."
      )
    }

    DataForPCA <- DataSubset[rows_used_for_pca, , drop = FALSE]
  }

  if (MissingDataStrategy == "stop") {
    rows_used_for_pca <- seq_len(nrow(DataSubset))
    DataForPCA <- DataSubset
  }

  if (MissingDataStrategy == "impute") {
    rows_eligible_for_imputation <- which(
      participant_missingness_final_vars <= MaxMissingForImputation
    )

    if (!ImputeRowsWithAllMissing) {
      rows_eligible_for_imputation <- setdiff(
        rows_eligible_for_imputation,
        rows_all_missing
      )
    }

    rows_used_for_pca <- rows_eligible_for_imputation
    rows_excluded_missing <- setdiff(seq_len(nrow(DataSubset)), rows_used_for_pca)
    rows_too_missing_to_impute <- which(
      participant_missingness_final_vars > MaxMissingForImputation
    )

    if (!ImputeRowsWithAllMissing) {
      rows_too_missing_to_impute <- union(
        rows_too_missing_to_impute,
        rows_all_missing
      )
    }

    rows_too_missing_to_impute <- sort(unique(rows_too_missing_to_impute))

    DataForPCA <- DataSubset[rows_used_for_pca, , drop = FALSE]

    rows_imputed <- rows_used_for_pca[
      !stats::complete.cases(DataSubset[rows_used_for_pca, , drop = FALSE])
    ]

    if (length(rows_excluded_missing) > 0 && !SuppressWarnings) {
      warning(
        length(rows_excluded_missing),
        " rows exceeded the allowed missingness threshold for imputation and will receive NA PCA scores. ",
        "Current MaxMissingForImputation = ",
        MaxMissingForImputation,
        "."
      )
    }

    if (anyNA(DataForPCA)) {
      if (ImputeMethod == "missRanger") {
        if (!requireNamespace("missRanger", quietly = TRUE)) {
          stop(
            "The missRanger package is required when ImputeMethod = 'missRanger'. ",
            "Install it with install.packages('missRanger') or use ImputeMethod = 'median'."
          )
        }

        set.seed(123456)

        original_names <- colnames(DataForPCA)
        safe_names <- make.names(original_names, unique = TRUE)

        colnames(DataForPCA) <- safe_names

        DataForPCA <- missRanger::missRanger(
          data = DataForPCA,
          num.trees = 100
        )

        colnames(DataForPCA) <- original_names
      }

      if (ImputeMethod == "median") {
        DataForPCA <- DataForPCA %>%
          dplyr::mutate(
            dplyr::across(
              dplyr::everything(),
              ~ {
                med <- stats::median(.x, na.rm = TRUE)

                if (is.na(med)) {
                  return(.x)
                }

                dplyr::if_else(is.na(.x), med, .x)
              }
            )
          )
      }
    }
  }

  if (nrow(DataForPCA) < 2) {
    stop(
      "Fewer than 2 rows are available for PCA after missing-data handling. ",
      "Consider using MissingDataStrategy = 'impute', increasing MaxMissingForImputation, ",
      "or reducing VarsToReduce."
    )
  }

  if (anyNA(DataForPCA)) {
    vars_still_missing <- names(which(colSums(is.na(DataForPCA)) > 0))

    stop(
      "Missing values remain after missing-data handling in these PCA variables: ",
      paste(vars_still_missing, collapse = ", "),
      ". PCA cannot proceed with missing values."
    )
  }

  # Drop zero or undefined variance variables after missing-data handling

  raw_sds <- vapply(
    DataForPCA,
    stats::sd,
    numeric(1),
    na.rm = TRUE
  )

  zero_sd_vars <- names(raw_sds)[raw_sds == 0 | is.na(raw_sds)]

  if (length(zero_sd_vars) > 0) {
    if (!SuppressWarnings) {
      warning(
        "Dropping PCA variables with zero or undefined standard deviation after preprocessing: ",
        paste(zero_sd_vars, collapse = ", ")
      )
    }

    VarsToReduce <- setdiff(VarsToReduce, zero_sd_vars)
    DataSubset <- DataSubset[, VarsToReduce, drop = FALSE]
    DataForPCA <- DataForPCA[, VarsToReduce, drop = FALSE]
  }

  if (length(VarsToReduce) < 2) {
    stop("Fewer than 2 valid PCA variables remain after preprocessing.")
  }

  # Center and scale once

  if (center || scale) {
    DataSubset_scaled <- scale(
      DataForPCA,
      center = center,
      scale = scale
    )

    scale_means <- attr(DataSubset_scaled, "scaled:center")
    scale_sds <- attr(DataSubset_scaled, "scaled:scale")

    if (is.null(scale_means)) {
      scale_means <- rep(0, length(VarsToReduce))
      names(scale_means) <- VarsToReduce
    }

    if (is.null(scale_sds)) {
      scale_sds <- rep(1, length(VarsToReduce))
      names(scale_sds) <- VarsToReduce
    }
  } else {
    DataSubset_scaled <- as.matrix(DataForPCA)

    scale_means <- rep(0, length(VarsToReduce))
    scale_sds <- rep(1, length(VarsToReduce))
    names(scale_means) <- VarsToReduce
    names(scale_sds) <- VarsToReduce
  }

  post_scale_bad_vars <- names(scale_sds)[scale_sds == 0 | is.na(scale_sds)]

  if (length(post_scale_bad_vars) > 0) {
    stop(
      "Internal PCA preprocessing error. These variables still have zero or undefined SD after filtering: ",
      paste(post_scale_bad_vars, collapse = ", ")
    )
  }

  DataSubset_scaled <- as.matrix(DataSubset_scaled)

  n_vars <- length(VarsToReduce)

  # Run PCA

  if (Mode == "classic") {
    fit1 <- psych::principal(
      r = DataSubset_scaled,
      nfactors = n_vars,
      rotate = "none",
      scores = FALSE
    )

    Vaccounted <- as.data.frame(t(fit1$Vaccounted))

    if (is.null(numComponents)) {
      nc_candidates <- which(Vaccounted$`Cumulative Proportion` > minThresh)

      if (length(nc_candidates) == 0) {
        nc <- ncol(Vaccounted)
      } else {
        nc <- min(nc_candidates)
      }
    } else {
      nc <- min(numComponents, n_vars)
    }

    Vaccounted$Component <- factor(
      rownames(Vaccounted),
      levels = rownames(Vaccounted)
    )

    p_scree <- ggplot2::ggplot(
      Vaccounted,
      ggplot2::aes(x = Component)
    ) +
      ggplot2::geom_line(
        ggplot2::aes(y = `Cumulative Var`, group = 1)
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = `Cumulative Proportion`)
      ) +
      ggplot2::geom_col(
        ggplot2::aes(y = `Proportion Var`)
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = minThresh),
        linetype = "dashed"
      ) +
      ggplot2::theme_bw()

    fit <- psych::principal(
      r = DataSubset_scaled,
      nfactors = nc,
      rotate = rotate,
      scores = TRUE
    )

    fit <- psych::fa.sort(fit)

    LoadingTable <- as.data.frame(
      xtable::xtable(unclass(fit$loadings))
    )

    LoadingTable$Variable <- rownames(LoadingTable)

    scores <- fit$scores
  }

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
        r = DataSubset_scaled,
        nfactors = scree_n,
        rotate = "none",
        scores = FALSE
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

    p_scree <- ggplot2::ggplot(
      Vaccounted,
      ggplot2::aes(x = Component)
    ) +
      ggplot2::geom_line(
        ggplot2::aes(y = `Cumulative Var`, group = 1)
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = `Cumulative Proportion`)
      ) +
      ggplot2::geom_col(
        ggplot2::aes(y = `Proportion Var`)
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = minThresh),
        linetype = "dashed"
      ) +
      ggplot2::theme_bw()

    if (is.null(numComponents)) {
      nc_candidates <- which(Vaccounted$`Cumulative Proportion` > minThresh)

      if (length(nc_candidates) == 0) {
        nc_from_thresh <- min(maxComponents, scree_n)
      } else {
        nc_from_thresh <- min(nc_candidates)
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
        r = DataSubset_scaled,
        nfactors = nc,
        rotate = rotate,
        scores = TRUE
      )

      fit <- psych::fa.sort(fit)

      LoadingTable <- as.data.frame(
        xtable::xtable(unclass(fit$loadings))
      )

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

      LoadingTable <- as.data.frame(
        xtable::xtable(unclass(loadings_final))
      )

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

  # Apply labels

  LoadingTable <- dplyr::left_join(
    LoadingTable,
    LabelCodebook,
    by = "Variable",
    multiple = "first"
  )

  LoadingTable <- LoadingTable %>%
    dplyr::mutate(
      Labels = dplyr::if_else(
        is.na(.data$Labels) | .data$Labels == "",
        .data$Variable,
        .data$Labels
      ),
      PlotLabel = make.unique(as.character(.data$Labels))
    )

  if (!Relabel) {
    LoadingTable <- LoadingTable %>%
      dplyr::mutate(
        Labels = .data$Variable,
        PlotLabel = .data$Variable
      )
  }

  # Build outputs

  ScoresUsedRows <- as.data.frame(scores)

  rc_score_cols <- colnames(ScoresUsedRows)

  Scores <- as.data.frame(
    matrix(
      NA_real_,
      nrow = nrow(Data),
      ncol = length(rc_score_cols)
    )
  )

  colnames(Scores) <- rc_score_cols

  Scores[rows_used_for_pca, rc_score_cols] <- ScoresUsedRows

  CombinedData <- cbind(Data, Scores)

  rc_cols <- names(LoadingTable)[grepl("^RC", names(LoadingTable))]

  meltedLoading <- LoadingTable %>%
    tidyr::pivot_longer(cols = dplyr::all_of(rc_cols))

  meltedLoading$name <- factor(
    meltedLoading$name,
    levels = rc_cols
  )

  meltedLoading <- meltedLoading %>%
    dplyr::group_by(.data$Variable) %>%
    dplyr::arrange(.data$value, .by_group = TRUE) %>%
    dplyr::mutate(
      TopContributor = abs(.data$value) > 0.4
    ) %>%
    dplyr::ungroup()

  if (is.null(VariableCategoryTable)) {
    meltedLoading$color <- "black"
  } else {
    meltedLoading <- meltedLoading %>%
      dplyr::left_join(
        VariableCategoryTable,
        by = "Variable"
      ) %>%
      dplyr::mutate(
        color = .data$Category
      )
  }

  meltedLoading$PlotLabel <- factor(
    meltedLoading$PlotLabel,
    levels = rev(unique(LoadingTable$PlotLabel))
  )

  if (is.null(VariableCategoryTable)) {
    p <- ggplot2::ggplot(
      meltedLoading,
      ggplot2::aes(
        x = PlotLabel,
        y = value,
        alpha = TopContributor,
        color = color
      )
    ) +
      ggplot2::coord_flip() +
      ggplot2::geom_segment(
        linewidth = 2,
        ggplot2::aes(
          x = PlotLabel,
          xend = PlotLabel,
          y = 0,
          yend = value
        )
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(name),
        nrow = 1
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.y = ggplot2::element_text(angle = 45),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = c("black"))
  } else {
    p <- ggplot2::ggplot(
      meltedLoading,
      ggplot2::aes(
        x = PlotLabel,
        y = value,
        alpha = TopContributor,
        color = color
      )
    ) +
      ggplot2::coord_flip() +
      ggplot2::geom_segment(
        linewidth = 2,
        ggplot2::aes(
          x = PlotLabel,
          xend = PlotLabel,
          y = 0,
          yend = value
        )
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(name),
        nrow = 1
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.y = ggplot2::element_text(angle = 45),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = classcolors)
  }

  ScaleParams <- list(
    means = scale_means,
    sds = scale_sds,
    center = center,
    scale = scale
  )

  Preprocessing <- list(
    OriginalVariables = VarsOriginal,
    DroppedNonNumeric = dropped_non_numeric,
    DroppedByVarianceFilter = dropped_by_variance_filter,
    DroppedZeroVariance = zero_sd_vars,
    FinalVariablesUsed = VarsToReduce,
    NumOriginalVariables = length(VarsOriginal),
    NumFinalVariables = length(VarsToReduce),
    VariableMissingness = variable_missingness,
    ParticipantMissingness = participant_missingness,
    ParticipantMissingnessFinalVars = participant_missingness_final_vars,
    HighMissingVariables = high_missing_vars,
    HighMissingParticipants = high_missing_participants,
    MissingDataStrategy = MissingDataStrategy,
    ImputeMethod = ifelse(MissingDataStrategy == "impute", ImputeMethod, NA_character_),
    MaxMissingForImputation = ifelse(
      MissingDataStrategy == "impute",
      MaxMissingForImputation,
      NA_real_
    ),
    ImputeRowsWithAllMissing = ImputeRowsWithAllMissing,
    CompleteCaseRows = complete_case_rows,
    IncompleteCaseRows = incomplete_case_rows,
    RowsUsedForPCA = rows_used_for_pca,
    RowsExcludedFromPCA = rows_excluded_missing,
    RowsImputed = rows_imputed,
    RowsTooMissingToImpute = rows_too_missing_to_impute,
    RowsAllMissing = rows_all_missing,
    NumRowsUsedForPCA = length(rows_used_for_pca),
    NumRowsExcludedFromPCA = length(rows_excluded_missing),
    NumRowsImputed = length(rows_imputed)
  )

  # Return result

  return(list(
    p_scree = p_scree,
    pcaresults = fit,
    LoadingTable = LoadingTable,
    Scores = Scores,
    CombinedData = CombinedData,
    Lollipop = p,
    ScaleParams = ScaleParams,
    VarsUsed = VarsToReduce,
    VarianceTable = VarianceTable,
    Preprocessing = Preprocessing,
    Mode = Mode,
    Backend = backend,
    Center = center,
    Scale = scale
  ))
}

#' Create PCA table and visualization
#'
#' Compatibility alias for [CreatePCAObject()]. Prefer `CreatePCAObject()` in
#' new code because this workflow returns a reusable PCA object, not only a
#' static table.
#'
#' @param ... Arguments passed to [CreatePCAObject()].
#'
#' @return The same object returned by [CreatePCAObject()].
#'
#' @examples
#' PCA <- CreatePCATable(
#'   Data = mtcars,
#'   VarsToReduce = names(mtcars),
#'   numComponents = 3
#' )
#'
#' @export
CreatePCATable <- function(...) {
  CreatePCAObject(...)
}
