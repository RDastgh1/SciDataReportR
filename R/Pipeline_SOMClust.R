#' SOM + latent profile clustering pipeline (with AHP and distance baselines)
#'
#' @description
#' End-to-end pipeline to:
#' - Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied
#'   Z-score object.
#' - Fit a Self-Organizing Map (SOM; kohonen) on complete cases.
#' - Generate aweSOM visualizations (Circular, Line, Cloud) with optional
#'   relabeling using variable labels from the original data frame.
#' - Cluster SOM codebook vectors using latent profile analysis (tidyLPA /
#'   mclust backend).
#' - In \code{method = "exploratory"}, fit a grid of models and select a
#'   recommended solution using an Analytic Hierarchy Process (AHP)-style
#'   index combining AIC, BIC, and Entropy.
#' - In \code{method = "finalize"}, fit a user-specified model and number
#'   of profiles.
#' - Map node-level clusters and posterior probabilities back to individuals.
#' - Store training variable summaries used later to quantify whether projected
#'   cohorts fall outside the original training range.
#'
#' This supports a train once, project many clinical phenotyping workflow: learn
#' phenotype structure in a training cohort, then project new cohorts into the
#' fixed phenotype space without reclustering.
#'
#' Missing data:
#' - SOM and clustering are fit only on rows with complete Z-scores.
#' - The returned \code{df_with_clusters} has the full original data with
#'   \code{.scidr_rowid} and a single cluster column appended; rows not used
#'   in SOM/LPA get NA.
#' - The returned \code{ProbFit$individual} is also full length, preserving
#'   one row per input row with NA posterior probabilities for rows excluded
#'   from SOM/LPA.
#'
#' Stable row id:
#' - \code{.scidr_rowid} is added to the input data and carried into
#'   \code{df_with_clusters} and \code{ProbFit$individual}.
#' - \code{ProbFit$individual$RowID} is set equal to \code{.scidr_rowid} so
#'   merges do not rely on row order.
#'
#' Z-score behavior:
#' - \code{ZScoreType = "Center and Scale"/"Center Only"/"Scale Only"} computes
#'   Z-scores from \code{df} via \code{CreateZScoreObject()}.
#' - \code{ZScoreType = "ZScoreObj"} projects Z-scores using an external
#'   \code{ZScoreObj} via \code{ProjectZScore()}.
#' - \code{ZScoreType = "PreZScored"} uses existing Z-score columns in \code{df}
#'   as-is and does not re-zscore.
#'
#' @param df Data frame containing the variables to be used in SOM and
#'   clustering.
#' @param variables Optional character vector of variable names. If NULL,
#'   numeric variables are auto-detected using
#'   \code{SciDataReportR::getNumVars(df, Ordinal = FALSE)}. In
#'   \code{ZScoreType = "PreZScored"}, this can also be NULL if you supply
#'   \code{ZScoreVars} or if Z-score columns can be auto-detected by prefix.
#' @param method One of \code{"exploratory"} (default) or \code{"finalize"}.
#'   In \code{"exploratory"}, a grid of models is fit and AHP chooses the
#'   recommended solution. In \code{"finalize"}, the user must specify
#'   \code{final_k} and \code{final_model}.
#' @param k_range Integer vector of numbers of clusters/profiles to consider
#'   in exploratory mode. Default \code{2:15}.
#' @param models Integer vector of model specifications for tidyLPA
#'   (mclust backend). Default \code{c(1, 2, 3)}.
#' @param final_k Integer; number of profiles for \code{method = "finalize"}.
#' @param final_model Integer; model specification for \code{method = "finalize"}
#'   (should be one of \code{models}).
#' @param ClusterName Name of the cluster column in the output. Defaults to
#'   \code{"Cluster"}. If this column already exists in \code{df}, it is
#'   overwritten (with a message).
#' @param ZScoreType One of:
#'   \itemize{
#'     \item \code{"Center and Scale"} (default)
#'     \item \code{"Center Only"}
#'     \item \code{"Scale Only"}
#'     \item \code{"ZScoreObj"} (use an existing ZScore object)
#'     \item \code{"PreZScored"} (use existing Z-score columns in df as-is)
#'   }
#' @param ZScoreObject Optional ZScoreObj (from \code{CreateZScoreObject()} or
#'   \code{ProjectZScore()}) to use when \code{ZScoreType = "ZScoreObj"}.
#' @param som_xdim,som_ydim Optional integers for SOM grid dimensions. If NULL,
#'   a square grid with side length \code{ceiling(n_complete^(1/3))} is used.
#' @param som_topo SOM topology for \code{kohonen::somgrid()}, default
#'   \code{"hexagonal"}.
#' @param som_neigh SOM neighbourhood function, default \code{"gaussian"}.
#' @param seed_som,seed_lpa Integer seeds for SOM and LPA steps
#'   (defaults 934521 and 93421).
#' @param Relabel Logical; if TRUE (default), aweSOM plots are relabeled
#'   using variable labels from the *original* \code{df} (via Hmisc or
#'   sjlabelled when available) by stripping the Z-score prefix.
#' @param ZScorePrefix Character prefix used for Z-score columns when
#'   \code{ZScoreType = "PreZScored"}. Default \code{"Z_"}.
#' @param ZScoreVars Optional character vector of Z-score column names to use
#'   when \code{ZScoreType = "PreZScored"}. If NULL, the function attempts to
#'   infer them from \code{variables} or by detecting columns starting with
#'   \code{ZScorePrefix}.
#' @param id_col Optional character scalar. If provided and present in \code{df},
#'   this column is carried into \code{ProbFit$individual} for convenience.
#' @param lpa_progress Logical; if TRUE, print short progress messages while
#'   fitting model/profile combinations.
#' @param lpa_em_itmax Integer; maximum number of EM iterations passed to
#'   \code{mclust::emControl()}. Use NULL to leave mclust defaults unchanged.
#' @param lpa_em_tol Numeric; EM convergence tolerance passed to
#'   \code{mclust::emControl()}. Use NULL to leave mclust defaults unchanged.
#' @param lpa_drop_zero_sd Logical; if TRUE, remove SOM code dimensions with
#'   near-zero standard deviation before LPA.
#' @param lpa_zero_sd_tol Numeric tolerance used when \code{lpa_drop_zero_sd = TRUE}.
#' @param lpa_timeout_seconds Optional timeout in seconds for individual LPA
#'   fits. Use NULL to disable timeouts.
#' @param skip_model_after_n_failures Optional integer; skip a model family
#'   after this many failures.
#' @param slow_fit_seconds Optional runtime threshold used to flag slow LPA
#'   fits in diagnostics.
#' @param min_nodes_per_cluster Optional minimum average SOM nodes per cluster
#'   considered before attempting a candidate profile count.
#' @param high_dist_quantile Numeric value between 0 and 1 used to define
#'   high SOM-distance flags from the training distance distribution. Default
#'   is \code{0.95}.
#' @param low_prob_threshold Numeric posterior probability threshold used to
#'   flag uncertain phenotype membership. Default is \code{0.70}.
#'
#' @details
#' The AHP-style index is computed by:
#' \enumerate{
#'   \item Scaling AIC, BIC, and Entropy across candidate solutions
#'         (AIC/BIC are negated so that lower values correspond to better
#'         fit; higher scaled scores are preferred).
#'   \item Taking the mean of the three scaled indices. The model with the
#'         highest AHP index is recommended.
#' }
#'
#' LPA model/profile combinations are fit one at a time so that failed or
#' warning-producing solutions are captured in diagnostics instead of blocking
#' the entire pipeline. Successful fits are retained and failed fits are listed
#' in \code{ModelInfo_MClust$diagnostics}.
#'
#' @return A list of class \code{"Pipeline_SOMClust"} with components:
#'   \itemize{
#'     \item \code{method}, \code{vars_used}, \code{ZScoreType},
#'           \code{ZScoreObject}, \code{ZScoreVars}, \code{ClusterName}
#'     \item \code{complete_rows}: logical vector (rows used for SOM/LPA)
#'     \item \code{df_with_clusters}: original \code{df} with \code{.scidr_rowid}
#'           and only the cluster column appended
#'     \item \code{fit_plot}: ggplot of AIC/BIC/Entropy vs k and model
#'     \item \code{ModelInfo_SOM}: list with \code{som_model},
#'           \code{som_codes}, \code{som_grid}, \code{training_variable_summary},
#'           \code{SOMFit} (distance diagnostics, baselines, and per-cluster
#'           flags), \code{plots} (aweSOM plots)
#'     \item \code{ModelInfo_MClust}: list with \code{lpa_models},
#'           \code{fit_table}, \code{AHP} information, and \code{diagnostics}
#'           for LPA warnings, failures, runtimes, and preprocessing
#'     \item \code{ProbFit}: list with \code{node} (node-level posterior
#'           probabilities), \code{individual} (full-length per-person mapping
#'           and probabilities, including \code{.scidr_rowid}), and probability plots
#'   }
#'
#' @references
#' Saaty TL. \emph{The Analytic Hierarchy Process}. McGraw-Hill, 1980.
#'
#' @export
CreateSOMClusterModel <- function(
    df,
    variables       = NULL,
    method          = c("exploratory", "finalize"),
    k_range         = 2:15,
    models          = c(1, 2, 3),
    final_k         = NULL,
    final_model     = NULL,
    ClusterName     = "Cluster",
    ZScoreType      = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj", "PreZScored"),
    ZScoreObject    = NULL,
    som_xdim        = NULL,
    som_ydim        = NULL,
    som_topo        = "hexagonal",
    som_neigh       = "gaussian",
    seed_som        = 934521L,
    seed_lpa        = 93421L,
    Relabel         = TRUE,
    ZScorePrefix    = "Z_",
    ZScoreVars      = NULL,
    id_col          = NULL,
    lpa_progress     = FALSE,
    lpa_em_itmax     = 100L,
    lpa_em_tol       = 1e-5,
    lpa_timeout_seconds = 120,
    lpa_drop_zero_sd = TRUE,
    lpa_zero_sd_tol  = 1e-8,
    skip_model_after_n_failures = 2L,
    slow_fit_seconds = 120,
    min_nodes_per_cluster = 5,
    high_dist_quantile = 0.95,
    low_prob_threshold = 0.70
) {

  method     <- match.arg(method)
  ZScoreType <- match.arg(ZScoreType)

  if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
    stop("SciDataReportR must be installed.")
  }
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required.")
  }
  if (!requireNamespace("aweSOM", quietly = TRUE)) {
    stop("Package 'aweSOM' is required.")
  }
  if (!requireNamespace("tidyLPA", quietly = TRUE)) {
    stop("Package 'tidyLPA' is required.")
  }
  if ((!is.null(lpa_em_itmax) || !is.null(lpa_em_tol)) &&
      !requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required when lpa_em_itmax or lpa_em_tol is supplied.")
  }
  if (!is.null(lpa_timeout_seconds) &&
      !requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required when lpa_timeout_seconds is supplied. Install it or set lpa_timeout_seconds = NULL.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr', 'tidyr', and 'ggplot2' are required.")
  }

  if (!is.logical(lpa_progress) || length(lpa_progress) != 1) {
    stop("lpa_progress must be TRUE or FALSE.")
  }
  if (!is.logical(lpa_drop_zero_sd) || length(lpa_drop_zero_sd) != 1) {
    stop("lpa_drop_zero_sd must be TRUE or FALSE.")
  }
  if (!is.numeric(lpa_zero_sd_tol) || length(lpa_zero_sd_tol) != 1 || lpa_zero_sd_tol < 0) {
    stop("lpa_zero_sd_tol must be a single non-negative numeric value.")
  }
  if (!is.null(lpa_em_itmax) &&
      (!is.numeric(lpa_em_itmax) || length(lpa_em_itmax) != 1 || lpa_em_itmax < 1)) {
    stop("lpa_em_itmax must be NULL or a single positive integer.")
  }
  if (!is.null(lpa_em_tol) &&
      (!is.numeric(lpa_em_tol) || length(lpa_em_tol) != 1 || lpa_em_tol <= 0)) {
    stop("lpa_em_tol must be NULL or a single positive numeric value.")
  }
  if (!is.null(lpa_timeout_seconds) &&
      (!is.numeric(lpa_timeout_seconds) ||
       length(lpa_timeout_seconds) != 1 ||
       lpa_timeout_seconds <= 0)) {
    stop("lpa_timeout_seconds must be NULL or a single positive numeric value.")
  }
  if (!is.null(skip_model_after_n_failures) &&
      (!is.numeric(skip_model_after_n_failures) ||
       length(skip_model_after_n_failures) != 1 ||
       skip_model_after_n_failures < 1)) {
    stop("skip_model_after_n_failures must be NULL or a single positive integer.")
  }
  if (!is.null(slow_fit_seconds) &&
      (!is.numeric(slow_fit_seconds) || length(slow_fit_seconds) != 1 || slow_fit_seconds <= 0)) {
    stop("slow_fit_seconds must be NULL or a single positive numeric value.")
  }
  if (!is.null(min_nodes_per_cluster) &&
      (!is.numeric(min_nodes_per_cluster) ||
       length(min_nodes_per_cluster) != 1 ||
       min_nodes_per_cluster <= 0)) {
    stop("min_nodes_per_cluster must be NULL or a single positive numeric value.")
  }
  if (!is.numeric(high_dist_quantile) || length(high_dist_quantile) != 1 ||
      is.na(high_dist_quantile) || high_dist_quantile <= 0 || high_dist_quantile >= 1) {
    stop("high_dist_quantile must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(low_prob_threshold) || length(low_prob_threshold) != 1 ||
      is.na(low_prob_threshold) || low_prob_threshold < 0 || low_prob_threshold > 1) {
    stop("low_prob_threshold must be a single numeric value between 0 and 1.")
  }

  # Stable row id ----------------------------------------------------------

  df_scidr <- df %>%
    dplyr::mutate(.scidr_rowid = dplyr::row_number())

  if (!is.null(id_col)) {
    if (!id_col %in% names(df_scidr)) {
      stop("id_col '", id_col, "' was provided but is not in df.")
    }
  }

  # Variables --------------------------------------------------------------

  ZScoreVars_used <- NULL
  vars_used       <- NULL

  if (ZScoreType == "PreZScored") {

    if (!is.null(ZScoreVars)) {
      ZScoreVars_used <- ZScoreVars
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)

    } else if (!is.null(variables)) {

      if (all(variables %in% names(df_scidr)) &&
          all(grepl(paste0("^", ZScorePrefix), variables))) {
        ZScoreVars_used <- variables
        vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
      } else {
        ZScoreVars_used <- paste0(ZScorePrefix, variables)
        vars_used <- variables
      }

    } else {

      ZScoreVars_used <- grep(paste0("^", ZScorePrefix), names(df_scidr), value = TRUE)
      if (length(ZScoreVars_used) == 0) {
        stop("ZScoreType = 'PreZScored' requires ZScoreVars, variables, or columns starting with '",
             ZScorePrefix, "'.")
      }
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
    }

    missing_z <- setdiff(ZScoreVars_used, names(df_scidr))
    if (length(missing_z) > 0) {
      stop("The following pre-zscored variables are not in df: ",
           paste(missing_z, collapse = ", "))
    }

    z_tmp <- df_scidr[, ZScoreVars_used, drop = FALSE]
    is_num_z <- vapply(z_tmp, is.numeric, logical(1))
    if (!all(is_num_z)) {
      warning("Dropping non-numeric pre-zscored variables: ",
              paste(names(z_tmp)[!is_num_z], collapse = ", "))
      ZScoreVars_used <- ZScoreVars_used[is_num_z]
      z_tmp <- z_tmp[, is_num_z, drop = FALSE]
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
    }

    if (length(ZScoreVars_used) == 0) {
      stop("No numeric pre-zscored variables available for SOM / clustering.")
    }

  } else {

    if (is.null(variables)) {
      variables <- SciDataReportR::getNumVars(df, Ordinal = FALSE)
    }

    missing_vars <- setdiff(variables, names(df))
    if (length(missing_vars) > 0) {
      stop("The following variables are not in df: ",
           paste(missing_vars, collapse = ", "))
    }

    is_num <- vapply(df[variables], is.numeric, logical(1))
    if (!all(is_num)) {
      warning("Dropping non-numeric variables: ",
              paste(variables[!is_num], collapse = ", "))
      variables <- variables[is_num]
    }
    if (length(variables) == 0) {
      stop("No numeric variables available for SOM / clustering.")
    }

    vars_used <- variables
  }



  # Training variable summary ---------------------------------------------

  training_variable_summary <- dplyr::tibble(
    Variable = vars_used
  ) %>%
    dplyr::mutate(
      present_in_training = Variable %in% names(df_scidr),
      is_numeric = vapply(Variable, function(v) {
        v %in% names(df_scidr) && is.numeric(df_scidr[[v]])
      }, logical(1)),
      n_training = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_integer_ else length(df_scidr[[v]])
      }, integer(1)),
      n_missing = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_integer_ else sum(is.na(df_scidr[[v]]))
      }, integer(1)),
      prop_missing = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_real_ else mean(is.na(df_scidr[[v]]))
      }, numeric(1)),
      mean = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]])) NA_real_ else mean(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      sd = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]])) NA_real_ else stats::sd(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      min = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]]) || all(is.na(df_scidr[[v]]))) NA_real_ else min(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]]) || all(is.na(df_scidr[[v]]))) NA_real_ else max(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1))
    )

  # Z-scores ---------------------------------------------------------------

  if (ZScoreType == "PreZScored") {

    z_df <- df_scidr[, ZScoreVars_used, drop = FALSE]
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows in pre-zscored columns; cannot fit SOM.")
    }

    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreObject_used <- NULL

  } else if (ZScoreType == "ZScoreObj") {

    if (is.null(ZScoreObject) || !"ZScoreObj" %in% class(ZScoreObject)) {
      stop("ZScoreType = 'ZScoreObj' requires a valid ZScoreObj.")
    }

    z_res <- SciDataReportR::ProjectZScore(
      df                 = df_scidr,
      variables          = vars_used,
      parameters         = ZScoreObject,
      ParameterInputType = "ZScoreObj",
      names_prefix       = "Z_"
    )
    ZScoreObject_used <- z_res

    z_df <- ZScoreObject_used$ZScores
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows after Z-score step; cannot fit SOM.")
    }

    z_df <- z_df[, names(z_df), drop = FALSE]
    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreVars_used <- names(z_df)

  } else {

    center_flag <- ZScoreType %in% c("Center and Scale", "Center Only")
    scale_flag  <- ZScoreType %in% c("Center and Scale", "Scale Only")

    ZScoreObject_used <- SciDataReportR::CreateZScoreObject(
      df           = df_scidr,
      variables    = vars_used,
      names_prefix = "Z_",
      center       = center_flag,
      scale        = scale_flag
    )

    z_df <- ZScoreObject_used$ZScores
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows after Z-score step; cannot fit SOM.")
    }

    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreVars_used <- names(z_df)
  }

  # SOM fitting ------------------------------------------------------------

  n_complete <- nrow(zmat)
  som_grid_initial_xdim <- som_xdim
  som_grid_initial_ydim <- som_ydim

  if (is.null(som_xdim) || is.null(som_ydim)) {
    side <- ceiling(n_complete^(1 / 3))
    som_xdim <- side
    som_ydim <- side
  }

  som_grid <- kohonen::somgrid(
    xdim              = som_xdim,
    ydim              = som_ydim,
    topo              = som_topo,
    neighbourhood.fct = som_neigh
  )

  set.seed(seed_som)
  som_model <- kohonen::som(
    X         = zmat,
    grid      = som_grid,
    keep.data = TRUE
  )

  codes_raw <- som_model$codes
  if (is.list(codes_raw)) {
    som_codes <- as.matrix(codes_raw[[1]])
  } else {
    som_codes <- as.matrix(codes_raw)
  }

  # SOM distances (distance to BMU) ----------------------------------------

  SOM_Node_full <- rep(NA_integer_, nrow(df_scidr))
  SOM_Node_full[complete_rows] <- som_model$unit.classif

  SOM_Dist_full <- rep(NA_real_, nrow(df_scidr))
  SOM_Dist_full[complete_rows] <- som_model$distances

  train_quant_err <- mean(som_model$distances, na.rm = TRUE)

  som_node_occupancy <- dplyr::tibble(
    NodeID = seq_len(nrow(som_codes))
  ) %>%
    dplyr::left_join(
      as.data.frame(table(som_model$unit.classif)) %>%
        dplyr::transmute(
          NodeID = as.integer(as.character(Var1)),
          n = as.integer(Freq)
        ),
      by = "NodeID"
    ) %>%
    dplyr::mutate(n = dplyr::if_else(is.na(n), 0L, n))

  som_occupancy_summary <- dplyr::tibble(
    n_nodes = nrow(som_node_occupancy),
    n_empty_nodes = sum(som_node_occupancy$n == 0),
    n_singleton_nodes = sum(som_node_occupancy$n == 1),
    min_occupancy = min(som_node_occupancy$n),
    median_occupancy = stats::median(som_node_occupancy$n),
    max_occupancy = max(som_node_occupancy$n)
  )

  # AweSOM plots -----------------------------------------------------------

  CircPlot  <- aweSOM::aweSOMplot(som = som_model, type = "Circular", palvar = "viridis")
  LinePlot  <- aweSOM::aweSOMplot(som = som_model, type = "Line",     palvar = "viridis")
  CloudPlot <- aweSOM::aweSOMplot(som = som_model, type = "Cloud")

  if (Relabel) {
    relabel_fun <- function(w) {
      if (!is.null(w$x$label)) {
        zvars     <- w$x$label
        base_vars <- sub(paste0("^", ZScorePrefix), "", zvars)
        base_vars <- sub("^Z_", "", base_vars)

        get_lab <- NULL
        if (requireNamespace("Hmisc", quietly = TRUE)) {
          get_lab <- function(v) Hmisc::label(df[[v]])
        } else if (requireNamespace("sjlabelled", quietly = TRUE)) {
          get_lab <- function(v) sjlabelled::get_label(df[[v]])
        }

        if (!is.null(get_lab)) {
          new_labels <- vapply(base_vars, function(bv) {
            if (!bv %in% names(df)) return(bv)
            lab <- get_lab(bv)
            if (is.null(lab) || !nzchar(lab)) bv else as.character(lab)
          }, character(1))
        } else {
          new_labels <- base_vars
        }

        w$x$label <- unname(new_labels)
      }
      w
    }

    CircPlot  <- relabel_fun(CircPlot)
    LinePlot  <- relabel_fun(LinePlot)
  }

  som_plots <- list(
    Circular = CircPlot,
    Line     = LinePlot,
    Cloud    = CloudPlot
  )

  # LPA / mclust on SOM codes ----------------------------------------------

  X <- som_codes

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }

  dropped_lpa_vars <- character(0)
  if (lpa_drop_zero_sd) {
    lpa_sd <- apply(X, 2, stats::sd, na.rm = TRUE)
    keep_lpa_vars <- is.finite(lpa_sd) & lpa_sd > lpa_zero_sd_tol
    dropped_lpa_vars <- colnames(X)[!keep_lpa_vars]

    if (length(dropped_lpa_vars) > 0) {
      warning(
        "Dropping near-zero SD SOM code dimensions before LPA: ",
        paste(dropped_lpa_vars, collapse = ", ")
      )
      X <- X[, keep_lpa_vars, drop = FALSE]
    }

    if (ncol(X) == 0) {
      stop("No SOM code dimensions with non-zero variance are available for LPA.")
    }
  }

  X <- as.data.frame(X)
  set.seed(seed_lpa)

  lpa_control <- NULL
  if (!is.null(lpa_em_itmax) || !is.null(lpa_em_tol)) {
    lpa_control <- mclust::emControl(
      itmax = if (is.null(lpa_em_itmax)) NULL else as.integer(lpa_em_itmax),
      tol   = if (is.null(lpa_em_tol)) NULL else lpa_em_tol
    )
  }

  # This helper is intentionally kept inside the function because it isolates
  # warning capture, runtime tracking, and failed-fit handling for one fragile step.
  safe_lpa_fit <- function(X, k, model, lpa_control = NULL, lpa_timeout_seconds = NULL) {

    warnings_vec <- character(0)
    error_msg    <- NA_character_
    start_time   <- Sys.time()

    fit_expr <- function() {
      if (is.null(lpa_control)) {
        tidyLPA::estimate_profiles(
          X,
          n_profiles = k,
          models     = model
        )
      } else {
        tidyLPA::estimate_profiles(
          X,
          n_profiles = k,
          models     = model,
          control    = lpa_control
        )
      }
    }

    fit <- withCallingHandlers(
      tryCatch(
        {
          if (is.null(lpa_timeout_seconds)) {
            fit_expr()
          } else {
            R.utils::withTimeout(
              fit_expr(),
              timeout = lpa_timeout_seconds,
              onTimeout = "error"
            )
          }
        },
        TimeoutException = function(e) {
          error_msg <<- paste0("Timed out after ", lpa_timeout_seconds, " seconds.")
          NULL
        },
        error = function(e) {
          error_msg <<- conditionMessage(e)
          NULL
        }
      ),
      warning = function(w) {
        warnings_vec <<- c(warnings_vec, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    slow_fit <- !is.null(slow_fit_seconds) && runtime_seconds > slow_fit_seconds

    fit_info <- NULL
    fit_error_msg <- NA_character_
    if (!is.null(fit)) {
      fit_info <- tryCatch(
        tidyLPA::get_fit(fit),
        error = function(e) {
          fit_error_msg <<- conditionMessage(e)
          NULL
        }
      )
    }

    if (is.null(fit_info)) {
      if (is.na(error_msg) && !is.na(fit_error_msg)) {
        error_msg <- fit_error_msg
      }
      status <- "failed"
    } else {
      status <- "success"
    }

    diagnostics <- dplyr::tibble(
      Model           = as.integer(model),
      Classes         = as.integer(k),
      status          = status,
      runtime_seconds = runtime_seconds,
      slow_fit        = slow_fit,
      n_warnings      = length(warnings_vec),
      warnings        = paste(unique(warnings_vec), collapse = " | "),
      error           = error_msg
    )

    if (!is.null(fit_info)) {
      fit_info <- as.data.frame(fit_info)
      fit_info$Model <- as.integer(model)
      fit_info$Classes <- as.integer(k)
    }

    list(
      fit         = fit,
      fit_info    = fit_info,
      diagnostics = diagnostics
    )
  }

  # This helper prevents one available model from getting all of the AHP weight
  # simply because there is no variation in a fit index.
  scale_for_ahp <- function(x, higher_is_better = TRUE) {
    if (!higher_is_better) {
      x <- -x
    }
    if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) {
      return(rep(0, length(x)))
    }
    as.numeric(scale(x))
  }

  if (method == "exploratory") {

    lpa_models <- list()
    fit_rows   <- list()
    diag_rows  <- list()

    model_grid <- expand.grid(
      Model   = models,
      Classes = k_range,
      KEEP.OUT.ATTRS = FALSE
    )

    total_fits <- nrow(model_grid)
    model_problem_counts <- stats::setNames(rep(0L, length(unique(models))), unique(models))

    for (i in seq_len(total_fits)) {

      this_model <- model_grid$Model[i]
      this_k     <- model_grid$Classes[i]
      fit_name   <- paste0("model_", this_model, "_class_", this_k)

      if (!is.null(min_nodes_per_cluster) && nrow(X) / this_k < min_nodes_per_cluster) {
        warning(
          "Model ", this_model, ", k = ", this_k,
          " has fewer than ", min_nodes_per_cluster,
          " SOM nodes per requested profile on average. Interpret cautiously."
        )
      }

      if (!is.null(skip_model_after_n_failures) &&
          model_problem_counts[as.character(this_model)] >= skip_model_after_n_failures) {

        if (lpa_progress) {
          message(
            "Skipping LPA fit ", i, "/", total_fits,
            ": model ", this_model,
            ", k = ", this_k,
            " after repeated failed or slow fits."
          )
        }

        diag_rows[[fit_name]] <- dplyr::tibble(
          Model           = as.integer(this_model),
          Classes         = as.integer(this_k),
          status          = "skipped_model_family",
          runtime_seconds = NA_real_,
          slow_fit        = NA,
          n_warnings      = 0L,
          warnings        = NA_character_,
          error           = paste0(
            "Skipped after ", skip_model_after_n_failures,
            " failed or slow fits for model ", this_model, "."
          )
        )
        next
      }

      if (lpa_progress) {
        message(
          "LPA fit ", i, "/", total_fits,
          ": model ", this_model,
          ", k = ", this_k
        )
      }

      fit_res <- safe_lpa_fit(
        X           = X,
        k           = this_k,
        model       = this_model,
        lpa_control = lpa_control,
        lpa_timeout_seconds = lpa_timeout_seconds
      )

      diag_rows[[fit_name]] <- fit_res$diagnostics

      if (fit_res$diagnostics$status != "success" || isTRUE(fit_res$diagnostics$slow_fit)) {
        model_problem_counts[as.character(this_model)] <- model_problem_counts[as.character(this_model)] + 1L
      }

      if (!is.null(fit_res$fit) && !is.null(fit_res$fit_info)) {
        lpa_models[[fit_name]] <- fit_res$fit
        fit_rows[[fit_name]]   <- fit_res$fit_info
      }
    }

    lpa_diagnostics <- dplyr::bind_rows(diag_rows)
    failed_fits <- lpa_diagnostics %>%
      dplyr::filter(status != "success")

    warning_fits <- lpa_diagnostics %>%
      dplyr::filter(status == "success", n_warnings > 0)

    if (length(fit_rows) == 0) {
      stop(
        "No LPA models were successfully estimated. Inspect ModelInfo_MClust$diagnostics if available; ",
        "try fewer variables, smaller k_range, excluding unstable models, or increasing lpa_em_itmax."
      )
    }

    fit_table <- dplyr::bind_rows(fit_rows) %>%
      dplyr::mutate(
        Model = as.integer(Model),
        Classes = as.integer(Classes)
      )

    if (!"AIC" %in% names(fit_table)) fit_table$AIC <- NA_real_
    if (!"BIC" %in% names(fit_table)) fit_table$BIC <- NA_real_
    if (!"Entropy" %in% names(fit_table)) fit_table$Entropy <- NA_real_

    fit_table <- fit_table %>%
      dplyr::mutate(
        AIC_scaled     = scale_for_ahp(AIC, higher_is_better = FALSE),
        BIC_scaled     = scale_for_ahp(BIC, higher_is_better = FALSE),
        Entropy_scaled = scale_for_ahp(Entropy, higher_is_better = TRUE),
        ahp_index      = rowMeans(
          cbind(AIC_scaled, BIC_scaled, Entropy_scaled),
          na.rm = TRUE
        )
      )

    if (all(is.na(fit_table$ahp_index))) {
      stop("LPA models were estimated, but AHP could not be computed because fit indices are missing.")
    }

    best_idx     <- which.max(fit_table$ahp_index)
    ahp_best_row <- fit_table[best_idx, ]
    best_model   <- as.integer(ahp_best_row$Model)
    best_k       <- as.integer(ahp_best_row$Classes)
    best_name    <- paste0("model_", best_model, "_class_", best_k)
    best_fit_name <- best_name
    best_lpa     <- lpa_models[[best_name]]

    recommendation_txt <- paste0(
      "AHP (AIC, BIC, Entropy) recommends Model ",
      best_model, " with k = ", best_k, " profiles."
    )

    mdata_wide <- fit_table

  } else {

    if (is.null(final_k) || is.null(final_model)) {
      stop("For method = 'finalize', final_k and final_model must be supplied.")
    }

    if (lpa_progress) {
      message(
        "LPA finalize fit: model ", final_model,
        ", k = ", final_k
      )
    }

    fit_res <- safe_lpa_fit(
      X           = X,
      k           = final_k,
      model       = final_model,
      lpa_control = lpa_control
    )

    lpa_diagnostics <- fit_res$diagnostics
    failed_fits <- lpa_diagnostics %>%
      dplyr::filter(status != "success")
    warning_fits <- lpa_diagnostics %>%
      dplyr::filter(status == "success", n_warnings > 0)

    if (is.null(fit_res$fit) || is.null(fit_res$fit_info)) {
      stop(
        "The requested final LPA model could not be estimated: model ",
        final_model, ", k = ", final_k, "."
      )
    }

    lpa_models <- fit_res$fit
    best_lpa   <- lpa_models
    best_fit_name <- paste0("model_", final_model, "_class_", final_k)

    fit_table <- as.data.frame(fit_res$fit_info) %>%
      dplyr::mutate(
        Classes        = final_k,
        Model          = final_model,
        AIC_scaled     = NA_real_,
        BIC_scaled     = NA_real_,
        Entropy_scaled = NA_real_,
        ahp_index      = NA_real_
      )

    ahp_best_row       <- fit_table[1, ]
    best_model         <- final_model
    best_k             <- final_k
    recommendation_txt <- paste0(
      "User-specified Model ", final_model,
      " with k = ", final_k, " profiles."
    )

    mdata_wide <- fit_table
  }

  # Fit plot ---------------------------------------------------------------

  mdata_cols <- intersect(c("AIC", "BIC", "Entropy"), names(mdata_wide))

  if (length(mdata_cols) == 0) {
    mdata <- dplyr::tibble(
      Classes = integer(0),
      Model   = integer(0),
      name    = character(0),
      value   = numeric(0)
    )
  } else {
    mdata <- mdata_wide %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(mdata_cols),
        names_to  = "name",
        values_to = "value"
      )
  }

  model_levels <- sort(unique(c(1, 2, 3, models, final_model)))
  model_levels <- model_levels[!is.na(model_levels)]
  model_labels <- as.character(model_levels)
  model_labels[model_levels == 1] <- "1:Equal variance, cov = 0"
  model_labels[model_levels == 2] <- "2:Varying variance, cov = 0"
  model_labels[model_levels == 3] <- "3:Equal variance, equal cov"

  mdata$Model <- factor(
    mdata$Model,
    levels = model_levels,
    labels = model_labels
  )

  pal_cols <- c(
    "1:Equal variance, cov = 0"      = "#50427B",
    "2:Varying variance, cov = 0"    = "#A5C660",
    "3:Equal variance, equal cov"    = "#F16A33"
  )
  pal_cols <- pal_cols[names(pal_cols) %in% levels(mdata$Model)]

  fit_plot <- ggplot2::ggplot(
    mdata,
    ggplot2::aes(x = Classes, y = value, color = Model)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~name, scales = "free_y", ncol = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(
      x = "Number of clusters",
      y = "Fit index value"
    )

  if (length(pal_cols) > 0) {
    fit_plot <- fit_plot + ggplot2::scale_color_manual(values = pal_cols)
  }

  # Node-level cluster and posteriors --------------------------------------

  best_data <- tidyLPA::get_data(best_lpa)

  node_df <- best_data %>%
    dplyr::mutate(NodeID = dplyr::row_number()) %>%
    dplyr::rename(Cluster = Class)

  prob_cols <- grep("^CPROB", names(node_df), value = TRUE)
  for (i in seq_along(prob_cols)) {
    names(node_df)[names(node_df) == prob_cols[i]] <- paste0("prob_", i)
  }
  prob_cols_new <- grep("^prob_", names(node_df), value = TRUE)

  if (length(prob_cols_new) > 0) {
    node_df <- node_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        max_prob      = max(dplyr::c_across(dplyr::all_of(prob_cols_new)), na.rm = TRUE),
        prob_assigned = dplyr::c_across(dplyr::all_of(prob_cols_new))[Cluster],
        uncertainty   = 1 - max_prob
      ) %>%
      dplyr::ungroup()
  } else {
    node_df <- node_df %>%
      dplyr::mutate(
        max_prob      = NA_real_,
        prob_assigned = NA_real_,
        uncertainty   = NA_real_
      )
  }

  node_df <- node_df %>%
    dplyr::mutate(
      NodeID = as.integer(NodeID),
      Cluster = as.integer(Cluster)
    )

  # Map nodes to individuals -----------------------------------------------

  SOM_Node_complete <- som_model$unit.classif
  som_cluster       <- node_df$Cluster
  patient_clust     <- som_cluster[SOM_Node_complete]

  Cluster_full <- rep(NA_integer_, nrow(df_scidr))
  Cluster_full[complete_rows] <- patient_clust

  individual_tbl <- dplyr::tibble(
    .scidr_rowid = df_scidr$.scidr_rowid,
    RowID        = df_scidr$.scidr_rowid,
    SOM_Node     = SOM_Node_full,
    SOM_Distance = SOM_Dist_full
  ) %>%
    dplyr::left_join(node_df, by = c("SOM_Node" = "NodeID")) %>%
    dplyr::mutate(
      Cluster = Cluster_full
    )

  if (!is.null(id_col) && id_col %in% names(df_scidr)) {
    individual_tbl[[id_col]] <- df_scidr[[id_col]]
  }

  if (nrow(individual_tbl) != nrow(df_scidr)) {
    stop("Internal row alignment error: ProbFit$individual does not match the number of rows in df.")
  }
  if (!identical(individual_tbl$.scidr_rowid, df_scidr$.scidr_rowid)) {
    stop("Internal row alignment error: .scidr_rowid shifted during posterior probability mapping.")
  }

  # df_with_clusters: only cluster label -----------------------------------

  df_with_clusters <- df_scidr

  if (ClusterName %in% names(df_with_clusters)) {
    message("Column '", ClusterName, "' already exists and will be overwritten.")
  }

  df_with_clusters[[ClusterName]] <- individual_tbl$Cluster

  # SOMFit diagnostics and distance baselines (training) -------------------

  som_fit_tbl    <- individual_tbl
  som_fit_non_na <- som_fit_tbl[!is.na(som_fit_tbl$SOM_Distance), ]

  dist_train <- som_fit_non_na$SOM_Distance
  dist_summary <- stats::quantile(
    dist_train,
    probs = c(0, 0.25, 0.5, 0.75, 0.95, 1),
    na.rm = TRUE
  )
  names(dist_summary) <- c("min", "q25", "median", "q75", "p95", "max")

  overall_mean <- mean(dist_train, na.rm = TRUE)
  overall_sd   <- stats::sd(dist_train, na.rm = TRUE)

  dist_by_cluster <- som_fit_non_na %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      mean = mean(SOM_Distance),
      sd   = stats::sd(SOM_Distance),
      p90  = stats::quantile(SOM_Distance, 0.90),
      p95  = stats::quantile(SOM_Distance, 0.95),
      p99  = stats::quantile(SOM_Distance, 0.99),
      .groups = "drop"
    )

  overall_p95 <- dist_summary["p95"]

  flag_by_cluster <- som_fit_non_na %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      n                 = dplyr::n(),
      prop_overall_high = mean(SOM_Distance > overall_p95),
      .groups = "drop"
    )

  overall_high_cutoff <- as.numeric(stats::quantile(
    dist_train,
    probs = high_dist_quantile,
    na.rm = TRUE
  ))

  dist_by_cluster_cutoffs <- som_fit_non_na %>%
    dplyr::filter(!is.na(.data$SOM_Distance), !is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      train_cluster_high_cutoff = as.numeric(stats::quantile(
        .data$SOM_Distance,
        probs = high_dist_quantile,
        na.rm = TRUE
      )),
      .groups = "drop"
    )

  som_fit_tbl <- som_fit_tbl %>%
    dplyr::mutate(
      SOMDist_z_overall = if (is.finite(overall_sd) && overall_sd > 0) {
        (.data$SOM_Distance - overall_mean) / overall_sd
      } else {
        NA_real_
      },
      SOMDist_percentile_train = vapply(.data$SOM_Distance, function(x) {
        if (is.na(x) || length(dist_train) == 0) {
          NA_real_
        } else {
          mean(dist_train <= x, na.rm = TRUE)
        }
      }, numeric(1)),
      Flag_SOMDist_overallHigh =
        !is.na(.data$SOM_Distance) &
        !is.na(overall_high_cutoff) &
        .data$SOM_Distance > overall_high_cutoff
    ) %>%
    dplyr::left_join(dist_by_cluster_cutoffs, by = "Cluster") %>%
    dplyr::mutate(
      Flag_SOMDist_clusterHigh =
        !is.na(.data$SOM_Distance) &
        !is.na(.data$train_cluster_high_cutoff) &
        .data$SOM_Distance > .data$train_cluster_high_cutoff,
      Projection_Fit_Class = dplyr::case_when(
        is.na(.data$SOM_Distance) ~ NA_character_,
        (.data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh) &
          (is.na(.data$max_prob) | .data$max_prob < low_prob_threshold) ~
          "Potential novel phenotype",
        .data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh ~
          "Poor SOM fit",
        !is.na(.data$max_prob) & .data$max_prob < low_prob_threshold ~
          "Uncertain membership",
        TRUE ~ "Good fit"
      )
    )

  som_fit_non_na <- som_fit_tbl[!is.na(som_fit_tbl$SOM_Distance), ]

  cluster_occupancy <- som_fit_non_na %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::count(.data$Cluster, name = "n") %>%
    dplyr::mutate(prop = .data$n / sum(.data$n))

  training_cluster_fit_summary <- som_fit_non_na %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data$SOM_Distance, na.rm = TRUE),
      median_distance = stats::median(.data$SOM_Distance, na.rm = TRUE),
      mean_probability = mean(.data$prob_assigned, na.rm = TRUE),
      median_probability = stats::median(.data$prob_assigned, na.rm = TRUE),
      mean_distance_percentile = mean(.data$SOMDist_percentile_train, na.rm = TRUE),
      median_distance_percentile = stats::median(.data$SOMDist_percentile_train, na.rm = TRUE),
      prop_high_distance = mean(.data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh, na.rm = TRUE),
      prop_low_probability = mean(.data$prob_assigned < low_prob_threshold, na.rm = TRUE),
      prop_poor_fit = mean(.data$Projection_Fit_Class == "Poor SOM fit", na.rm = TRUE),
      prop_potential_novel = mean(.data$Projection_Fit_Class == "Potential novel phenotype", na.rm = TRUE),
      .groups = "drop"
    )

  potential_novel_cases <- som_fit_tbl %>%
    dplyr::filter(.data$Projection_Fit_Class == "Potential novel phenotype")

  p_dist_hist <- ggplot2::ggplot(
    som_fit_non_na,
    ggplot2::aes(x = SOM_Distance)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "SOM distances (training data)",
      x     = "Distance to BMU",
      y     = "Count"
    )

  p_dist_box <- ggplot2::ggplot(
    som_fit_non_na,
    ggplot2::aes(x = factor(Cluster), y = SOM_Distance)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "SOM distances by cluster (training data)",
      x     = "Cluster",
      y     = "Distance to BMU"
    )

  p_dist_ecdf <- ggplot2::ggplot(
    som_fit_non_na,
    ggplot2::aes(x = .data$SOM_Distance)
  ) +
    ggplot2::stat_ecdf() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Training SOM distance ECDF",
      x = "Distance to BMU",
      y = "Empirical cumulative probability"
    )

  p_cluster_fit_summary <- ggplot2::ggplot(
    training_cluster_fit_summary,
    ggplot2::aes(x = factor(.data$Cluster), y = .data$prop_high_distance)
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "High-distance training cases by cluster",
      x = "Cluster",
      y = "Proportion"
    )

  PhenotypeReference <- list(
    training_variable_summary = training_variable_summary,
    distance_reference = list(
      mean = overall_mean,
      median = stats::median(dist_train, na.rm = TRUE),
      sd = overall_sd,
      p90 = as.numeric(stats::quantile(dist_train, 0.90, na.rm = TRUE)),
      p95 = as.numeric(stats::quantile(dist_train, 0.95, na.rm = TRUE)),
      p99 = as.numeric(stats::quantile(dist_train, 0.99, na.rm = TRUE)),
      high_dist_quantile = high_dist_quantile,
      high_dist_cutoff = overall_high_cutoff
    ),
    node_occupancy = som_node_occupancy,
    cluster_occupancy = cluster_occupancy,
    projection_thresholds = list(
      high_dist_quantile = high_dist_quantile,
      low_prob_threshold = low_prob_threshold
    ),
    Health = list(
      n_training = nrow(df_scidr),
      n_complete = sum(complete_rows),
      n_excluded = sum(!complete_rows),
      prop_excluded = mean(!complete_rows),
      mean_distance = overall_mean,
      median_distance = stats::median(dist_train, na.rm = TRUE),
      node_utilization = mean(som_node_occupancy$n > 0),
      empty_nodes = sum(som_node_occupancy$n == 0),
      singleton_nodes = sum(som_node_occupancy$n == 1),
      mean_node_occupancy = mean(som_node_occupancy$n),
      median_node_occupancy = stats::median(som_node_occupancy$n)
    ),
    FitDiagnostics = list(
      overall_fit_summary = dplyr::tibble(
        n = nrow(som_fit_non_na),
        mean_distance = overall_mean,
        median_distance = stats::median(dist_train, na.rm = TRUE),
        sd_distance = overall_sd,
        high_distance_cutoff = overall_high_cutoff,
        prop_high_distance = mean(som_fit_non_na$Flag_SOMDist_overallHigh, na.rm = TRUE),
        prop_low_probability = mean(som_fit_non_na$prob_assigned < low_prob_threshold, na.rm = TRUE),
        prop_poor_fit = mean(som_fit_non_na$Projection_Fit_Class == "Poor SOM fit", na.rm = TRUE),
        prop_potential_novel = mean(som_fit_non_na$Projection_Fit_Class == "Potential novel phenotype", na.rm = TRUE)
      ),
      cluster_fit_summary = training_cluster_fit_summary,
      potential_novel_cases = potential_novel_cases
    ),
    plots = list(
      distance_ecdf = p_dist_ecdf,
      cluster_fit_summary_plot = p_cluster_fit_summary
    )
  )

  SOMFit <- list(
    train_quant_error = train_quant_err,
    distance_summary  = dist_summary,
    overall_mean      = overall_mean,
    overall_sd        = overall_sd,
    dist_by_cluster   = dist_by_cluster,
    flag_by_cluster   = flag_by_cluster,
    overall_p95       = overall_p95,
    node_occupancy    = som_node_occupancy,
    occupancy_summary = som_occupancy_summary,
    table             = som_fit_tbl,
    plots             = list(
      distance_hist           = p_dist_hist,
      distance_by_cluster_box = p_dist_box
    )
  )

  ModelInfo_SOM <- list(
    som_model = som_model,
    som_codes = som_codes,
    som_grid  = som_grid,
    training_variable_summary = training_variable_summary,
    PhenotypeReference = PhenotypeReference,
    som_grid_info = list(
      som_xdim_initial = som_grid_initial_xdim,
      som_ydim_initial = som_grid_initial_ydim,
      som_xdim_used    = som_xdim,
      som_ydim_used    = som_ydim,
      n_nodes_used     = som_xdim * som_ydim
    ),
    SOMFit    = SOMFit,
    plots     = som_plots
  )

  # ProbFit diagnostics (probabilities) ------------------------------------

  node_non_na <- node_df

  node_MaxProbBoxplot <- ggplot2::ggplot(
    node_non_na,
    ggplot2::aes(x = factor(Cluster), y = max_prob)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Node-level max posterior probability by class",
      x     = "Cluster",
      y     = "Max posterior probability"
    )

  node_ProbAssignedDensity <- ggplot2::ggplot(
    node_non_na,
    ggplot2::aes(x = prob_assigned)
  ) +
    ggplot2::geom_density(fill = "grey40", alpha = 0.7) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Density of node prob_assigned by class",
      x     = "prob_assigned",
      y     = "Density"
    )

  indiv_non_na <- som_fit_tbl[!is.na(som_fit_tbl$max_prob), ]
  indiv_non_na$Cluster <- factor(indiv_non_na$Cluster)

  individual_MaxProbBoxplot <- ggplot2::ggplot(
    indiv_non_na,
    ggplot2::aes(x = Cluster, y = max_prob)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Individual max posterior probability (training)",
      x     = "Cluster",
      y     = "Max posterior probability"
    )

  individual_ProbAssignedDensity <- ggplot2::ggplot(
    indiv_non_na,
    ggplot2::aes(x = prob_assigned)
  ) +
    ggplot2::geom_density(fill = "darkblue", alpha = 0.7) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Density of prob_assigned by class (training)",
      x     = "Posterior probability for assigned class",
      y     = "Density"
    )

  ProbFit <- list(
    node       = node_df,
    individual = som_fit_tbl,
    plots      = list(
      node_MaxProbBoxplot            = node_MaxProbBoxplot,
      node_ProbAssignedDensity       = node_ProbAssignedDensity,
      individual_MaxProbBoxplot      = individual_MaxProbBoxplot,
      individual_ProbAssignedDensity = individual_ProbAssignedDensity
    )
  )

  lpa_preprocess <- list(
    n_rows_original          = nrow(df_scidr),
    n_complete_rows          = sum(complete_rows),
    n_excluded_rows          = sum(!complete_rows),
    prop_excluded_rows       = mean(!complete_rows),
    n_som_nodes              = nrow(som_codes),
    n_lpa_variables_original = ncol(som_codes),
    n_lpa_variables_used     = ncol(X),
    dropped_zero_sd_vars     = dropped_lpa_vars,
    lpa_drop_zero_sd         = lpa_drop_zero_sd,
    lpa_zero_sd_tol          = lpa_zero_sd_tol,
    lpa_em_itmax             = lpa_em_itmax,
    lpa_em_tol               = lpa_em_tol,
    lpa_timeout_seconds      = lpa_timeout_seconds,
    skip_model_after_n_failures = skip_model_after_n_failures,
    slow_fit_seconds         = slow_fit_seconds,
    min_nodes_per_cluster    = min_nodes_per_cluster
  )

  ModelInfo_MClust <- list(
    lpa_models    = lpa_models,
    fit_table     = fit_table,
    best_fit_name = best_fit_name,
    AHP         = list(
      ahp_index      = if ("ahp_index" %in% names(fit_table)) fit_table$ahp_index else NA_real_,
      ahp_best_row   = ahp_best_row,
      recommendation = recommendation_txt
    ),
    diagnostics = list(
      lpa_fit_diagnostics = lpa_diagnostics,
      failed_fits         = failed_fits,
      warning_fits        = warning_fits,
      lpa_preprocess      = lpa_preprocess
    )
  )

  out <- list(
    method           = method,
    vars_used        = vars_used,
    ZScoreType       = ZScoreType,
    ZScoreObject     = ZScoreObject_used,
    ZScoreVars       = ZScoreVars_used,
    ClusterName      = ClusterName,
    complete_rows    = complete_rows,
    df_with_clusters = df_with_clusters,
    fit_plot         = fit_plot,
    ModelInfo_SOM    = ModelInfo_SOM,
    ModelInfo_MClust = ModelInfo_MClust,
    ProbFit          = ProbFit
  )

  class(out) <- c("Pipeline_SOMClust", class(out))
  out
}

#' SOM + latent profile clustering pipeline (with AHP and distance baselines)
#'
#' Compatibility alias for [CreateSOMClusterModel()]. Prefer
#' `CreateSOMClusterModel()` in new code because this function fits a reusable
#' SOM clustering model.
#'
#' @param ... Arguments passed to [CreateSOMClusterModel()].
#' @return The same `Pipeline_SOMClust` object returned by
#'   [CreateSOMClusterModel()].
#' @export
Pipeline_SOMClust <- function(...) {
  CreateSOMClusterModel(...)
}
