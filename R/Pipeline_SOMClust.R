#' SOM + latent profile clustering pipeline (with AHP and distance baselines)
#'
#' @description
#' End-to-end pipeline to:
#' - Standardize variables using SciDataReportR::CalcZScore() or a supplied
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
#'
#' Missing data:
#' - SOM and clustering are fit only on rows with complete Z-scores.
#' - The returned \code{df_with_clusters} has the full original data with
#'   a single cluster column appended; rows not used in SOM/LPA get NA.
#'
#' @param df Data frame containing the variables to be used in SOM and
#'   clustering.
#' @param variables Optional character vector of variable names. If NULL,
#'   numeric variables are auto-detected using
#'   \code{SciDataReportR::getNumVars(df, Ordinal = FALSE)}.
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
#'   }
#' @param ZScoreObject Optional ZScoreObj (from \code{CalcZScore()} or
#'   \code{Project_ZScore()}) to use when \code{ZScoreType = "ZScoreObj"}.
#' @param som_xdim,som_ydim Optional integers for SOM grid dimensions. If NULL,
#'   a square grid with side length \code{ceiling(n_complete^(1/3))} is used.
#' @param som_topo SOM topology for \code{kohonen::somgrid()}, default
#'   \code{"hexagonal"}.
#' @param som_neigh SOM neighbourhood function, default \code{"gaussian"}.
#' @param seed_som,seed_lpa Integer seeds for SOM and LPA steps
#'   (defaults 934521 and 93421).
#' @param Relabel Logical; if TRUE (default), aweSOM plots are relabeled
#'   using variable labels from \code{df} (via Hmisc or sjlabelled when
#'   available).
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
#' @return A list of class \code{"Pipeline_SOMClust"} with components:
#'   \itemize{
#'     \item \code{method}, \code{vars_used}, \code{ZScoreType},
#'           \code{ZScoreObject}, \code{ClusterName}
#'     \item \code{complete_rows}: logical vector (rows used for SOM/LPA)
#'     \item \code{df_with_clusters}: original \code{df} with only the
#'           cluster column appended
#'     \item \code{fit_plot}: ggplot of AIC/BIC/Entropy vs k and model
#'     \item \code{ModelInfo_SOM}: list with \code{som_model},
#'           \code{som_codes}, \code{som_grid}, \code{SOMFit} (distance
#'           diagnostics, baselines, and per-cluster flags), \code{plots}
#'           (aweSOM plots)
#'     \item \code{ModelInfo_MClust}: list with \code{lpa_models},
#'           \code{fit_table}, and \code{AHP} information
#'     \item \code{ProbFit}: list with \code{node} (node-level posterior
#'           probabilities), \code{individual} (per-person mapping and
#'           probabilities), and probability plots
#'   }
#'
#' @references
#' Saaty TL. \emph{The Analytic Hierarchy Process}. McGraw–Hill, 1980.
#'
#' @export
Pipeline_SOMClust <- function(
    df,
    variables    = NULL,
    method       = c("exploratory", "finalize"),
    k_range      = 2:15,
    models       = c(1, 2, 3),
    final_k      = NULL,
    final_model  = NULL,
    ClusterName  = "Cluster",
    ZScoreType   = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj"),
    ZScoreObject = NULL,
    som_xdim     = NULL,
    som_ydim     = NULL,
    som_topo     = "hexagonal",
    som_neigh    = "gaussian",
    seed_som     = 934521L,
    seed_lpa     = 93421L,
    Relabel      = TRUE
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
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr', 'tidyr', and 'ggplot2' are required.")
  }

  # Variables --------------------------------------------------------------

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

  # Z-scores ---------------------------------------------------------------

  if (ZScoreType == "ZScoreObj") {
    if (is.null(ZScoreObject) || !"ZScoreObj" %in% class(ZScoreObject)) {
      stop("ZScoreType = 'ZScoreObj' requires a valid ZScoreObj.")
    }
    z_res <- SciDataReportR::Project_ZScore(
      df                 = df,
      variables          = variables,
      parameters         = ZScoreObject,
      ParameterInputType = "ZScoreObj",
      names_prefix       = "Z_"
    )
    ZScoreObject_used <- z_res
  } else {
    center_flag <- ZScoreType %in% c("Center and Scale", "Center Only")
    scale_flag  <- ZScoreType %in% c("Center and Scale", "Scale Only")

    ZScoreObject_used <- SciDataReportR::CalcZScore(
      df           = df,
      variables    = variables,
      names_prefix = "Z_",
      center       = center_flag,
      scale        = scale_flag
    )
  }

  z_df <- ZScoreObject_used$ZScores
  complete_rows <- stats::complete.cases(z_df)
  if (!any(complete_rows)) {
    stop("No complete rows after Z-score step; cannot fit SOM.")
  }

  zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])

  # SOM fitting ------------------------------------------------------------

  n_complete <- nrow(zmat)
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

  SOM_Node_full <- rep(NA_integer_, nrow(df))
  SOM_Node_full[complete_rows] <- som_model$unit.classif

  SOM_Dist_full <- rep(NA_real_, nrow(df))
  SOM_Dist_full[complete_rows] <- som_model$distances

  train_quant_err <- mean(som_model$distances, na.rm = TRUE)

  # AweSOM plots -----------------------------------------------------------

  CircPlot  <- aweSOM::aweSOMplot(som = som_model, type = "Circular", palvar = "viridis")
  LinePlot  <- aweSOM::aweSOMplot(som = som_model, type = "Line",     palvar = "viridis")
  CloudPlot <- aweSOM::aweSOMplot(som = som_model, type = "Cloud")

  if (Relabel) {
    relabel_fun <- function(w) {
      if (!is.null(w$x$label)) {
        vars <- w$x$label
        get_lab <- NULL
        if (requireNamespace("Hmisc", quietly = TRUE)) {
          get_lab <- function(v) Hmisc::label(df[[v]])
        } else if (requireNamespace("sjlabelled", quietly = TRUE)) {
          get_lab <- function(v) sjlabelled::get_label(df[[v]])
        }
        if (!is.null(get_lab)) {
          new_labels <- vapply(vars, function(v) {
            if (!v %in% names(df)) return(v)
            lab <- get_lab(v)
            if (is.null(lab) || !nzchar(lab)) v else as.character(lab)
          }, character(1))
          w$x$label <- unname(new_labels)
        }
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
  set.seed(seed_lpa)

  if (method == "exploratory") {
    lpa_models <- suppressWarnings(
      suppressMessages(
        tidyLPA::estimate_profiles(
          X,
          n_profiles = k_range,
          models     = models
        )
      )
    )

    comp <- tidyLPA::compare_solutions(
      lpa_models,
      statistics = c("AIC", "BIC", "Entropy")
    )

    fit_table <- comp$fit

    fit_table <- fit_table %>%
      dplyr::mutate(
        AIC_scaled     = as.numeric(scale(-AIC)),
        BIC_scaled     = as.numeric(scale(-BIC)),
        Entropy_scaled = as.numeric(scale(Entropy)),
        ahp_index      = (AIC_scaled + BIC_scaled + Entropy_scaled) / 3
      )

    best_idx     <- which.max(fit_table$ahp_index)
    ahp_best_row <- fit_table[best_idx, ]
    best_model   <- as.integer(as.character(ahp_best_row$Model))
    best_k       <- as.integer(ahp_best_row$Classes)

    recommendation_txt <- paste0(
      "AHP (AIC, BIC, Entropy) recommends Model ",
      best_model, " with k = ", best_k, " profiles."
    )

    best_lpa <- suppressWarnings(
      suppressMessages(
        tidyLPA::estimate_profiles(
          X,
          n_profiles = best_k,
          models     = best_model
        )
      )
    )

    mdata_wide <- as.data.frame(comp$fits)

  } else {
    if (is.null(final_k) || is.null(final_model)) {
      stop("For method = 'finalize', final_k and final_model must be supplied.")
    }

    lpa_models <- suppressWarnings(
      suppressMessages(
        tidyLPA::estimate_profiles(
          X,
          n_profiles = final_k,
          models     = final_model
        )
      )
    )

    best_lpa <- lpa_models

    fit <- tidyLPA::get_fit(best_lpa)
    fit_table <- fit %>%
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

    mdata_wide <- data.frame(
      Classes = final_k,
      Model   = final_model,
      AIC     = fit$AIC,
      BIC     = fit$BIC,
      Entropy = fit$Entropy
    )
  }

  # Fit plot ---------------------------------------------------------------

  mdata <- mdata_wide %>%
    tidyr::pivot_longer(
      cols      = c("AIC", "BIC", "Entropy"),
      names_to  = "name",
      values_to = "value"
    )

  mdata$Model <- factor(
    mdata$Model,
    levels = c(1, 2, 3),
    labels = c(
      "1:Equal variance, cov = 0",
      "2:Varying variance, cov = 0",
      "3:Equal variance, equal cov"
    )
  )

  pal_cols <- c(
    "1:Equal variance, cov = 0"      = "#50427B",
    "2:Varying variance, cov = 0"    = "#A5C660",
    "3:Equal variance, equal cov"    = "#F16A33"
  )

  fit_plot <- ggplot2::ggplot(
    mdata,
    ggplot2::aes(x = Classes, y = value, color = Model)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~name, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = pal_cols) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(
      x = "Number of clusters",
      y = "Fit index value"
    )

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

  node_df <- node_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_prob      = max(dplyr::c_across(dplyr::all_of(prob_cols_new))),
      prob_assigned = dplyr::c_across(dplyr::all_of(prob_cols_new))[Cluster],
      uncertainty   = 1 - max_prob
    ) %>%
    dplyr::ungroup()

  # Map nodes to individuals -----------------------------------------------

  SOM_Node_complete <- som_model$unit.classif
  som_cluster       <- node_df$Cluster
  patient_clust     <- som_cluster[SOM_Node_complete]

  Cluster_full <- rep(NA_integer_, nrow(df))
  Cluster_full[complete_rows] <- patient_clust

  individual_tbl <- dplyr::tibble(
    RowID        = seq_len(nrow(df)),
    SOM_Node     = SOM_Node_full,
    SOM_Distance = SOM_Dist_full
  ) %>%
    dplyr::left_join(node_df, by = c("SOM_Node" = "NodeID"))

  # df_with_clusters: only cluster label -----------------------------------

  df_with_clusters <- df

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

  SOMFit <- list(
    train_quant_error = train_quant_err,
    distance_summary  = dist_summary,
    overall_mean      = overall_mean,
    overall_sd        = overall_sd,
    dist_by_cluster   = dist_by_cluster,
    flag_by_cluster   = flag_by_cluster,
    overall_p95       = overall_p95,
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

  indiv_non_na <- individual_tbl[!is.na(individual_tbl$max_prob), ]
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
    individual = individual_tbl,
    plots      = list(
      node_MaxProbBoxplot            = node_MaxProbBoxplot,
      node_ProbAssignedDensity       = node_ProbAssignedDensity,
      individual_MaxProbBoxplot      = individual_MaxProbBoxplot,
      individual_ProbAssignedDensity = individual_ProbAssignedDensity
    )
  )

  ModelInfo_MClust <- list(
    lpa_models = lpa_models,
    fit_table  = fit_table,
    AHP        = list(
      ahp_index      = if ("ahp_index" %in% names(fit_table)) fit_table$ahp_index else NA_real_,
      ahp_best_row   = ahp_best_row,
      recommendation = recommendation_txt
    )
  )

  out <- list(
    method           = method,
    vars_used        = variables,
    ZScoreType       = ZScoreType,
    ZScoreObject     = ZScoreObject_used,
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
