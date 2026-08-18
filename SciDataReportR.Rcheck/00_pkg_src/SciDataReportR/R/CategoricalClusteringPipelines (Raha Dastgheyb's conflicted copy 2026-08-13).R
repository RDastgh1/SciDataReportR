.DropUnusedClusterLevels <- function(data, variables) {
  for (variable in intersect(variables, names(data))) {
    if (is.factor(data[[variable]])) {
      data[[variable]] <- droplevels(data[[variable]])
    }
  }
  data
}

# Categories the frozen reduction never saw cannot be projected. Set them to NA
# so the row is handled as an incomplete case, which is how the rest of the
# package treats out-of-support projection input, instead of failing outright.
#' Fit a projectable latent class model for categorical measures
#'
#' @description Use latent class analysis when the clustering variables are
#' nominal or ordinal questionnaire, symptom, diagnosis, or assay-call items.
#' It estimates class-specific response probabilities and posterior membership.
#' @references Lazarsfeld PF, Henry NW. *Latent Structure Analysis*. Houghton
#' Mifflin; 1968. Linzer DA, Lewis JB. poLCA: An R package for polytomous
#' variable latent class analysis. *J Stat Softw.* 2011;42(10):1-29.
#' @param data Data frame containing categorical variables.
#' @param variables Variables included in the latent class model.
#' @param k_range Candidate class counts.
#' @param final_k Optional class count; when supplied, fit only this solution.
#' @param ClusterVariableName Output class column name.
#' @param nrep Number of random starts per candidate.
#' @param seed Random seed controlling latent-class random starts.
#' @param method Either `"exploratory"`/`"explore"` or `"finalize"`.
#' @param stability_resamples Number of bootstrap stability refits.
#' @param stability_seed Seed controlling bootstrap resampling.
#' @param stability_progress Whether to print bootstrap progress messages.
#' @return A fitted latent-class model with BIC, AIC, log likelihood, entropy,
#'   class-size and bootstrap stability metrics in `ModelInfo$fit_table`. Log
#'   likelihood and entropy are higher-is-better; AIC and BIC are lower-is-better.
#'   Figures sit beside what they describe: `fit_plot` reviews candidates;
#'   `ModelInfo$plots` holds `response_probabilities`, `item_profiles`,
#'   `bic`, `entropy`, `posterior_map`, and `categorical_composition`; `ProbFit$plots` holds posterior-confidence
#'   figures; and `Stability$plots` holds cluster-recovery and complementary
#'   recovery.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Categorical <- paste0("CatVar", 1:3)
#' review <- CreateClusterModel_LatentClass(
#'   df_Training, vars_Categorical, k_range = 2:4, nrep = 5,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$bic
#' model <- CreateClusterModel_LatentClass(
#'   df_Training, vars_Categorical, method = "finalize", final_k = 4,
#'   nrep = 5, stability_resamples = 2
#' )
#' model$ModelInfo$plots$response_probabilities
#' model$ModelInfo$plots$item_profiles
#' model$ProbFit$plots$confidence_density
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @export
CreateClusterModel_LatentClass <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10,
    final_k = NULL, ClusterVariableName = "Cluster", nrep = 20L, seed = 93421L,
    stability_resamples = 0L, stability_seed = seed + 1L,
    stability_progress = FALSE) {
  if (!requireNamespace("poLCA", quietly = TRUE)) stop("Package 'poLCA' is required.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (is.null(variables)) variables <- names(data)[vapply(data, function(x) is.factor(x) || is.character(x) || is.logical(x) || is.ordered(x), logical(1))]
  if (!length(variables) || length(setdiff(variables, names(data)))) stop("variables must be present categorical columns.")
  df_scidr <- data
  training <- df_scidr[variables]
  levels_used <- lapply(training, function(x) levels(factor(x)))
  training[] <- Map(function(x, lev) as.integer(factor(x, levels = lev)), training, levels_used)
  complete <- stats::complete.cases(training); if (!any(complete)) stop("No complete rows available for latent class analysis.")
  formula <- stats::as.formula(paste0("cbind(", paste(variables, collapse = ","), ") ~ 1"))
  if (!is.null(final_k) && length(method) > 1L) method <- "finalize"
  method <- .ClusterMethod(method)
  ks <- if (method == "finalize") { if (is.null(final_k)) stop("final_k is required for method = 'finalize'."); final_k } else k_range
  fits <- lapply(ks, function(k) {
    # poLCA uses random starting values. Retrying once with a fixed alternate
    # seed makes a projectable workflow robust to an occasional singular start
    # without changing the selected model when the primary fit succeeds.
    attempted <- unlist(lapply(c(seed + k, seed + 1000L + k), function(start_seed) {
      lapply(unique(c(nrep, 1L)), function(starts) {
        set.seed(start_seed)
        tryCatch(
          poLCA::poLCA(formula, training[complete, , drop = FALSE],
            nclass = k, nrep = starts, verbose = FALSE),
          error = function(e) e
        )
      })
    }), recursive = FALSE)
    successful <- Filter(function(x) !inherits(x, "error"), attempted)
    if (length(successful)) successful[[1]] else attempted[[length(attempted)]]
  })
  good <- which(!vapply(fits, inherits, logical(1), what = "error")); if (!length(good)) stop("No latent-class candidate could be estimated.")
  rows <- dplyr::bind_rows(lapply(good, function(i) {
    posterior <- fits[[i]]$posterior
    class_sizes <- tabulate(fits[[i]]$predclass, nbins = ks[[i]])
    entropy <- 1 - (-sum(pmax(posterior, .Machine$double.eps) *
      log(pmax(posterior, .Machine$double.eps))) /
      (nrow(posterior) * log(ncol(posterior))))
    dplyr::tibble(Classes = ks[[i]], BIC = fits[[i]]$bic,
      AIC = fits[[i]]$aic, LogLik = fits[[i]]$llik, Entropy = entropy,
      MinClassN = min(class_sizes),
      SizeBalance = min(class_sizes) / max(class_sizes))
  })) %>% dplyr::arrange(.data$BIC)
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(good, function(i) {
      reference <- rep(NA_integer_, nrow(df_scidr)); reference[complete] <- fits[[i]]$predclass
      .ClusterBootstrapStability(df_scidr, reference, function(boot, original) {
        fitted <- CreateClusterModel_LatentClass(boot, variables, method = "finalize",
          final_k = ks[[i]], nrep = nrep, seed = seed, stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i, list(Classes = ks[[i]]),
      progress = stability_progress, preserve_levels = variables)
    })
    stability_summary <- dplyr::bind_rows(lapply(stabilities, `[[`, "summary"))
    rows <- dplyr::left_join(rows, stability_summary, by = "Classes")
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, stability_seed)
  }
  review <- .ClusterAHP(rows, higher = c("Entropy", "MinClassN", "SizeBalance", "ReproducibilityScore"), lower = c("BIC", "AIC"), setting = "latent class count")
  rows <- review$fit_table
  recommended_k <- if (any(rows$Recommended)) rows$Classes[rows$Recommended][[1]] else rows$Classes[[1]]
  best_i <- good[[which(vapply(good, function(i) ks[[i]] == recommended_k, logical(1)))[1]]]
  best <- fits[[best_i]]
  prep <- list(data = df_scidr, complete_rows = complete, variables = variables); ind <- dplyr::tibble(PosteriorMax = rep(NA_real_, nrow(df_scidr)), Uncertainty = rep(NA_real_, nrow(df_scidr)))
  ind$PosteriorMax[complete] <- apply(best$posterior, 1, max); ind$Uncertainty[complete] <- 1 - ind$PosteriorMax[complete]
  for (i in seq_len(ncol(best$posterior))) { ind[[paste0("prob_", i)]] <- NA_real_; ind[[paste0("prob_", i)]][complete] <- best$posterior[, i] }
  response_df <- dplyr::bind_rows(lapply(names(best$probs), function(variable) {
    probabilities <- as.data.frame(as.table(best$probs[[variable]]))
    names(probabilities) <- c("Class", "Response", "Probability")
    probabilities$Response <- sub(
      "^.*: ", "", as.character(probabilities$Response))
    probabilities$Variable <- .ClusterVariableLabel(df_scidr, variable)
    probabilities
  }))
  response_df$Class <- .ClusterFactor(
    sub("^.*?([0-9]+)$", "\\1", as.character(response_df$Class)), NULL)
  p_response <- ggplot2::ggplot(response_df,
    ggplot2::aes(x = .data$Response, y = .data$Class, fill = .data$Probability)) +
    ggplot2::geom_tile() + ggplot2::facet_wrap(~Variable, scales = "free_x") +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Class-conditional response probabilities",
      subtitle = "How each latent class answers each item",
      x = "Response", y = "Latent class", fill = "Probability")
  p_item_profiles <- ggplot2::ggplot(
    response_df,
    ggplot2::aes(
      x = .data$Response, y = .data$Probability, color = .data$Class,
      group = .data$Class)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~Variable, scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    .ClusterColourScale(response_df$Class) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Latent class item profiles",
      x = "Response", y = "Response probability", color = "Latent class")
  p_posterior_map <- if (ncol(best$posterior) >= 2L) {
    posterior_plot <- as.data.frame(best$posterior)
    posterior_plot$Cluster <- best$predclass
    PlotClusterMap(
      posterior_plot, names(posterior_plot)[[1]], names(posterior_plot)[[2]],
      "Cluster", title = "Training cluster review map",
      subtitle = "Posterior membership space; corners are unambiguous cases",
      xlab = "Posterior probability: class 1",
      ylab = "Posterior probability: class 2")
  } else NULL
  ModelInfo <- list(
    lca_model = best, fit_table = rows, levels = levels_used,
    AHP = review$AHP, response_probabilities = response_df,
    final_k = ncol(best$posterior))
  base <- .ClusterOutput(prep, best$predclass, ClusterVariableName, ind,
    ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual
  # Latent class analysis has no distance metric, so classification uncertainty
  # plays the role distance plays elsewhere: it is what the frozen reference
  # thresholds and the projection triage are built on.
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, "Uncertainty", "Classification uncertainty")
  ModelInfo$plots <- .DropNullPlots(list(
    response_probabilities = p_response,
    item_profiles = p_item_profiles,
    bic = .ClusterCandidateCurve(rows, "BIC", "Bayesian information criterion"),
    entropy = .ClusterCandidateCurve(rows, "Entropy", "Classification entropy"),
    categorical_composition = PlotClusterComposition(
      base$DataWithClusters, variables, base$DataWithClusters[[ClusterVariableName]]),
    categorical_composition_by_cluster = PlotClusterComposition(
      base$DataWithClusters, variables, base$DataWithClusters[[ClusterVariableName]],
      facet_by = "cluster"),
    categorical_enrichment = PlotClusterComposition(
      base$DataWithClusters, variables, base$DataWithClusters[[ClusterVariableName]],
      style = "enrichment")))
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(
    individual, "PosteriorMax", "Posterior probability for the assigned class")
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(method = method, vars_used = variables, ClusterVariableName = ClusterVariableName, Preprocessing = list(Scaling = "Categorical", levels = levels_used), fit_plot = PlotClusterFitReview(rows)), base)
  out$Stability <- Stability
  out$ModelInfo_LatentClass <- out$ModelInfo
  out$Specification <- .ClusterSpecification(
    "LatentClass", variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, dplyr::select(rows, dplyr::any_of("Classes")),
    list(k = ncol(best$posterior)), complete,
    list(projection = "frozen class prevalence and item-response probabilities"))
  class(out) <- c("Pipeline_LatentClass", class(out)); out
}

#' Project data onto a latent class model
#' @inheritParams ProjectCluster
#' @return A full-length projection containing posterior class probabilities,
#'   uncertainty, and assignments. `ProjectionFit` triages every projected case
#'   against the frozen training uncertainty reference and holds the projected
#'   uncertainty triad, the `posterior_map`, the training-versus-projected
#'   fit-class summaries and method-specific projection diagnostics.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_LatentClass(
#'   df_Training, paste0("CatVar", 1:3), method = "finalize",
#'   final_k = 4, nrep = 5, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$plots$posterior_map
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_LatentClass <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  vars <- object$vars_used; if (length(setdiff(vars, names(new_df)))) stop("new_df is missing required clustering variables.")
  out_of_support <- dplyr::bind_rows(lapply(vars, function(variable) {
    values <- as.character(new_df[[variable]])
    bad <- !is.na(values) & !values %in% object$Preprocessing$levels[[variable]]
    if (!any(bad)) return(dplyr::tibble())
    dplyr::tibble(RowID = which(bad), Variable = variable, Value = values[bad],
      Issue = "category absent from training")
  }))
  prep <- list(data = new_df, complete_rows = stats::complete.cases(new_df[vars]), variables = vars)
  coded <- prep$data[vars]; coded[] <- Map(function(x, lev) as.integer(factor(x, levels = lev)), coded, object$Preprocessing$levels)
  valid <- prep$complete_rows & stats::complete.cases(coded); fit <- object$ModelInfo$lca_model; posterior <- matrix(NA_real_, nrow(coded), length(fit$P))
  if (any(valid)) for (row_i in which(valid)) {
    likelihood <- fit$P
    for (v in vars) likelihood <- likelihood * fit$probs[[v]][, coded[[v]][row_i]]
    posterior[row_i, ] <- likelihood / sum(likelihood)
  }
  cl <- rep(NA_integer_, nrow(prep$data)); cl[valid] <- max.col(posterior[valid, , drop = FALSE])
  ind <- dplyr::tibble(
    PosteriorMax = rep(NA_real_, nrow(prep$data)),
    Uncertainty = rep(NA_real_, nrow(prep$data)))
  if (any(valid)) ind$PosteriorMax[valid] <- apply(
    posterior[valid, , drop = FALSE], 1, max)
  ind$Uncertainty[valid] <- 1 - ind$PosteriorMax[valid]
  for (i in seq_len(ncol(posterior))) ind[[paste0("prob_", i)]] <- posterior[, i]
  p_posterior_map <- if (ncol(posterior) >= 2L) {
    posterior_plot <- as.data.frame(posterior)
    posterior_plot$Cluster <- cl
    PlotClusterMap(
      posterior_plot, names(posterior_plot)[[1]], names(posterior_plot)[[2]],
      "Cluster", title = "Projected cluster review map",
      subtitle = "Frozen posterior membership space",
      xlab = "Posterior probability: class 1",
      ylab = "Posterior probability: class 2")
  } else NULL
  base <- .ClusterOutput(
    prep, cl[prep$complete_rows], ClusterVariableName, ind, object$ModelInfo)
  individual <- base$ProbFit$individual
  ProjectionFit <- .ClusterProjectionFit(
    individual, object$ModelInfo$FitDiagnostics$distance_reference,
    "Uncertainty", "Classification uncertainty",
    training_cluster = object$ProbFit$individual$Cluster)
  ProjectionFit$out_of_support <- out_of_support
  ProjectionFit$plots$posterior_map <- p_posterior_map
  out <- c(list(vars_used = vars, ClusterVariableName = ClusterVariableName,
    Preprocessing = object$Preprocessing), base)
  out$ProbFit <- .ClusterProbFit(
    individual, "PosteriorMax", "Posterior probability for the assigned class",
    prefix = "Projected")
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- ProjectionFit$fit_class
  out$ProbFit$individual$Projection_Fit_Class <- ProjectionFit$fit_class
  class(out) <- c("Project_LatentClass", class(out)); out
}

#' Fit MCA followed by Mclust for nominal categorical data
#'
#' @description Use this pipeline for many nominal categorical variables when a
#' lower-dimensional category space is useful before mixture clustering.
#' @references Greenacre M. *Correspondence Analysis in Practice*. 3rd ed.
#' Chapman and Hall/CRC; 2017. Scrucca L et al. mclust 5. *J Stat Softw.* 2016;71(11):1-29.
#' @inheritParams CreateClusterModel_MClust
#' @param mca_variance_threshold Cumulative MCA inertia percentage retained.
#' @return A fitted MCA plus Mclust model containing the frozen MCA object,
#'   mixture-model fit and bootstrap stability metrics, and posterior
#'   probabilities. `ModelInfo$plots` leads with the reduction layer's `scree`
#'   and `loadings`, followed by the mixture-model structure figures in MCA
#'   score space and `categorical_composition` restated on the original items.
#'   BIC/ICL/AIC, entropy, and uncertainty use the same interpretation as Mclust.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Categorical <- paste0("CatVar", 1:3)
#' review <- CreateClusterModel_MCA_MClust(
#'   df_Training, vars_Categorical, k_range = 2:4, models = 1,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$scree
#' model <- CreateClusterModel_MCA_MClust(
#'   df_Training, vars_Categorical, method = "finalize", final_k = 4,
#'   final_model = 1, stability_resamples = 2
#' )
#' model$ModelInfo$plots$separation_map
#' model$ModelInfo$plots$categorical_composition
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @export
CreateClusterModel_MCA_MClust <- function(data, variables,
    method = c("exploratory", "finalize"), k_range = 2:10,
    models = c(1L, 2L, 3L), final_k = NULL,
    final_model = NULL, ClusterVariableName = "Cluster", seed = 93421L,
    mca_variance_threshold = 75, stability_resamples = 0L,
    stability_seed = seed + 1L, stability_progress = FALSE) {
  if (!requireNamespace("FactoMineR", quietly = TRUE)) stop("Package 'FactoMineR' is required.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  if (!is.null(final_k) && length(method) > 1L) method <- "finalize"
  method <- .ClusterMethod(method)
  models <- .ResolveMclustModels(models)
  if (method == "finalize") final_model <- .ResolveMclustModels(final_model, "final_model")
  # FactoMineR validates projection levels against the categories it actually
  # saw, so a declared-but-unobserved level would make every later projection
  # fail. Freeze the MCA on observed categories only.
  data <- .DropUnusedClusterLevels(data, variables)
  mca <- CreateMCAObject(
    data = data, VarsToReduce = variables,
    minThresh = mca_variance_threshold, Relabel = FALSE)
  scores <- as.data.frame(mca$Scores$coord)
  score_vars <- names(scores)
  review_model <- CreateClusterModel_MClust(
    scores, score_vars, method = method, k_range = k_range,
    models = models, final_k = final_k, final_model = final_model,
    ZScoreType = "None", ClusterVariableName = ClusterVariableName, seed = seed,
    stability_resamples = 0L)
  fit_table <- review_model$ModelInfo$fit_table
  Stability <- NULL
  if (stability_resamples > 0L) {
    candidate_rows <- dplyr::distinct(
      fit_table, .data$Model, .data$Classes)
    stabilities <- lapply(seq_len(nrow(candidate_rows)), function(i) {
      candidate_k <- candidate_rows$Classes[[i]]
      candidate_model <- candidate_rows$Model[[i]]
      reference_fit <- CreateClusterModel_MCA_MClust(
        data, variables, method = "finalize", k_range = candidate_k,
        models = candidate_model, final_k = candidate_k,
        final_model = candidate_model, ClusterVariableName = ClusterVariableName, seed = seed,
        mca_variance_threshold = mca_variance_threshold,
        stability_resamples = 0L)
      reference <- reference_fit$ProbFit$individual$Cluster
      .ClusterBootstrapStability(data, reference, function(boot, original) {
        fitted <- CreateClusterModel_MCA_MClust(
          boot, variables, method = "finalize", k_range = candidate_k,
          models = candidate_model, final_k = candidate_k,
          final_model = candidate_model, ClusterVariableName = ClusterVariableName,
          seed = seed, mca_variance_threshold = mca_variance_threshold,
          stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i,
      list(Model = candidate_model, Classes = candidate_k),
      progress = stability_progress, preserve_levels = variables)
    })
    Stability <- .CombineClusterStabilities(
      stabilities, stability_resamples, stability_seed)
    fit_table <- dplyr::left_join(
      fit_table, Stability$summary, by = c("Model", "Classes"))
  }
  review <- .ClusterAHP(
    fit_table,
    higher = c("BIC", "ICL", "Entropy", "MinClusterN", "SizeBalance", "ReproducibilityScore"),
    lower = "MaxUncertainty", setting = "MCA plus Gaussian-mixture model")
  selected_k <- review$AHP$ahp_best_row$Classes[[1]]
  selected_model <- review$AHP$ahp_best_row$Model[[1]]
  model <- CreateClusterModel_MClust(
    scores, score_vars, method = "finalize", final_k = selected_k,
    final_model = selected_model, ZScoreType = "None", ClusterVariableName = ClusterVariableName,
    seed = seed, stability_resamples = 0L)
  model$method <- method
  model$ModelInfo$fit_table <- review$fit_table
  model$ModelInfo$AHP <- review$AHP
  model$Stability <- Stability
  model$fit_plot <- PlotClusterFitReview(review$fit_table)
  df_scidr <- data
  model$DataWithClusters <- df_scidr
  model$DataWithClusters[[ClusterVariableName]] <- model$ProbFit$individual$Cluster
  # The inner mixture model built its structure figures in MCA score space,
  # which is where this pipeline clusters. Add the reduction layer's own
  # figures and restate composition on the original items, because that is what
  # a reader has to interpret.
  score_plots <- model$ModelInfo$plots[
    !names(model$ModelInfo$plots) %in% "profiles"]
  model$ModelInfo$plots <- .DropNullPlots(c(
    list(scree = mca$p_scree, loadings = mca$Lollipop),
    score_plots,
    list(
      categorical_composition = PlotClusterComposition(
        df_scidr, variables, model$ProbFit$individual$Cluster),
      categorical_composition_by_cluster = PlotClusterComposition(
        df_scidr, variables, model$ProbFit$individual$Cluster,
        facet_by = "cluster"),
      categorical_enrichment = PlotClusterComposition(
        df_scidr, variables, model$ProbFit$individual$Cluster,
        style = "enrichment"))))
  if (!is.null(Stability)) {
    model$Stability$plots <- .ClusterStabilityPlots(Stability)
  }
  model$vars_used <- variables
  model$Preprocessing <- list(
    Scaling = "MCA", MCAObject = mca, MCAVariables = score_vars)
  model$MCA <- mca
  model$ModelInfo_MCA <- list(MCAObject = mca, MCAVariables = score_vars)
  model$ModelInfo_MClust <- model$ModelInfo
  model$Specification <- .ClusterSpecification(
    "MCA + Mclust", variables, list(seed = seed, stability_seed = stability_seed),
    model$Preprocessing, model$ModelInfo$fit_table, list(),
    stats::complete.cases(data[variables]))
  class(model) <- c("Pipeline_MCA_MClust", class(model))
  model
}

#' Project data onto an MCA + Mclust model
#' @inheritParams ProjectCluster
#' @return A full-length projection through the frozen MCA and Mclust layers,
#'   with posterior membership in `ProbFit` and the projection triage in
#'   `ProjectionFit`. Incomplete categorical rows are retained with `NA`
#'   assignments.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_MCA_MClust(
#'   df_Training, paste0("CatVar", 1:3), method = "finalize",
#'   final_k = 4, final_model = 1, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' projected$ProjectionFit$plots$training_vs_projected_cluster_occupancy
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_MCA_MClust <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  vars <- object$vars_used
  if (length(setdiff(vars, names(new_df)))) {
    stop("new_df is missing required clustering variables.")
  }
  df_scidr <- new_df
  reference_levels <- lapply(
    object$Preprocessing$MCAObject$mcaresults$call$X[vars], levels)
  for (variable in vars) {
    values <- as.character(df_scidr[[variable]])
    values[!is.na(values) & !values %in% reference_levels[[variable]]] <- NA_character_
    df_scidr[[variable]] <- factor(values, levels = reference_levels[[variable]])
  }
  complete <- stats::complete.cases(df_scidr[vars])
  prep <- list(data = df_scidr, complete_rows = complete, variables = vars)
  if (!any(complete)) {
    individual <- dplyr::tibble(
      PosteriorMax = rep(NA_real_, nrow(df_scidr)),
      Uncertainty = rep(NA_real_, nrow(df_scidr)))
    result <- .ClusterOutput(
      prep, integer(), ClusterVariableName, individual, object$ModelInfo)
    ProjectionFit <- NULL
    ProbFit <- result$ProbFit
  } else {
    scores <- as.data.frame(FactoMineR::predict.MCA(
      object$Preprocessing$MCAObject$mcaresults,
      newdata = df_scidr[complete, vars, drop = FALSE])$coord)
    scores <- scores[, object$Preprocessing$MCAVariables, drop = FALSE]
    score_model <- object
    score_model$vars_used <- names(scores)
    score_model$Preprocessing <- list(ZScoreType = "None", Scaling = "None")
    class(score_model) <- "Pipeline_MClust"
    projected_scores <- ProjectCluster(
      score_model, scores, ClusterVariableName)
    reduced <- projected_scores$ProbFit$individual
    indicator_names <- setdiff(names(reduced), "Cluster")
    individual <- .PadIndividual(
      reduced, indicator_names, nrow(df_scidr), complete)
    result <- .ClusterOutput(
      prep, reduced$Cluster, ClusterVariableName, individual, object$ModelInfo)
    ProbFit <- projected_scores$ProbFit
    ProbFit$individual <- result$ProbFit$individual
    ProjectionFit <- projected_scores$ProjectionFit
    ProjectionFit$plots$categorical_composition <- PlotClusterComposition(
      df_scidr, vars, result$ProbFit$individual$Cluster)
  }
  projected <- c(list(
    vars_used = vars, ClusterVariableName = ClusterVariableName,
    Preprocessing = object$Preprocessing), result)
  projected$ProbFit <- ProbFit
  projected$ProjectionFit <- ProjectionFit
  projected$DataWithClusters$Projection_Fit_Class <-
    projected$ProbFit$individual$Projection_Fit_Class
  class(projected) <- c("Project_MCAMclust", class(projected))
  projected
}
