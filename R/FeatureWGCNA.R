#' Build a weighted feature correlation network
#'
#' @param data A data frame with samples in rows and features in columns.
#' @param sample_id Optional unique sample-ID column retained for trait joins.
#' @param variables Explicit feature columns to include.
#' @param network_type One of `"signed"`, `"unsigned"`, or `"signed_hybrid"`.
#' @param soft_power Optional soft-thresholding power.
#' @param sft_rsq Target scale-free topology fit for automatic power selection.
#' @param min_module_size Minimum dynamic-tree-cut module size.
#' @param deep_split Dynamic tree-cut sensitivity from 0 through 4.
#' @param merge_cut_height Eigengene dissimilarity for module merging.
#' @param keep_tom Store the topological-overlap matrix.
#' @param seed Optional random seed.
#'
#' @return A `feature_wgcna_obj` with module assignments, eigengenes, diagnostics,
#'   hub statistics, and the matrix used for fitting.
#' @export
BuildFeatureWGCNA <- function(
    data,
    variables,
    sample_id = NULL,
    network_type = c("signed", "unsigned", "signed_hybrid"),
    soft_power = NULL,
    sft_rsq = 0.85,
    min_module_size = 30,
    deep_split = 2,
    merge_cut_height = 0.25,
    keep_tom = FALSE,
    seed = NULL
) {
  network_type <- match.arg(network_type)
  if (!requireNamespace("WGCNA", quietly = TRUE) || !requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("`BuildFeatureWGCNA()` requires optional packages `WGCNA` and `dynamicTreeCut`.", call. = FALSE)
  }
  if (!is.data.frame(data) || missing(variables) || is.null(variables)) {
    stop("`data` must be a data frame and `variables` must explicitly name feature columns.", call. = FALSE)
  }
  variables <- as.character(variables)
  if (anyNA(variables) || any(!nzchar(variables)) || anyDuplicated(variables)) {
    stop("`variables` must contain unique, non-missing feature names.", call. = FALSE)
  }
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0 || length(variables) < 3) {
    stop("`variables` must contain at least three columns present in `data`. Missing: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(sft_rsq) || length(sft_rsq) != 1 || sft_rsq <= 0 || sft_rsq > 1 ||
      !is.numeric(min_module_size) || length(min_module_size) != 1 || min_module_size < 2 ||
      !is.numeric(deep_split) || length(deep_split) != 1 || !is.finite(deep_split) || deep_split != as.integer(deep_split) || !deep_split %in% 0:4 ||
      !is.numeric(merge_cut_height) || length(merge_cut_height) != 1 || merge_cut_height <= 0 || merge_cut_height >= 1) {
    stop("`sft_rsq`, `min_module_size`, `deep_split`, and `merge_cut_height` must be valid scalar WGCNA parameters.", call. = FALSE)
  }
  if (!is.null(soft_power) && (!is.numeric(soft_power) || length(soft_power) != 1 || !is.finite(soft_power) || soft_power <= 0)) {
    stop("`soft_power` must be NULL or a positive numeric value.", call. = FALSE)
  }
  dat_expr <- as.matrix(data[, variables, drop = FALSE])
  if (!is.null(sample_id)) {
    if (!is.character(sample_id) || length(sample_id) != 1 || !sample_id %in% names(data) || anyNA(data[[sample_id]]) || anyDuplicated(data[[sample_id]])) {
      stop("`sample_id` must name a unique, non-missing column in `data`.", call. = FALSE)
    }
    if (sample_id %in% variables) stop("`sample_id` cannot also be a WGCNA feature.", call. = FALSE)
    attr(dat_expr, "sample_ids") <- as.character(data[[sample_id]])
  }
  storage.mode(dat_expr) <- "double"
  if (anyNA(dat_expr)) stop("WGCNA requires complete data in selected features; impute or remove missing values before fitting.", call. = FALSE)
  if (any(!is.finite(dat_expr))) stop("Selected feature values must be finite.", call. = FALSE)
  constant <- variables[vapply(seq_along(variables), function(i) stats::sd(dat_expr[, i]) == 0, logical(1))]
  if (length(constant) > 0) stop("WGCNA cannot use constant features: ", paste(constant, collapse = ", "), call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  wgcna_type <- c(signed = "signed", unsigned = "unsigned", signed_hybrid = "signed hybrid")[[network_type]]
  tom_type <- if (network_type == "unsigned") "unsigned" else "signed"
  invisible(utils::capture.output(suppressMessages(sft <- WGCNA::pickSoftThreshold(dat_expr, powerVector = 1:20, networkType = wgcna_type, RsquaredCut = sft_rsq, verbose = 0))))
  auto_selected <- is.null(soft_power)
  chosen_power <- if (auto_selected) sft$powerEstimate else soft_power
  if (is.na(chosen_power)) {
    chosen_power <- if (network_type == "unsigned") 6L else 12L
    warning("No soft-thresholding power reached the requested fit; using ", chosen_power, ". Review `PlotFeatureWGCNASoftThreshold()`.", call. = FALSE)
  }
  invisible(utils::capture.output(suppressMessages({
    adjacency <- WGCNA::adjacency(dat_expr, power = chosen_power, type = wgcna_type)
    tom <- WGCNA::TOMsimilarity(adjacency, TOMType = tom_type, verbose = 0)
    dissimilarity <- 1 - tom
    feature_tree <- stats::hclust(stats::as.dist(dissimilarity), method = "average")
    dynamic <- dynamicTreeCut::cutreeDynamic(feature_tree, distM = dissimilarity, deepSplit = deep_split, pamRespectsDendro = FALSE, minClusterSize = min_module_size)
    merged <- WGCNA::mergeCloseModules(dat_expr, WGCNA::labels2colors(dynamic), cutHeight = merge_cut_height, verbose = 0)
    module_colors <- merged$colors
    eigengenes <- WGCNA::moduleEigengenes(dat_expr, module_colors)$eigengenes
    eigengenes <- eigengenes[, order(names(eigengenes)), drop = FALSE]
  })))
  module_membership <- stats::cor(dat_expr, eigengenes, use = "pairwise.complete.obs")
  own_me <- paste0("ME", module_colors)
  kme <- module_membership[cbind(seq_along(variables), match(own_me, colnames(module_membership)))]
  modules <- tibble::tibble(Variable = variables, Module = module_colors, ModuleMembership = as.numeric(kme))
  structure(list(modules = modules, module_sizes = table(module_colors), eigengenes = eigengenes,
                 soft_threshold = list(power = chosen_power, auto_selected = auto_selected, sft_rsq_target = sft_rsq, power_estimate = sft$powerEstimate, fit_indices = sft$fitIndices),
                 network_type = network_type, params = list(min_module_size = min_module_size, deep_split = deep_split, merge_cut_height = merge_cut_height, seed = seed),
                 datExpr = dat_expr, gene_tree = feature_tree, module_colors = module_colors,
                 tom = if (isTRUE(keep_tom)) tom else NULL), class = "feature_wgcna_obj")
}

#' Analyze feature module relationships with sample traits
#'
#' @param wgcna_obj A `feature_wgcna_obj`.
#' @param trait_data Trait data with a unique sample identifier.
#' @param sample_id Sample-ID column shared with the input data.
#' @param trait_sample_id Sample-ID column in `trait_data`.
#' @param outcome Trait columns to analyze.
#' @param covariates Optional covariate columns.
#'
#' @return A `feature_module_trait_obj` with association results and join audit.
#' @export
AnalyzeFeatureModuleTraits <- function(wgcna_obj, trait_data, sample_id, trait_sample_id = sample_id, outcome = NULL, covariates = NULL) {
  if (!inherits(wgcna_obj, "feature_wgcna_obj") || !is.data.frame(trait_data)) stop("Supply a `feature_wgcna_obj` and data-frame `trait_data`.", call. = FALSE)
  sample_ids <- attr(wgcna_obj$datExpr, "sample_ids")
  if (is.null(sample_ids)) stop("This WGCNA object has no sample IDs. Rebuild it with `sample_id` support.", call. = FALSE)
  if (!trait_sample_id %in% names(trait_data) || anyNA(trait_data[[trait_sample_id]]) || anyDuplicated(trait_data[[trait_sample_id]])) stop("`trait_data` must contain unique, non-missing trait sample IDs.", call. = FALSE)
  audit <- tibble::tibble(SampleID = sample_ids, Matched = sample_ids %in% trait_data[[trait_sample_id]])
  if (any(!audit$Matched) || any(!trait_data[[trait_sample_id]] %in% sample_ids)) stop("Sample-ID join is incomplete. Matrix-only: ", sum(!audit$Matched), "; trait-only: ", sum(!trait_data[[trait_sample_id]] %in% sample_ids), ".", call. = FALSE)
  traits <- trait_data[match(sample_ids, trait_data[[trait_sample_id]]), , drop = FALSE]
  if (is.null(outcome)) outcome <- setdiff(names(traits), c(trait_sample_id, covariates))
  if (length(setdiff(c(outcome, covariates), names(traits))) > 0) stop("Requested outcome or covariate columns are absent from `trait_data`.", call. = FALSE)
  mes <- wgcna_obj$eigengenes[, colnames(wgcna_obj$eigengenes) != "MEgrey", drop = FALSE]
  rows <- list(); index <- 1L
  for (trait_name in outcome) for (module_name in colnames(mes)) {
    df_model <- data.frame(Eigengene = mes[, module_name], Trait = traits[[trait_name]], traits[, covariates, drop = FALSE], check.names = FALSE)
    df_model <- df_model[stats::complete.cases(df_model), , drop = FALSE]
    observed <- unique(df_model$Trait)
    family <- if (length(observed) == 2) "binary" else if (is.numeric(df_model$Trait)) "continuous" else if (length(observed) > 2) "multicategory" else "unsupported"
    result <- tibble::tibble(Module = module_name, Trait = trait_name, OutcomeFamily = family, Effect = NA_real_, EffectType = NA_character_, PValue = NA_real_, N = nrow(df_model), Note = NA_character_)
    if (family == "continuous" && nrow(df_model) > 2) { fit <- stats::lm(stats::reformulate(c("Trait", covariates), response = "Eigengene"), data = df_model); result$Effect <- stats::coef(fit)["Trait"]; result$EffectType <- "linear beta"; result$PValue <- summary(fit)$coefficients["Trait", "Pr(>|t|)"] }
    if (family == "binary" && nrow(df_model) > 3) { df_model$Trait <- factor(df_model$Trait); fit <- stats::glm(stats::reformulate(c("Eigengene", covariates), response = "Trait"), data = df_model, family = stats::binomial()); term <- grep("Eigengene", rownames(summary(fit)$coefficients), value = TRUE)[1]; result$Effect <- stats::coef(fit)[term]; result$EffectType <- "log odds ratio"; result$PValue <- summary(fit)$coefficients[term, "Pr(>|z|)"] }
    if (family == "multicategory" && nrow(df_model) > 3) { fit <- stats::lm(stats::reformulate(c("Trait", covariates), response = "Eigengene"), data = df_model); result$EffectType <- "omnibus F"; result$PValue <- stats::anova(fit)["Trait", "Pr(>F)"]; result$Effect <- stats::anova(fit)["Trait", "F value"] }
    if (family == "unsupported") result$Note <- "Trait has fewer than two observed values."
    rows[[index]] <- result; index <- index + 1L
  }
  results <- dplyr::bind_rows(rows) %>% dplyr::group_by(.data$Trait) %>% dplyr::mutate(FDR = stats::p.adjust(.data$PValue, method = "BH")) %>% dplyr::ungroup()
  structure(list(results = results, eigengenes = mes, trait_data = traits, join_audit = audit), class = "feature_module_trait_obj")
}

#' Plot feature WGCNA soft-threshold diagnostics
#' @param wgcna_obj A `feature_wgcna_obj`.
#' @return A ggplot object.
#' @export
PlotFeatureWGCNASoftThreshold <- function(wgcna_obj) {
  if (!inherits(wgcna_obj, "feature_wgcna_obj")) stop("`wgcna_obj` must be a `feature_wgcna_obj`.", call. = FALSE)
  fit <- wgcna_obj$soft_threshold$fit_indices
  df_plot <- dplyr::bind_rows(tibble::tibble(Power = fit$Power, Value = -sign(fit$slope) * fit$SFT.R.sq, Metric = "Scale-free topology fit"), tibble::tibble(Power = fit$Power, Value = fit$mean.k., Metric = "Mean connectivity"))
  ggplot2::ggplot(df_plot, ggplot2::aes(.data$Power, .data$Value)) + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~Metric, scales = "free_y") + ggplot2::theme_minimal()
}

#' Plot a feature WGCNA dendrogram
#' @param wgcna_obj A `feature_wgcna_obj`.
#' @return Invisibly returns `wgcna_obj` after drawing base graphics.
#' @export
PlotFeatureWGCNADendrogram <- function(wgcna_obj) {
  if (!inherits(wgcna_obj, "feature_wgcna_obj") || !requireNamespace("WGCNA", quietly = TRUE)) stop("A `feature_wgcna_obj` and `WGCNA` are required.", call. = FALSE)
  WGCNA::plotDendroAndColors(wgcna_obj$gene_tree, wgcna_obj$module_colors, groupLabels = "Module", dendroLabels = FALSE, addGuide = TRUE, hang = 0.03)
  invisible(wgcna_obj)
}

#' Plot feature module-trait associations
#' @param module_trait_obj A `feature_module_trait_obj`.
#' @return A ggplot heatmap.
#' @export
PlotFeatureModuleTraitHeatmap <- function(module_trait_obj) {
  if (!inherits(module_trait_obj, "feature_module_trait_obj")) stop("`module_trait_obj` must be a `feature_module_trait_obj`.", call. = FALSE)
  ggplot2::ggplot(module_trait_obj$results, ggplot2::aes(.data$Trait, .data$Module, fill = .data$Effect)) + ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(.data$FDR), "", formatC(.data$FDR, digits = 2, format = "f"))), size = 3) + ggplot2::scale_fill_gradient2(na.value = "grey90") + ggplot2::theme_minimal()
}
