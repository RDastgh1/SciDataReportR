SOMKMeansExampleData <- function() {
  set.seed(521)
  df_Example <- data.frame(x = c(rnorm(24, -3, .5), rnorm(24, 3, .5), rnorm(24, 0, .5)),
    y = c(rnorm(24, -2, .5), rnorm(24, 2, .5), rnorm(24, 5, .5)))
  attr(df_Example$x, "label") <- "First clinical marker"
  attr(df_Example$y, "label") <- "Second clinical marker"
  df_Example
}

test_that("node whitening equals the frozen regularized Mahalanobis metric", {
  set.seed(289)
  codes <- cbind(x = rnorm(25), y = rnorm(25))
  codes <- cbind(codes, z = codes[, 1] + codes[, 2])
  whitening <- .FitSOMNodeWhitening(codes, .05)
  transformed <- .ApplySOMNodeWhitening(codes, whitening)
  expected_covariance <- .95 * stats::cov(codes) +
    .05 * diag(mean(diag(stats::cov(codes))), 3)
  expect_equal(whitening$covariance, expected_covariance)
  expect_equal(rowSums(transformed^2), unname(stats::mahalanobis(
    codes, whitening$center, whitening$covariance)))
  for (i in seq_len(nrow(codes))) {
    expect_equal(as.numeric(as.matrix(stats::dist(transformed))[, i]^2),
      unname(stats::mahalanobis(codes, codes[i, ], whitening$covariance)))
  }
  expect_equal(sweep(transformed %*% whitening$inverse_transform, 2L,
    whitening$center, "+"), codes, ignore_attr = TRUE)
  expect_error(.FitSOMNodeWhitening(codes, 0), "positive definite")
  expect_error(.FitSOMNodeWhitening(matrix(1, 5, 3), .05), "no variation")
  expect_error(.FitSOMNodeWhitening(codes, 1.1), "between 0 and 1")
  expect_error(.FitSOMNodeWhitening(codes, NA_real_), "between 0 and 1")
  ordinary <- .FitSOMNodeWhitening(codes[, 1:2], 0)
  expect_equal(stats::cov(.ApplySOMNodeWhitening(codes[, 1:2], ordinary)),
    diag(2), ignore_attr = TRUE)
})

test_that("signed centroid margins preserve inherited labels", {
  X <- rbind(c(.1, 0), c(2, 0), c(1, 0), c(NA_real_, 0))
  centers <- rbind(c(0, 0), c(2, 0))
  fit <- .SOMKMeansAssignmentFit(X, centers, c(2L, 2L, 1L, NA_integer_))
  expect_equal(fit$DistanceToCentroid[1:3], c(1.9, 0, 1))
  expect_equal(fit$AssignmentMargin[1:3], c(-1.8, 2, 0))
  expect_equal(fit$Flag_ProjectionDisagreement[1:3], c(TRUE, FALSE, FALSE))
  expect_true(is.na(fit$DistanceToCentroid[[4]]))
})

test_that("SOM K-means freezes the map, node transform, and standard outputs", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  df_Example <- SOMKMeansExampleData()
  df_Example$x[[4]] <- NA_real_
  model <- suppressWarnings(CreateClusterModel_SOM_KMeans(df_Example, c("x", "y"),
    method = "finalize", final_k = 3, som_xdim = 3, som_ydim = 3,
    som_epochs = 10L, seed_som = 102L, seed_kmeans = 45L, Relabel = FALSE))
  expect_s3_class(model, "Pipeline_SOM_KMeans")
  expect_identical(model$ModelInfo, model$ModelInfo_KMeans)
  expect_identical(model$Preprocessing$Whitening, model$ModelInfo$Whitening)
  expect_null(model$Specification$seeds$seed_mixture)
  expect_equal(model$Specification$seeds$seed_kmeans, 45L)
  expect_true(all(c("ModelInfo_SOM", "Preprocessing", "Specification", "ProbFit",
    "fit_plot", "DataWithClusters") %in% names(model)))
  expect_equal(sum(model$ModelInfo$cluster_sizes$ParticipantN), 71L)
  expect_equal(sum(model$ModelInfo$cluster_sizes$NodeN), 9L)
  expect_true(is.na(model$DataWithClusters$Cluster[[4]]))
  expect_equal(model$DataWithClusters$Cluster,
    unname(model$ModelInfo$node_cluster[model$ProbFit$individual$SOM_Node]))
  codes <- model$ModelInfo_SOM$som_codes
  whitened <- .ApplySOMNodeWhitening(codes, model$ModelInfo$Whitening)
  set.seed(48L)
  reference <- stats::kmeans(whitened, 3, nstart = 50)
  expect_equal(model$ModelInfo$node_cluster, reference$cluster)
  expect_equal(model$ModelInfo$fit_table$WSS, reference$tot.withinss)
  expect_equal(model$ModelInfo$fit_table$Silhouette,
    mean(cluster::silhouette(reference$cluster, stats::dist(whitened))[, "sil_width"]))
  expect_equal(model$ModelInfo$centers,
    do.call(rbind, lapply(split(seq_len(nrow(codes)), reference$cluster),
      function(rows) colMeans(codes[rows, , drop = FALSE]))), ignore_attr = TRUE)
  # Changing node clustering must not change the underlying SOM or BMUs.
  baseline <- suppressWarnings(CreateClusterModel_SOM_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3, som_xdim = 3, som_ydim = 3,
    som_epochs = 10L, seed_som = 102L, Relabel = FALSE))
  expect_equal(model$ModelInfo_SOM$som_codes, baseline$ModelInfo_SOM$som_codes)
  expect_equal(model$ProbFit$individual$SOM_Node, baseline$ProbFit$individual$SOM_Node)
  expect_equal(model$ModelInfo$plots$silhouette$labels$x, "SOM node")
  expect_equal(model$ModelInfo$plots$silhouette$labels$y, "Silhouette width")
  expect_equal(baseline$ModelInfo$plots$silhouette$labels$x, "SOM node")
  expect_equal(baseline$ModelInfo$plots$silhouette$labels$y, "Silhouette width")
  for (plot in model$ModelInfo$plots) expect_silent(ggplot2::ggplot_build(plot))
  heatmap <- model$ModelInfo$plots$centre_heatmap$data
  expect_true(all(c("First clinical marker", "Second clinical marker") %in% as.character(heatmap$Variable)))
  original <- model
  df_New <- df_Example[c(1, 4, 35, 60), ]
  df_New$x[[1]] <- 100
  projected <- suppressMessages(suppressWarnings(ProjectCluster(model, df_New)))
  expect_identical(model, original)
  expect_s3_class(projected, "Project_SOM_KMeans")
  expect_equal(projected$ProbFit$individual$Cluster,
    unname(model$ModelInfo$node_cluster[projected$ProbFit$individual$SOM_Node]))
  expect_true(is.na(projected$DataWithClusters$Cluster[[2]]))
  expect_equal(nrow(projected$ProjectionFit$individual), 4L)
  expect_true(all(c("SOM_Distance", "SOMDist_percentile_train", "Flag_FarFromMap",
    "DistanceToCentroid", "NodeDistanceToCentroid", "AssignmentMargin",
    "NodeAssignmentMargin", "Flag_ProjectionDisagreement") %in%
    names(projected$ProjectionFit$individual)))
  expect_true(nrow(projected$ProjectionFit$out_of_support) > 0L)
  expect_true(all(is.na(projected$ProbFit$individual$max_prob)))
  expect_equal(projected$Specification$projection$membership_probability, "unavailable")
  expect_equal(projected$ProjectionFit$centroid_fit$policy$high_distance_cutoff,
    model$ModelInfo$FitDiagnostics$distance_reference$high_distance_cutoff)
  X <- .SOMKMeansFeatures(model, df_New)
  for (i in c(1, 3, 4)) {
    cluster <- projected$ProbFit$individual$Cluster[[i]]
    expect_equal(projected$ProbFit$individual$DistanceToCentroid[[i]]^2,
      unname(stats::mahalanobis(X[i, , drop = FALSE],
        model$ModelInfo$centers[cluster, ], model$ModelInfo$Whitening$covariance)))
  }
})

test_that("SOM K-means exploration and stability refit the full pipeline", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  df_Example <- SOMKMeansExampleData()
  fitted_covariances <- list()
  FitWhitening <- .FitSOMNodeWhitening
  local_mocked_bindings(.FitSOMNodeWhitening = function(codes, covariance_shrinkage) {
    result <- FitWhitening(codes, covariance_shrinkage)
    fitted_covariances[[length(fitted_covariances) + 1L]] <<- result$covariance
    result
  }, .package = "SciDataReportR")
  review <- suppressWarnings(CreateClusterModel_SOM_KMeans(df_Example, c("x", "y"),
    k_range = 2:3, som_xdim = 3, som_ydim = 3, som_epochs = 5L,
    stability_resamples = 2L, stability_cores = 1L, Relabel = FALSE))
  expect_null(review$DataWithClusters)
  expect_null(review$ProbFit)
  expect_equal(review$ModelInfo$fit_table$Classes, 2:3)
  expect_error(ProjectCluster(review, df_Example), "finalized")
  expect_equal(length(fitted_covariances), 5L)
  expect_false(isTRUE(all.equal(fitted_covariances[[1]], fitted_covariances[[2]])))
  expect_equal(nrow(review$Stability$replicates), 4L)
  expect_equal(nrow(review$Stability$failures), 0L)
  expect_equal(review$Stability$settings$whitening, "refit_from_subsample_SOM_nodes")
  expect_equal(review$Stability$settings$som_grid$policy, "fixed_to_reference_fit")
  expect_error(CreateClusterModel_SOM_KMeans(df_Example, final_k = 3), "finalize")
  expect_error(CreateClusterModel_SOM_KMeans(df_Example, method = "finalize"), "final_k")
  expect_error(CreateClusterModel_SOM_KMeans(df_Example, nstart = 0), "nstart")
  expect_error(CreateClusterModel_SOM_KMeans(df_Example, assignment = "soft"), "managed")
  expect_error(CreateClusterModel_SOM_KMeans(df_Example, mixture_bootstrap_n = 10), "managed")
})

test_that("shared review labels describe training PCs and observation units", {
  space <- .ClusterReviewSpace(as.matrix(SOMKMeansExampleData()))
  expect_match(.ClusterReviewAxisLabels(space)[[1]], "Training PC1.*variance")
  expect_match(.ClusterReviewAxisLabels(space)[[2]], "Training PC2.*variance")
  expect_equal(.ClusterReviewAxisLabels(list()), c("Training PC1", "Training PC2"))
  silhouette <- cluster::silhouette(rep(1:3, each = 24), stats::dist(SOMKMeansExampleData()))
  expect_equal(PlotClusterSilhouette(silhouette)$labels$x, "Participant")
  expect_equal(PlotClusterSilhouette(silhouette, observation_unit = "SOM node")$labels$x, "SOM node")
  expect_equal(PlotClusterSilhouette(silhouette, observation_unit = "SOM node")$labels$y, "Silhouette width")
})
