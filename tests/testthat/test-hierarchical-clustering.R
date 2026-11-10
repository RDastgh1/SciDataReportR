HierarchicalExampleData <- function() {
  set.seed(482)
  df_Example <- data.frame(
    x = c(rnorm(18, -3, .4), rnorm(18, 3, .4), rnorm(18, 0, .4)),
    y = c(rnorm(18, -2, .4), rnorm(18, 2, .4), rnorm(18, 5, .4)))
  attr(df_Example$x, "label") <- "Marker X"
  attr(df_Example$y, "label") <- "Marker Y"
  df_Example
}

test_that("hierarchical training matches hclust and preserves the clustering contract", {
  df_Example <- HierarchicalExampleData()
  model <- CreateClusterModel_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3)
  reference <- stats::hclust(stats::dist(scale(df_Example)), method = "ward.D2")
  expect_s3_class(model, "Pipeline_Hierarchical")
  expect_equal(model$DataWithClusters$Cluster, unname(stats::cutree(reference, k = 3)))
  expect_equal(model$ModelInfo$hclust_model$height, reference$height)
  expect_identical(model$ModelInfo, model$ModelInfo_Hierarchical)
  expect_true(all(c("DataWithClusters", "ProbFit", "Specification", "Preprocessing",
    "ModelInfo", "fit_plot") %in% names(model)))
  expect_true(all(c("dendrogram", "silhouette", "silhouette_by_k", "centre_heatmap",
    "centre_profile", "profiles") %in% names(model$ModelInfo$plots)))
  expect_true(all(vapply(model$ModelInfo$plots, inherits, logical(1), "ggplot")))
  expect_equal(model$ModelInfo$FitDiagnostics$projection_agreement,
    mean(!model$ProbFit$individual$Flag_ProjectionDisagreement))
  expect_false(any(grepl("prob_assigned|max_prob", names(model$ProbFit$individual))))
  expect_true(all(model$ModelInfo$medoid_row_ids %in% seq_len(nrow(df_Example))))
  for (plot in model$ModelInfo$plots) expect_silent(ggplot2::ggplot_build(plot))
})

test_that("candidate cuts share one tree and exploration cannot be projected", {
  model <- CreateClusterModel_Hierarchical(HierarchicalExampleData(), c("x", "y"),
    k_range = 2:4)
  expect_equal(model$ModelInfo$fit_table$Classes, 2:4)
  expect_null(model$DataWithClusters)
  expect_null(model$ProbFit)
  expect_error(ProjectCluster(model, HierarchicalExampleData()), "finalized")
  expect_equal(sum(model$ModelInfo$fit_table$Recommended), 1L)
})

test_that("hierarchical input and lifecycle choices are validated", {
  df_Example <- HierarchicalExampleData()
  expect_error(CreateClusterModel_Hierarchical(df_Example, final_k = 3), "finalize")
  expect_error(CreateClusterModel_Hierarchical(df_Example, method = "finalize"), "final_k")
  expect_error(CreateClusterModel_Hierarchical(df_Example, method = "finalize",
    final_k = c(2, 3)), "single")
  expect_error(CreateClusterModel_Hierarchical(df_Example, k_range = c(1, 3)), "Cluster counts")
  expect_error(CreateClusterModel_Hierarchical(df_Example, k_range = 2.5), "Cluster counts")
  expect_error(CreateClusterModel_Hierarchical(df_Example, distance = "manhattan"), "Ward")
  expect_error(CreateClusterModel_Hierarchical(df_Example, linkage = "bogus"), "arg")
  df_Example$x[[1]] <- Inf
  expect_error(CreateClusterModel_Hierarchical(df_Example, ZScoreType = "None"), "finite")
})

test_that("Manhattan linkage uses actual within-cluster medoids", {
  df_Example <- HierarchicalExampleData()
  model <- CreateClusterModel_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3, distance = "manhattan", linkage = "average",
    ZScoreType = "None")
  reference <- stats::hclust(stats::dist(df_Example, method = "manhattan"), method = "average")
  expect_equal(model$DataWithClusters$Cluster, unname(stats::cutree(reference, 3)))
  expect_true(is.na(model$ModelInfo$fit_table$CalinskiHarabasz))
  for (group in 1:3) {
    rows <- which(model$DataWithClusters$Cluster == group)
    expected <- rows[[which.min(rowSums(as.matrix(stats::dist(df_Example[rows, ], method = "manhattan"))))]]
    expect_equal(model$ModelInfo$medoid_row_ids[[group]], expected)
  }
})

test_that("projection freezes preprocessing and retains fit, margins, ranges, and missing rows", {
  df_Example <- HierarchicalExampleData()
  model <- CreateClusterModel_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3)
  original <- model
  df_New <- df_Example[c(4, 25, 44), ]
  df_New$x[[1]] <- 100
  df_New$x[[3]] <- NA_real_
  projected <- suppressWarnings(ProjectCluster(model, df_New))
  expect_identical(model, original)
  expect_equal(nrow(projected$DataWithClusters), 3L)
  expect_true(is.na(projected$DataWithClusters$Cluster[[3]]))
  expect_true(projected$ProjectionFit$individual$Flag_HighDistance[[1]])
  expect_equal(projected$ProjectionFit$individual$DistancePercentileTrain[[1]], 1)
  expect_true(all(projected$ProbFit$individual$AssignmentMargin[1:2] >= 0))
  expect_true("x" %in% projected$ProjectionFit$out_of_support$Variable)
  expect_true(all(c("individual", "summary", "by_cluster", "policy", "plots") %in%
    names(projected$ProjectionFit)))
  expect_equal(projected$ProjectionFit$policy$high_distance_cutoff,
    model$ModelInfo$FitDiagnostics$distance_reference$high_distance_cutoff)
  expect_equal(projected$ProbFit$individual$Projection_Fit_Class,
    projected$DataWithClusters$Projection_Fit_Class)
  df_Missing <- data.frame(x = c(NA_real_, NA_real_), y = c(NA_real_, NA_real_))
  missing <- suppressWarnings(ProjectCluster(model, df_Missing))
  expect_true(all(is.na(missing$DataWithClusters$Cluster)))
})

test_that("hierarchical stability refits the full pipeline", {
  model <- CreateClusterModel_Hierarchical(HierarchicalExampleData(), c("x", "y"),
    method = "finalize", final_k = 3, stability_resamples = 2L,
    seed_stability = 53L, stability_cores = 1L)
  expect_equal(nrow(model$Stability$replicates), 2L)
  expect_equal(nrow(model$Stability$failures), 0L)
  expect_equal(model$Stability$settings$memory_method, "hierarchical")
  expect_equal(model$Specification$stability, model$Stability$settings)
})

test_that("SOM hierarchical labels follow frozen nodes and projection retains map fit", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  df_Example <- HierarchicalExampleData()
  model <- suppressWarnings(CreateClusterModel_SOM_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3, som_xdim = 3, som_ydim = 3,
    som_epochs = 10L, Relabel = FALSE))
  expect_s3_class(model, "Pipeline_SOM_Hierarchical")
  reference <- stats::cutree(stats::hclust(stats::dist(model$ModelInfo_SOM$som_codes),
    method = "ward.D2"), 3)
  expect_equal(model$ModelInfo$node_cluster, unname(reference))
  expect_equal(model$DataWithClusters$Cluster,
    unname(reference[model$ProbFit$individual$SOM_Node]))
  expect_equal(sum(model$ModelInfo$cluster_sizes$ParticipantN), nrow(df_Example))
  expect_equal(sum(model$ModelInfo$cluster_sizes$NodeN), 9L)
  expect_true(all(is.na(model$ProbFit$individual$max_prob)))
  expect_true("MinParticipantN" %in% names(model$ModelInfo$fit_table))
  df_New <- df_Example[1:5, ]
  df_New$x[[1]] <- 100
  projected <- suppressMessages(suppressWarnings(ProjectCluster(model, df_New)))
  expect_equal(projected$ProbFit$individual$Cluster,
    model$ModelInfo$node_cluster[projected$ProjectionFit$individual$SOM_Node])
  expect_true(projected$ProjectionFit$individual$Flag_FarFromMap[[1]])
  expect_equal(projected$DataWithClusters$Projection_Fit_Class[[1]], "Poor fit to training structure")
  expect_false(any(projected$DataWithClusters$Projection_Fit_Class == "Uncertain membership"))
  expect_true(all(c("SOM_Distance", "SOMDist_percentile_train", "NodeAssignmentMargin") %in%
    names(projected$ProjectionFit$individual)))
  expect_true(all(c("summary", "by_cluster", "out_of_support", "plots") %in% names(projected$ProjectionFit)))
  expect_equal(projected$Specification$projection$membership_probability, "unavailable")
})

test_that("SOM hierarchical exploration and stability use the shared lifecycle", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  df_Example <- HierarchicalExampleData()
  review <- suppressWarnings(CreateClusterModel_SOM_Hierarchical(df_Example, c("x", "y"),
    k_range = 2:3, som_xdim = 3, som_ydim = 3, som_epochs = 5L, Relabel = FALSE))
  expect_null(review$DataWithClusters)
  expect_equal(review$ModelInfo$fit_table$Classes, 2:3)
  expect_error(ProjectCluster(review, df_Example), "finalized")
  expect_error(CreateClusterModel_SOM_Hierarchical(df_Example, assignment = "soft"), "managed")
  model <- suppressWarnings(CreateClusterModel_SOM_Hierarchical(df_Example, c("x", "y"),
    method = "finalize", final_k = 3, som_xdim = 3, som_ydim = 3,
    som_epochs = 5L, Relabel = FALSE, stability_resamples = 2L, stability_cores = 1L))
  expect_equal(nrow(model$Stability$replicates), 2L)
  expect_equal(nrow(model$Stability$failures), 0L)
  expect_equal(model$Stability$settings$som_grid$policy, "fixed_to_reference_fit")
})
