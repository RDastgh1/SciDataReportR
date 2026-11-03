test_that("soft SOM probabilities match independent distance-weighted calculations", {
  df_Nodes <- tibble::tibble(NodeID = 1:3, Cluster = c(1L, 2L, 2L),
    prob_1 = c(0.51, 0, 0), prob_2 = c(0.49, 1, 1))
  codes <- matrix(c(0, 1, 2), ncol = 1)
  participants <- matrix(c(0.1, 0.4, 1.2), ncol = 1)
  result <- .SOMSoftProjection(participants, codes, df_Nodes, n_nodes = 3L)

  # Calculate each participant separately so row/nearest-node indexing errors
  # cannot hide behind a second implementation of the vectorized calculation.
  expected <- t(vapply(as.vector(participants), function(value) {
    weights <- 1 / (value - as.vector(codes))^2
    weights <- weights / sum(weights)
    c(sum(weights * df_Nodes$prob_1), sum(weights * df_Nodes$prob_2))
  }, numeric(2)))
  expect_equal(unname(result$probability), expected)
  expect_identical(result$cluster, c(1L, 2L, 2L))
  # The first two participants share a best-matching node, but soft assignment
  # changes the second participant to class 2 because of neighbouring nodes.
  expect_identical(df_Nodes$Cluster[c(1, 1, 2)], c(1L, 1L, 2L))

  reversed <- .SOMSoftProjection(participants[3:1, , drop = FALSE], codes,
    df_Nodes, n_nodes = 3L)
  expect_equal(reversed$probability, result$probability[3:1, , drop = FALSE])
  nearest_only <- .SOMSoftProjection(participants, codes, df_Nodes, n_nodes = 1L)
  expect_identical(nearest_only$cluster, c(1L, 1L, 2L))
  at_node <- .SOMSoftProjection(codes[1, , drop = FALSE], codes, df_Nodes)
  expect_equal(unname(at_node$probability), matrix(c(0.51, 0.49), nrow = 1))
})

test_that("soft SOM ties use the actual maximum and resolve exact ties deterministically", {
  codes <- matrix(c(0, 1, 2), ncol = 1)
  df_Nodes <- tibble::tibble(NodeID = 1:3, Cluster = c(1L, 2L, 1L),
    prob_1 = c(0.5, 0.499999, 0.5), prob_2 = c(0.5, 0.500001, 0.5))
  for (seed in 1:5) {
    set.seed(seed)
    result <- .SOMSoftProjection(codes, codes, df_Nodes)
    expect_identical(result$cluster, c(1L, 2L, 1L))
  }
})

test_that("hard and soft modes share the fit but select their own participant assignments", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  models <- lapply(c("hard", "soft"), function(assignment) {
    suppressWarnings(CreateClusterModel_SOM_MClust(
      df_Training, paste0("Var", 1:6), method = "finalize",
      final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
      assignment = assignment, stability_resamples = 0L, Relabel = FALSE,
      min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30))
  })
  expect_equal(models[[1]]$ModelInfo_SOM$som_codes, models[[2]]$ModelInfo_SOM$som_codes)
  expect_equal(models[[1]]$ProbFit$node, models[[2]]$ProbFit$node)
  expect_equal(models[[1]]$ModelInfo_MClust$fit_table,
    models[[2]]$ModelInfo_MClust$fit_table)

  for (index in seq_along(models)) {
    model <- models[[index]]
    individual <- model$ProbFit$individual
    selected <- if (index == 1L) individual$Cluster_Hard else individual$Cluster_Soft
    expect_identical(individual$Cluster, selected)
    expect_identical(model$DataWithClusters$Cluster, selected)
    expect_identical(individual$Flag_NodeDisagreement,
      individual$Cluster != individual$NodeCluster)
    expect_identical(models[[1]]$ProbFit$individual$Cluster_Soft, individual$Cluster_Soft)
    expect_identical(models[[1]]$ProbFit$individual$Cluster_Hard, individual$Cluster_Hard)

    projected <- suppressWarnings(ProjectCluster(model, df_Training))
    expect_identical(projected$ProbFit$individual$Cluster, selected)
    expect_identical(projected$DataWithClusters$Cluster, selected)
    expect_identical(projected$ProbFit$individual$Flag_NodeDisagreement,
      individual$Flag_NodeDisagreement)
  }

  # A controlled boundary posterior forces disagreement through the public
  # projection API, rather than relying on the simulated clusters overlapping.
  boundary_model <- models[[2]]
  first_node <- boundary_model$ProbFit$individual$SOM_Node[1]
  boundary_model$ProbFit$node$prob_1 <- 0
  boundary_model$ProbFit$node$prob_1[first_node] <- 0.5 + 1e-12
  boundary_model$ProbFit$node$prob_2 <- 1 - boundary_model$ProbFit$node$prob_1
  boundary_model$ProbFit$node$prob_3 <- 0
  boundary_model$ProbFit$node$Cluster <- 2L
  boundary_model$ProbFit$node$Cluster[first_node] <- 1L
  df_Projection <- df_Training[1:2, ]
  df_Projection[2, paste0("Var", 1:6)] <- NA_real_
  projected <- suppressWarnings(ProjectCluster(boundary_model, df_Projection))
  individual <- projected$ProbFit$individual
  expect_identical(individual$Cluster_Hard, c(1L, NA_integer_))
  expect_identical(individual$Cluster_Soft, c(2L, NA_integer_))
  expect_identical(individual$Cluster, individual$Cluster_Soft)
  expect_identical(individual$Flag_NodeDisagreement, c(TRUE, NA))
  boundary_model$Specification$projection$assignment <- "hard"
  projected <- suppressWarnings(ProjectCluster(boundary_model, df_Projection))
  expect_identical(projected$ProbFit$individual$Cluster, individual$Cluster_Hard)
  expect_identical(projected$ProbFit$individual$Flag_NodeDisagreement, c(FALSE, NA))
})
