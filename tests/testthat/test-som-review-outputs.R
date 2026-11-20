# Review, map and assignment outputs for SOM + Mclust. Every model is fitted on
# the embedded SimulatedPhenotypeData only.

.ReviewModels <- new.env(parent = emptyenv())

ReviewFixture <- function(name = "default", ...) {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  if (is.null(.ReviewModels[[name]])) {
    set.seed(2611)
    .ReviewModels[[name]] <- suppressMessages(suppressWarnings(
      CreateClusterModel_SOM_MClust(EmbeddedTestData(), paste0("Var", 1:3),
        k_range = 2:3, models = c(1, 6), som_xdim = 5, som_ydim = 5,
        Relabel = FALSE, mixture_bootstrap_n = 0, ...)))
  }
  .ReviewModels[[name]]
}

test_that("SOM review shows the model-selection criteria and marks the recommendation", {
  model <- ReviewFixture("stability", stability_resamples = 2)
  plot <- model$fit_plot
  expect_identical(levels(plot$data$Metric),
    c("AIC", "BIC", "Entropy", "BoundaryPct", "StabilityARI_Mean"))
  expect_false("ICL" %in% plot$data$Metric)
  expect_identical(
    plot$facet$params$labeller(data.frame(Metric = levels(plot$data$Metric)))[[1]],
    c("AIC (lower is better)", "BIC (lower is better)", "Entropy (higher is better)",
      "Boundary % (lower is better)", "Mean ARI across refits (higher is better)"))
  boundary <- plot$data[plot$data$Metric == "BoundaryPct", ]
  expect_equal(boundary$PlotValue, 100 * boundary$Value)

  ring <- Filter(function(layer) identical(layer$aes_params$shape, 21), plot$layers)
  expect_length(ring, 1L)
  expect_true(all(ring[[1]]$data$Recommended))
  expect_identical(nrow(ring[[1]]$data), nlevels(plot$data$Metric))
  expect_match(plot$labels$caption, "Ringed", fixed = TRUE)
  over_parameterized <- any(model$ModelInfo_MClust$fit_table$SupportStatus == "over-parameterized")
  expect_identical(grepl("Cross", plot$labels$caption, fixed = TRUE), over_parameterized)
  expect_identical(plot$labels$colour, "Model")

  # One figure whichever way it is requested.
  expect_equal(PlotClusterFitReview(model)$data, plot$data)
})

test_that("a SOM fit table that lost its recorded attributes keeps the SOM review", {
  model <- ReviewFixture()
  bare <- model$ModelInfo_MClust$fit_table
  attr(bare, "ReviewMetrics") <- NULL
  attr(bare, "MetricDirection") <- NULL
  plot <- PlotClusterFitReview(bare)
  expect_identical(levels(plot$data$Metric), c("AIC", "BIC", "Entropy", "BoundaryPct"))
  expect_identical(.ClusterMetricDirection("BIC", bare), "lower")
})

test_that("review marks candidates excluded by min_cluster_prop", {
  fit_table <- data.frame(Classes = 2:4, Model = 1L,
    ModelLabel = "Model 1: equal variance, zero covariance", BIC = c(30, 20, 25),
    SupportStatus = "supported", Eligible = c(TRUE, TRUE, FALSE),
    Recommended = c(FALSE, TRUE, FALSE), MixtureNodes = 25L)
  plot <- PlotClusterFitReview(fit_table, metrics = "BIC")
  expect_identical(as.character(plot$data$SupportStatus),
    c("supported", "supported", "ineligible"))
  expect_match(plot$labels$caption, "Star", fixed = TRUE)
  expect_identical(.ClusterMetricDirection("BIC", fit_table), "lower")
})

test_that("model-selection table lists each candidate's decision criteria", {
  model <- ReviewFixture("stability", stability_resamples = 2)
  fits <- model$ModelInfo_MClust$fit_table
  selection <- ClusterModelSelectionTable(model)
  expect_identical(names(selection), c("Model", "k", "AIC (lower is better)",
    "BIC (lower is better)", "Entropy (higher is better)",
    "Boundary % (lower is better)", "Mean ARI across refits (higher is better)",
    "Cluster sizes", "Eligible", "AHP pick"))
  expect_identical(nrow(selection), nrow(fits))
  expect_equal(selection$k, sort(fits$Classes))
  expect_true(all(grepl("%$", selection$`Boundary % (lower is better)`)))
  picked <- fits[fits$Recommended, , drop = FALSE]
  expect_identical(selection$Model[selection$`AHP pick` != ""], picked$ModelName)

  # Sizes cover every complete participant, and the recommended candidate's
  # sizes are the clusters of the returned model.
  n_complete <- sum(!is.na(model$ProbFit$individual$SOM_Node))
  totals <- vapply(strsplit(fits$ClusterSizes, " / "),
    function(sizes) sum(as.integer(sizes)), numeric(1))
  expect_true(all(totals == n_complete))
  final_sizes <- paste(sort(as.integer(table(model$DataWithClusters$Cluster)),
    decreasing = TRUE), collapse = " / ")
  expect_identical(picked$ClusterSizes, final_sizes)
})

test_that("model-selection table works for other pipelines", {
  df_Training <- EmbeddedTestData(variables = paste0("Var", 1:4))
  model <- CreateClusterModel_KMeans(df_Training, names(df_Training), k_range = 2:4)
  selection <- ClusterModelSelectionTable(model)
  expect_equal(selection$k, 2:4)
  expect_true("Silhouette (higher is better)" %in% names(selection))
  expect_identical(sum(selection$`AHP pick` != ""), 1L)
})

test_that("border nodes have a neighbour in another cluster", {
  skip_if_not_installed("kohonen")
  pts <- kohonen::somgrid(3, 2, "hexagonal")$pts
  adjacent <- .SOMGridAdjacency(pts)
  # kohonen offsets the bottom row, so node 1 touches 2, 4 and 5.
  expect_identical(unname(which(adjacent[1, ])), c(2L, 4L, 5L))
  clusters <- c(1, 1, 2, 1, 1, 2)
  expect_identical(.SOMBorderNodes(clusters, adjacent),
    c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE))

  # Borders are drawn on the hexagon edges shared by the two nodes.
  segments <- .SOMBorderSegments(pts, clusters, adjacent)
  expect_identical(nrow(segments), 3L)
  expect_equal(sqrt((segments$xend - segments$x)^2 + (segments$yend - segments$y)^2),
    rep(1 / sqrt(3), 3))
  shapes <- .SOMNodeShapes(pts)
  expect_equal(sqrt((shapes$x - pts[shapes$NodeID, 1])^2 + (shapes$y - pts[shapes$NodeID, 2])^2),
    rep(1 / sqrt(3), nrow(shapes)))
  VertexCount <- function(x, y) sum(abs(shapes$x - x) < 1e-9 & abs(shapes$y - y) < 1e-9)
  expect_true(all(mapply(VertexCount, c(segments$x, segments$xend),
    c(segments$y, segments$yend)) >= 2))

  model <- ReviewFixture()
  node <- model$ProbFit$node
  expect_identical(node$Flag_BorderNode,
    .SOMBorderNodes(node$Cluster, .SOMGridAdjacency(model$ModelInfo_SOM$som_grid$pts)))
  expect_true(any(node$Flag_BorderNode))
})

test_that("PlotSOMMap draws the map for every fill", {
  model <- ReviewFixture()
  n_nodes <- nrow(model$ProbFit$node)
  for (fill in c("Cluster", "People", "Boundary", "UMatrix", "Var1", "Z_Var2")) {
    built <- ggplot2::ggplot_build(PlotSOMMap(model, fill = fill))
    expect_identical(length(unique(built$data[[1]]$group)), n_nodes, info = fill)
  }

  plot <- PlotSOMMap(model)
  pts <- model$ModelInfo_SOM$som_grid$pts
  adjacent <- .SOMGridAdjacency(pts)
  clusters <- model$ProbFit$node$Cluster
  differing <- sum(adjacent & upper.tri(adjacent) & outer(clusters, clusters, "!="))
  segment_layers <- Filter(function(layer) inherits(layer$geom, "GeomSegment"), plot$layers)
  expect_identical(nrow(segment_layers[[1]]$data), differing)

  people <- PlotSOMMap(model, fill = "People")$data
  expect_equal(people$Value[!duplicated(people$NodeID)], model$ProbFit$node$NodeN)
  boundary <- PlotSOMMap(model, fill = "Boundary")$data
  individual <- model$ProbFit$individual
  expected <- tapply(as.numeric(individual$Flag_BoundaryParticipant),
    factor(individual$SOM_Node, levels = seq_len(n_nodes)), mean)
  expect_equal(boundary$Value[!duplicated(boundary$NodeID)], as.numeric(expected))

  # Highlighted nodes get a white halo and a black outline, drawn last.
  outlines <- utils::tail(PlotSOMMap(model, highlight_nodes = c(1, 2))$layers, 2)
  expect_identical(unname(vapply(outlines, function(layer) layer$aes_params$colour, "")),
    c("white", "black"))
  expect_true(all(vapply(outlines, function(layer) {
    setequal(unique(layer$data$NodeID), c(1, 2))
  }, logical(1))))
  expect_error(PlotSOMMap(model, fill = "Unknown"), "Var1")
  expect_error(PlotSOMMap(model, highlight_nodes = n_nodes + 1), "highlight_nodes")
  expect_error(PlotSOMMap(model$ModelInfo_MClust$fit_table), "CreateClusterModel_SOM_MClust")
})

test_that("ExplainSOMAssignment reproduces the model's assignments", {
  model <- ReviewFixture()
  rows <- c(1, 7, 40)
  explanation <- ExplainSOMAssignment(model, row = rows)
  stored <- model$ProbFit$individual[rows, ]
  soft_cols <- grep("^SoftProb_", names(stored), value = TRUE)
  expect_equal(as.matrix(explanation$person[soft_cols]), as.matrix(stored[soft_cols]),
    ignore_attr = TRUE)
  expect_identical(explanation$person$Cluster_Hard, stored$Cluster_Hard)
  expect_identical(explanation$person$Cluster_Soft, stored$Cluster_Soft)
  expect_identical(explanation$person$Cluster, stored$Cluster)
  expect_identical(explanation$person$Flag_BoundaryParticipant, stored$Flag_BoundaryParticipant)
  expect_equal(as.integer(explanation$person$SOM_Node), as.integer(stored$SOM_Node))
  expect_equal(explanation$person$SOM_Distance, stored$SOM_Distance)

  nodes <- explanation$nodes
  expect_equal(as.numeric(tapply(nodes$Weight, nodes$Row, sum)), rep(1, length(rows)))
  expect_true(all(tapply(nodes$Distance, nodes$Row, function(d) !is.unsorted(d))))
  # Soft probabilities are the weight-averaged Mclust probabilities of the nodes.
  prob_cols <- grep("^Prob_", names(nodes), value = TRUE)
  blended <- t(vapply(rows, function(r) {
    colSums(nodes$Weight[nodes$Row == r] * as.matrix(nodes[nodes$Row == r, prob_cols]))
  }, numeric(length(prob_cols))))
  expect_equal(blended, as.matrix(stored[soft_cols]), ignore_attr = TRUE)

  expect_output(print(explanation), "best-matching node")
  expect_error(ExplainSOMAssignment(model, row = 0), "row")
})

test_that("ExplainSOMAssignment uses the stored Gaussian bandwidth and frozen scaling", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  model <- suppressMessages(suppressWarnings(CreateClusterModel_SOM_MClust(
    EmbeddedTestData(), paste0("Var", 1:3), method = "finalize", final_model = 1,
    final_k = 3, som_xdim = 5, som_ydim = 5, soft_kernel = "gaussian",
    Relabel = FALSE)))
  explanation <- ExplainSOMAssignment(model, row = 1:10)
  stored <- model$ProbFit$individual[1:10, ]
  soft_cols <- grep("^SoftProb_", names(stored), value = TRUE)
  expect_equal(as.matrix(explanation$person[soft_cols]), as.matrix(stored[soft_cols]),
    ignore_attr = TRUE)
  expect_identical(explanation$settings$soft_kernel, "gaussian")
  expect_equal(explanation$settings$soft_bandwidth,
    model$Specification$projection$soft_bandwidth)

  df_New <- EmbeddedTestData("Projection")
  projected <- suppressMessages(suppressWarnings(ProjectCluster(model, data = df_New)))
  new_explanation <- ExplainSOMAssignment(model, row = 1:5, data = df_New)
  expect_equal(as.matrix(new_explanation$person[soft_cols]),
    as.matrix(projected$ProbFit$individual[1:5, soft_cols]), ignore_attr = TRUE)
  expect_identical(new_explanation$person$Cluster_Soft,
    projected$ProbFit$individual$Cluster_Soft[1:5])

  df_Missing <- df_New[1:2, ]
  df_Missing$Var1[2] <- NA
  expect_error(ExplainSOMAssignment(model, row = 2, data = df_Missing), "missing")
})
