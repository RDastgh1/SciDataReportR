test_that("PlotClusterSplits reports a clean split as high parent purity", {
  skip_if_not_installed("ggalluvial")
  # k = 2 -> 3 splits cluster 2 into two; k = 3 -> 4 reshuffles half of cluster 1
  assignments <- data.frame(
    k2 = rep(c(1, 2), each = 100),
    k3 = c(rep(1, 100), rep(2, 50), rep(3, 50)),
    k4 = c(rep(1, 50), rep(4, 25), rep(2, 25), rep(2, 50), rep(3, 50))
  )
  splits <- PlotClusterSplits(assignments)
  transitions <- splits$transitions

  expect_s3_class(splits$plot, "ggplot")
  expect_equal(nrow(transitions), 2)
  expect_equal(transitions$ParentPurity[1], 1)
  expect_equal(transitions$SplitRetained[1], 0.5)
  expect_lt(transitions$ParentPurity[2], 1)
  expect_equal(transitions$SmallestNewCluster[1], 0.25)
  expect_no_error(ggplot2::ggplot_build(splits$plot))

  # Flows are drawn between strata (an alluvial plot, not only stacked bars)
  flow_grob <- ggplot2::layer_grob(splits$plot, i = 1)[[1]]
  expect_gt(length(flow_grob$children), 0)
})

test_that("PlotClusterSplits reads every model family from a SOM Mclust exploration", {
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled, c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory", k_range = 2:4, models = c(1, 2),
    som_xdim = 5, som_ydim = 5, min_nodes_per_cluster = NULL,
    stability_resamples = 0,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  splits <- PlotClusterSplits(model)
  expect_setequal(unique(splits$transitions$Model),
    c("Model 1: equal variance, zero covariance", "Model 2: varying variance, zero covariance"))
  expect_true(all(splits$transitions$ParentPurity > 0 & splits$transitions$ParentPurity <= 1))
  expect_no_error(ggplot2::ggplot_build(splits$plot))
  expect_equal(unique(PlotClusterSplits(model, models = 1)$transitions$Model),
    "Model 1: equal variance, zero covariance")
})
