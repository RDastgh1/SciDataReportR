test_that("CreateClusterModel_SOM_MClust returns the complete model-review output set", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory",
    k_range = 2:4,
    models = 1,
    lpa_timeout_seconds = 30
  )

  expect_s3_class(model$fit_plot, "ggplot")
  expect_s3_class(model$ModelInfo_SOM$plots$Circular, "htmlwidget")
  expect_s3_class(model$ModelInfo_SOM$plots$Line, "htmlwidget")
  expect_s3_class(model$ModelInfo_SOM$plots$Cloud, "htmlwidget")

  som_fit_plots <- model$ModelInfo_SOM$SOMFit$plots
  reference_plots <- model$ModelInfo_SOM$ProjectionReference$plots
  probability_plots <- model$ProbFit$plots

  expect_true(all(vapply(som_fit_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(reference_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(probability_plots, inherits, logical(1), what = "ggplot")))
  expect_false(any(grepl("MaxProbBoxplot", names(probability_plots))))
  expect_true(all(c(
    "node_ProbAssignedDensity", "individual_ProbAssignedDensity"
  ) %in% names(probability_plots)))
  expect_equal(
    probability_plots$node_ProbAssignedDensity$scales$get_scales("x")$limits,
    c(0, 1)
  )
  expect_equal(
    probability_plots$individual_ProbAssignedDensity$scales$get_scales("x")$limits,
    c(0, 1)
  )

  fit_table <- model$ModelInfo_MClust$fit_table
  expect_false("CAICCLC" %in% names(fit_table))
  expect_true(all(c("CAIC", "CLC") %in% names(fit_table)))
  expect_true(all(c("MinProfileNodeN", "MaxProfileNodeN",
    "MinProfileNodeProportion", "MaxProfileNodeProportion",
    "BLRTStatistic", "BLRTPValue") %in% names(fit_table)))
  expect_true(all(fit_table$MinProfileNodeN == as.integer(fit_table$MinProfileNodeN)))
  expect_true(all(fit_table$MaxProfileNodeN == as.integer(fit_table$MaxProfileNodeN)))
  expect_false(any(c("n_min", "n_max", "BLRT_val", "BLRT_p") %in% names(fit_table)))
  expect_true("BLRTPValue" %in% unique(model$fit_plot$data$name))
  expected_best <- fit_table[which.max(fit_table$ahp_index), , drop = FALSE]
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Model, expected_best$Model)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, expected_best$Classes)
})

test_that("CreateClusterModel_SOM_MClust handles a single exploratory candidate", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory",
    k_range = 2,
    models = 1,
    lpa_timeout_seconds = 30
  )

  expect_equal(nrow(model$ModelInfo_MClust$fit_table), 1)
  expect_equal(model$ModelInfo_MClust$fit_table$ahp_index, 0)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, 2)
})

test_that("CreateClusterModel_SOM_MClust adds deterministic bootstrap stability to fit review", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  fit_model <- function() {
    CreateClusterModel_SOM_MClust(
      data = df_Labelled,
      variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
      method = "exploratory",
      k_range = 2,
      models = 1,
      stability_resamples = 2,
      stability_seed = 1234,
      lpa_timeout_seconds = 30
    )
  }

  model_one <- fit_model()
  model_two <- fit_model()
  stability <- model_one$ModelInfo_MClust$Stability
  fit_table <- model_one$ModelInfo_MClust$fit_table

  expect_equal(stability$settings$resamples, 2)
  expect_equal(nrow(stability$replicates), 2)
  expect_true(all(stability$replicates$Status == "success"))
  expect_true(all(stability$cluster_recovery$Jaccard >= 0))
  expect_true(all(stability$cluster_recovery$Jaccard <= 1))
  expect_true(all(c(
    "StabilitySuccessRate", "StabilityARI_Mean", "StabilityARI_P05",
    "StabilityJaccard_Mean", "StabilityJaccard_Min",
    "ReproducibilityScore", "Reproducibility_scaled"
  ) %in% names(fit_table)))
  expect_false("CAICCLC" %in% names(fit_table))
  expect_true(all(c("CAIC", "CLC") %in% names(fit_table)))
  expect_false(any(grepl("MaxProbBoxplot", names(model_one$ProbFit$plots))))
  expect_true("ReproducibilityScore" %in% unique(model_one$fit_plot$data$name))
  expect_equal(
    model_one$ModelInfo_MClust$Stability$replicates$ARI,
    model_two$ModelInfo_MClust$Stability$replicates$ARI
  )
  expect_equal(mclust::adjustedRandIndex(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  negative_ari <- mclust::adjustedRandIndex(
    c(1, 1, 1, 2), c(1, 2, 2, 2)
  )
  expect_lt(negative_ari, 0)
  expect_equal(mean(c(negative_ari, 0.8)), (negative_ari + 0.8) / 2)
})
