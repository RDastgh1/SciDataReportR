test_that("CreateSOMClusterModel returns the complete model-review output set", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateSOMClusterModel(
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
  reference_plots <- model$ModelInfo_SOM$PhenotypeReference$plots
  probability_plots <- model$ProbFit$plots

  expect_true(all(vapply(som_fit_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(reference_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(probability_plots, inherits, logical(1), what = "ggplot")))

  fit_table <- model$ModelInfo_MClust$fit_table
  expected_best <- fit_table[which.max(fit_table$ahp_index), , drop = FALSE]
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Model, expected_best$Model)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, expected_best$Classes)
})

test_that("CreateSOMClusterModel handles a single exploratory candidate", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateSOMClusterModel(
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

test_that("CreateSOMClusterModel adds deterministic bootstrap stability to fit review", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  fit_model <- function() {
    CreateSOMClusterModel(
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
  expect_true("ReproducibilityScore" %in% unique(model_one$fit_plot$data$name))
  expect_equal(
    model_one$ModelInfo_MClust$Stability$replicates$ARI,
    model_two$ModelInfo_MClust$Stability$replicates$ARI
  )
  expect_equal(mclust::adjustedRandIndex(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
})
