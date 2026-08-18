test_that("AssemblePlots combines categorical and continuous PlotAssociations outputs", {
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  p_Categorical <- PlotAssociations(df_Labelled, "Diagnosis", "Genotype")
  p_Continuous <- PlotAssociations(df_Labelled, "age", "AXL")

  actual <- AssemblePlots(
    list(
      "Diagnosis and genotype" = p_Categorical,
      "Age and AXL" = p_Continuous
    ),
    ncol = 2,
    CollectLegend = FALSE,
    Labels = c("A", "B")
  )

  expect_s3_class(p_Categorical, "ggplot")
  expect_s3_class(p_Continuous, "ggplot")
  expect_s3_class(actual, "ggplot")
})
