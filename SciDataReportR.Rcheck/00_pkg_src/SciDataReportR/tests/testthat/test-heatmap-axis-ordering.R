test_that("heatmap ordering preserves requested axes by default", {
  df_Tiles <- data.frame(
    Row = rep(c("r1", "r2"), each = 3),
    Column = rep(c("c1", "c2", "c3"), 2),
    Value = c(1, 2, 3, 3, 2, 1)
  )

  order <- SciDataReportR:::.OrderHeatmapAxes(
    df_Tiles, row_id = "Row", column_id = "Column", value = "Value",
    row_order = c("r2", "r1"), column_order = c("c3", "c2", "c1")
  )

  expect_identical(order$rows, c("r2", "r1"))
  expect_identical(order$columns, c("c3", "c2", "c1"))
})

test_that("correlation heatmap returns clustered AxisOrder", {
  df_Test <- data.frame(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(1, 2, 3, 4, 5, 6),
    c = c(6, 5, 4, 3, 2, 1)
  )

  result <- PlotCorrelationsHeatmap(
    df_Test, predictor_vars = c("a", "b", "c"),
    outcome_vars = c("a", "b", "c"), cluster_rows = TRUE,
    cluster_columns = TRUE
  )

  expect_setequal(result$AxisOrder$rows, c("a", "b", "c"))
  expect_setequal(result$AxisOrder$columns, c("a", "b", "c"))
  expect_true(result$AxisOrder$cluster_rows)
  expect_true(result$AxisOrder$cluster_columns)
})

test_that("axis clustering handles masked and missing tile values", {
  df_Tiles <- data.frame(
    Row = rep(c("r1", "r2", "r3"), each = 3),
    Column = rep(c("c1", "c2", "c3"), 3),
    Value = c(NA, 0.8, -0.2, 0.8, NA, NA, -0.2, NA, NA)
  )

  order <- SciDataReportR:::.OrderHeatmapAxes(
    df_Tiles, row_id = "Row", column_id = "Column", value = "Value",
    cluster_rows = TRUE, cluster_columns = TRUE
  )

  expect_setequal(order$rows, c("r1", "r2", "r3"))
  expect_setequal(order$columns, c("c1", "c2", "c3"))
})

test_that("symmetric heatmaps support display-only triangles", {
  set.seed(20260911)
  df_Test <- data.frame(
    Group = rep(c("Reference", "Comparison"), each = 20),
    x = rnorm(40),
    y = rnorm(40),
    z = rnorm(40),
    binary_a = rep(c("No", "Yes"), 20),
    binary_b = rep(c("Low", "High"), each = 20),
    binary_c = rep(c("Absent", "Present", "Absent", "Present"), 10)
  )

  correlation <- PlotCorrelationsHeatmap(
    df_Test,
    predictor_vars = c("x", "y", "z"),
    outcome_vars = c("x", "y", "z"),
    triangle = "upper",
    cluster_rows = TRUE
  )
  expect_equal(nrow(correlation$Unadjusted$plot$data), 3)
  expect_identical(correlation$AxisOrder$rows, correlation$AxisOrder$columns)
  expect_identical(correlation$AxisOrder$triangle_applied, "upper")
  expect_equal(dim(correlation$Unadjusted$r), c(3, 3))

  rectangular <- PlotCorrelationsHeatmap(
    df_Test,
    predictor_vars = c("x", "y"),
    outcome_vars = "z",
    triangle = "upper"
  )
  expect_equal(nrow(rectangular$Unadjusted$plot$data), 2)
  expect_identical(rectangular$AxisOrder$triangle_applied, "full")

  comparison <- PlotCorrelationComparisons(
    df_Test,
    predictor_vars = c("x", "y", "z"),
    outcome_vars = c("x", "y", "z"),
    group_var = "Group",
    triangle = "lower"
  )
  expect_equal(nrow(comparison$Unadjusted$plot$data), 3)
  expect_equal(nrow(comparison$Results), 6)
  expect_true(all(c(
    "TestAvailable", "ComparisonStatus", "InferenceApproximate"
  ) %in% names(comparison$Results)))
  expect_identical(comparison$Metadata$TriangleApplied, "lower")

  phi <- PlotPhiHeatmap(
    df_Test,
    CatVars = c("binary_a", "binary_b", "binary_c"),
    triangle = "upper"
  )
  expect_equal(nrow(phi$Unadjusted$plot$data), 3)

  directional <- PlotDirectionalHeatmaps(
    df_Test,
    variables = c("x", "y", "z"),
    triangle = "lower"
  )
  expect_equal(
    nrow(directional$Unadjusted$plot$layers[[length(directional$Unadjusted$plot$layers)]]$data),
    3
  )
})
