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
