DirectionalRectangularData <- function() {
  set.seed(711)
  df_Test <- data.frame(a = rnorm(70), b = rnorm(70), c = rnorm(70),
    g = factor(rep(c("no", "yes"), 35)), h = factor(rep(c("low", "high"), each = 35)))
  attr(df_Test$a, "label") <- "Measure A"
  df_Test
}

test_that("directional rectangular axes preserve historical orientation and statistics", {
  df_Test <- DirectionalRectangularData()
  vars_Columns <- c("g", "a")
  vars_Rows <- c("h", "b", "c")
  square <- PlotDirectionalHeatmaps(df_Test, variables = c(vars_Columns, vars_Rows))
  for (scope in c("matrix", "per_outcome", "per_predictor")) {
    result <- PlotDirectionalHeatmaps(df_Test, predictor_vars = vars_Columns,
      outcome_vars = vars_Rows, fdr_scope = scope, triangle = "upper")
    df_Result <- result$Unadjusted$data
    expect_identical(result$AxisOrder$columns, vars_Columns)
    expect_identical(result$AxisOrder$rows, vars_Rows)
    expect_equal(nrow(df_Result), 6)
    expect_equal(nrow(unique(df_Result[c("XVar", "YVar")])), 6)
    for (i in seq_len(nrow(df_Result))) {
      df_Expected <- square$Unadjusted$data
      df_Expected <- df_Expected[df_Expected$XVar == df_Result$XVar[i] &
        df_Expected$YVar == df_Result$YVar[i], ]
      expect_equal(df_Result$correlation[i], df_Expected$correlation)
      raw_p <- if (df_Result$test[i] == "pearson") df_Expected$P else df_Expected$p_value
      expect_equal(df_Result$P[i], raw_p)
    }
    for (block in unique(df_Result$test)) {
      df_Block <- df_Result[df_Result$test == block, ]
      expect_equal(df_Block$P_adj, ApplyFDRCorrection(df_Block$P, fdr_scope = scope,
        outcome_ids = ifelse(block == "Point Correlation" & df_Block$XVar == "a", df_Block$XVar, df_Block$YVar),
        predictor_ids = ifelse(block == "Point Correlation" & df_Block$XVar == "a", df_Block$YVar, df_Block$XVar)))
    }
    expect_identical(result$p, result$Unadjusted)
    expect_identical(result$p_fdr, result$FDRCorrected)
    expect_s3_class(result$Unadjusted$plot, "ggplot")
    expect_s3_class(result$FDRCorrected$plot, "ggplot")
  }
})

test_that("legacy positional and named axes use rectangular mode", {
  df_Test <- DirectionalRectangularData()
  expected <- PlotDirectionalHeatmaps(df_Test, predictor_vars = c("g", "a"), outcome_vars = c("h", "b", "c"))
  positional <- PlotDirectionalHeatmaps(df_Test, c("g", "a"), c("h", "b", "c"))
  expect_equal(positional$Unadjusted$data, expected$Unadjusted$data)
  withr::local_options(lifecycle_verbosity = "warning")
  expect_warning(expect_warning(legacy <- PlotDirectionalHeatmaps(df_Test,
    xVars = c("g", "a"), yVars = c("h", "b", "c")), "xVars"), "yVars")
  expect_equal(legacy$Unadjusted$data, expected$Unadjusted$data)
  expect_warning(ignored <- PlotDirectionalHeatmaps(df_Test, variables = "c",
    predictor_vars = c("g", "a"), outcome_vars = c("h", "b", "c")), "ignored")
  expect_equal(ignored$Unadjusted$data, expected$Unadjusted$data)
  square <- PlotDirectionalHeatmaps(df_Test, c("a", "b"), FALSE)
  expect_false(square$Relabel)
  expect_equal(nrow(square$Unadjusted$data), 4)
})

test_that("rectangles show overlapping axes with blank self-pairs and retain labels", {
  df_Test <- DirectionalRectangularData()
  result <- PlotDirectionalHeatmaps(df_Test, predictor_vars = c("a", "g", "a"),
    outcome_vars = c("a", "g"), triangle = "lower")
  df_Result <- result$Unadjusted$data
  expect_equal(nrow(df_Result), 4)
  expect_true(all(is.na(df_Result$correlation[df_Result$XVar == df_Result$YVar])))
  expect_true(all(is.na(df_Result$P[df_Result$XVar == df_Result$YVar])))
  expect_true(all(df_Result$XLabel[df_Result$XVar == "a"] == "Measure A"))
  expect_identical(result$Triangle, "full")
  names(df_Test)[4] <- "binary group"
  expect_no_error(PlotDirectionalHeatmaps(df_Test, predictor_vars = "binary group", outcome_vars = "a"))
})


test_that("rectangular exclusions and empty axes remain inspectable", {
  df_Test <- DirectionalRectangularData()
  df_Test$constant <- 1
  expect_warning(result <- PlotDirectionalHeatmaps(df_Test,
    predictor_vars = "constant", outcome_vars = "a"), "fewer than two usable values")
  expect_equal(nrow(result$Unadjusted$data), 0)
  expect_true("constant" %in% result$Excluded$Variable)
  expect_error(PlotDirectionalHeatmaps(df_Test, predictor_vars = "absent", outcome_vars = "a"), "not found")
})
