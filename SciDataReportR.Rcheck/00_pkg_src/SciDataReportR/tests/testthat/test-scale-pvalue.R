test_that("p-value scales target the requested aesthetic", {
  color_scale <- scale_color_pvalue()
  fill_scale <- scale_fill_pvalue(palette = "viridis")

  expect_identical(color_scale$aesthetics, "colour")
  expect_identical(fill_scale$aesthetics, "fill")
  expect_identical(color_scale$name, "P-value")
  expect_identical(fill_scale$name, "P-value")
})


test_that("p-value scales use a taller default colorbar", {
  default_scale <- scale_color_pvalue()
  standard_scale <- scale_color_pvalue(guide = "colourbar")

  expect_s3_class(default_scale$guide, "GuideColourbar")
  expect_equal(
    default_scale$guide$params$theme$legend.key.height,
    grid::unit(70, "mm")
  )
  expect_identical(standard_scale$guide, "colourbar")
  expect_equal(
    default_scale$breaks,
    c(1, 0.1, 0.05, 0.01, 0.001, 1e-5, 1e-8)
  )
})


test_that("p-value scales retain raw guide configuration", {
  label_function <- function(x) paste0("p=", x)
  scale <- scale_color_pvalue(
    direction = 1,
    limits = c(1e-6, 0.5),
    breaks = c(0.5, 0.05, 0.001, 1e-6),
    labels = label_function
  )

  expect_equal(
    sort(scale$trans$inverse(scale$limits)),
    c(1e-6, 0.5)
  )
  expect_equal(scale$breaks, c(0.5, 0.05, 0.001, 1e-6))
  expect_equal(
    scale$get_labels(scale$transform(c(0.5, 0.05, 0.001, 1e-6))),
    label_function(c(0.5, 0.05, 0.001, 1e-6))
  )
})


test_that("p-value scale arguments fail fast", {
  expect_error(
    scale_color_pvalue(palette = "magma"),
    "should be one of"
  )
  expect_error(scale_color_pvalue(direction = 0), "either -1 or 1")
  expect_error(scale_fill_pvalue(direction = c(-1, 1)), "either -1 or 1")
  expect_error(scale_fill_pvalue(limits = c(0.05)), "two finite")
  expect_error(scale_fill_pvalue(limits = c(0, 1)), "positive values")
  expect_error(scale_fill_pvalue(limits = c(1, 0.05)), "in increasing order")
  expect_error(scale_fill_pvalue(limits = c(0.05, Inf)), "two finite")
})


test_that("threshold warping preserves raw p-value evidence ordering", {
  p_values <- c(1, 0.1, 0.05, 0.01, 0.001, 1e-5, 1e-8, 0)
  warped <- SciDataReportR:::.warp_pvalues(p_values)

  expect_equal(warped[1:7], c(0, 0.16, 0.32, 0.56, 0.75, 0.90, 1))
  expect_true(all(diff(warped) >= 0))
  expect_equal(warped[8], 1)
})


test_that("threshold warping is invertible", {
  p_values <- c(1, 0.7, 0.1, 0.075, 0.05, 0.02, 0.01, 0.001, 1e-5, 1e-8)
  warped <- SciDataReportR:::.warp_pvalues(p_values)

  expect_equal(
    SciDataReportR:::.unwarp_pvalues(warped),
    p_values,
    tolerance = 1e-12
  )
  expect_true(is.na(SciDataReportR:::.warp_pvalues(NA_real_)))
  expect_true(is.na(SciDataReportR:::.unwarp_pvalues(NA_real_)))
})


test_that("non-significant colors are progressively desaturated", {
  positions <- SciDataReportR:::.warp_pvalues(c(1, 0.5, 0.05))
  colors <- SciDataReportR:::.desaturate_pvalue_colors(
    colors = rep("#FF0000", 3),
    positions = positions
  )
  saturation <- grDevices::rgb2hsv(grDevices::col2rgb(colors))["s", ]

  expect_lt(saturation[1], saturation[2])
  expect_lt(saturation[2], saturation[3])
  expect_equal(saturation[3], 1, tolerance = 0.01)
})


test_that("palette direction and palette choice affect mapped colors", {
  p_values <- c(1, 0.05, 0.01, 0.001, 1e-8)
  inferno <- scale_color_pvalue()
  reversed <- scale_color_pvalue(direction = 1)
  viridis <- scale_color_pvalue(palette = "viridis")
  transformed <- inferno$transform(p_values)

  expect_false(identical(inferno$map(transformed), reversed$map(transformed)))
  expect_false(identical(inferno$map(transformed), viridis$map(transformed)))
  expect_identical(inferno$map(inferno$transform(NA_real_)), "grey80")
  expect_identical(
    inferno$map(inferno$transform(0)),
    inferno$map(inferno$transform(1e-8))
  )
})


test_that("trained color and fill guides use warped raw p-value breaks", {
  df_Guide <- data.frame(
    x = seq_len(7),
    y = seq_len(7),
    PValue = c(1, 0.1, 0.05, 0.01, 0.001, 1e-5, 1e-8)
  )
  expected_positions <- c(0, 0.16, 0.32, 0.56, 0.75, 0.90, 1)
  expected_labels <- c("1", "0.1", "0.05", "0.01", "0.001", "1e-05", "1e-08")

  color_build <- ggplot2::ggplot_build(
    ggplot2::ggplot(
      df_Guide,
      ggplot2::aes(x = x, y = y, color = PValue)
    ) +
      ggplot2::geom_point() +
      scale_color_pvalue()
  )
  fill_build <- ggplot2::ggplot_build(
    ggplot2::ggplot(
      df_Guide,
      ggplot2::aes(x = x, y = y, fill = PValue)
    ) +
      ggplot2::geom_tile() +
      scale_fill_pvalue()
  )
  color_key <- color_build$plot$guides$params[[1]]$key
  fill_key <- fill_build$plot$guides$params[[1]]$key

  expect_equal(color_key$.value, expected_positions, tolerance = 0.003)
  expect_equal(fill_key$.value, expected_positions, tolerance = 0.003)
  expect_identical(color_key$.label, expected_labels)
  expect_identical(fill_key$.label, expected_labels)
})


test_that("ANOVA relationship matrices use the appropriate continuous scale", {
  set.seed(42)
  df_Anova <- data.frame(
    Group = factor(rep(c("Control", "Case"), each = 20)),
    Outcome = c(stats::rnorm(20), stats::rnorm(20, mean = 0.8))
  )

  effect_result <- PlotAnovaRelationshipsMatrix(
    df_Anova,
    CatVars = "Group",
    ContVars = "Outcome",
    Relabel = FALSE,
    Parametric = TRUE
  )
  pvalue_result <- PlotAnovaRelationshipsMatrix(
    df_Anova,
    CatVars = "Group",
    ContVars = "Outcome",
    Relabel = FALSE,
    Parametric = FALSE
  )

  effect_scale <- effect_result$Unadjusted$plot$scales$get_scales("colour")
  raw_scale <- pvalue_result$Unadjusted$plot$scales$get_scales("colour")
  adjusted_scale <- pvalue_result$FDRCorrected$plot$scales$get_scales("colour")
  raw_mapping <- rlang::get_expr(
    pvalue_result$Unadjusted$plot$layers[[1]]$mapping$colour
  )
  adjusted_mapping <- rlang::get_expr(
    pvalue_result$FDRCorrected$plot$layers[[1]]$mapping$colour
  )

  expect_identical(effect_scale$name, "Effect Size")
  expect_identical(raw_scale$name, "P-value")
  expect_identical(adjusted_scale$name, "FDR-adjusted P-value")
  expect_identical(raw_mapping, quote(.data[["p"]]))
  expect_identical(adjusted_mapping, quote(.data[["p.adj"]]))
})
