SigInteractionFixture <- function(levels = c("A", "B"), continuous = FALSE) {
  set.seed(611)
  n <- 180 * length(levels)
  df_Test <- data.frame(x = rnorm(n), z = rnorm(n),
    group = factor(rep(levels, each = 180), levels = levels))
  df_Test$moderator <- if (continuous) rnorm(n, 2, 0.8) else df_Test$group
  df_Test$y <- 0.3 * df_Test$x + df_Test$x *
    (if (continuous) df_Test$moderator else as.numeric(df_Test$group)) * 1.5 +
    0.7 * df_Test$z + rnorm(n, sd = 0.3)
  df_Test$null <- rnorm(n)
  df_Test
}

test_that("significant selection reads Results and Metadata only", {
  df_Test <- SigInteractionFixture()
  result <- PlotInteractionEffectsMatrix(df_Test, "moderator", c("y", "null"), "x", "z")
  rows <- result$Results$Term == "Interaction"
  # Deterministic selection boundaries, independent of random null p-values.
  result$Results$PValue[rows & result$Results$Outcome == "null"] <- 0.8
  result$Results$PValueFDR[rows] <- 0.8
  result$Unadjusted$plot <- NULL
  plots <- plotSigInteractions(df_Test, result)
  expect_named(plots, "y__x")
  expect_s3_class(plots[[1]], "ggplot")
  expect_equal(attr(plots, "Results")$Outcome, "y")
  expect_message(empty <- plotSigInteractions(df_Test, result, PVar = "PValueFDR"), "No interactions")
  expect_length(empty, 0)
  expect_error(plotSigInteractions(df_Test, list()), "Results and Metadata")
  result$Results$PValue[rows] <- 0.001
  expect_warning(limited <- plotSigInteractions(df_Test, result, max_plots = 1), "max_plots")
  expect_length(limited, 1)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- plotSigInteractions(df_Test, result, ncol = 2)
    expect_s3_class(attr(combined, "Layout"), "patchwork")
  }
})

test_that("categorical adjusted lines and ribbons match matrix slopes", {
  for (levels in list(c("A", "B"), c("A", "B", "C"))) {
    df_Test <- SigInteractionFixture(levels)
    df_Test$y <- df_Test$y + df_Test$z * as.numeric(df_Test$group)
    result <- PlotInteractionEffectsMatrix(df_Test, "moderator", "y", "x", "z",
                                          covariate_interactions = TRUE)
    plots <- plotSigInteractions(df_Test, result)
    plot <- plots[[1]]
    df_Lines <- attr(plot, "LineData")
    df_Slopes <- result$Results[result$Results$Term != "Interaction", ]
    for (level in levels) {
      df_Line <- df_Lines[df_Lines$Group == level, ]
      slope <- diff(range(df_Line$fit)) / diff(range(df_Line$predictor))
      expected <- df_Slopes$EstimateStd[df_Slopes$Level == level]
      expect_equal(slope, abs(expected))
      expect_true(all(df_Line$lower <= df_Line$fit & df_Line$fit <= df_Line$upper))
    }
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    expect_s3_class(plot$Results, "tbl_df")
    expect_type(plot$Metadata, "list")
    expect_equal(plot$Metadata$Display$SD[["Predictor"]], sd(df_Test$x))
    raw <- plotSigInteractions(df_Test, result, adjust = FALSE, standardize = FALSE)[[1]]
    expect_equal(raw$data$predictor, df_Test$x)
    expect_equal(raw$data$outcome, df_Test$y)
    expect_equal(attr(raw, "interaction_p"), result$Results$PValue[result$Results$Term == "Interaction"])
  }
})

test_that("continuous moderators use endpoint slopes and continuous point colours", {
  df_Test <- SigInteractionFixture(continuous = TRUE)
  result <- PlotInteractionEffectsMatrix(df_Test, "moderator", "y", "x", "z")
  plot <- plotSigInteractions(df_Test, result)[[1]]
  df_Lines <- attr(plot, "LineData")
  expect_setequal(df_Lines$Group, c("Mean - 1 SD", "Mean + 1 SD"))
  df_Slopes <- result$Results[result$Results$Term != "Interaction", ]
  for (level in df_Slopes$Level) {
    df_Line <- df_Lines[df_Lines$Group == level, ]
    expect_equal(diff(df_Line$fit[1:2]) / diff(df_Line$predictor[1:2]),
                 df_Slopes$EstimateStd[df_Slopes$Level == level])
  }
  expect_true(is.numeric(plot$data$moderator))
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  single <- PlotInteractionEffectsContinuous(df_Test, "moderator", "y", "x", "z")
  expect_s3_class(single, "ggplot")
  expect_equal(nrow(single$Results), 3)
  expect_true(is.list(single$Metadata))
})

test_that("interaction captions use displayed units outside the panel", {
  df_Test <- SigInteractionFixture()
  attr(df_Test$moderator, "label") <- "Moderator label"
  source <- PlotInteractionEffectsMatrix(df_Test, "moderator", "y", "x", "z")
  row <- source$Results[source$Results$Term == "Interaction", ]
  for (scaled in c(TRUE, FALSE)) {
    plot <- plotSigInteractions(df_Test, source, standardize = scaled)[[1]]
    column <- if (scaled) "EstimateStd" else "Estimate"
    expect_match(plot$labels$caption, as.character(signif(row[[column]], 3)), fixed = TRUE)
    expect_false(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
    expect_equal(plot$scales$get_scales("colour")$name, "Moderator label")
    expect_no_error(ggplot2::ggplot_build(plot))
    slopes <- plot$Results[plot$Results$Term != "Interaction", ]
    expect_equal(diff(slopes[[column]]), row[[column]], tolerance = 1e-8)
  }
})

test_that("empty interaction notes render in narrow and wide panels", {
  df_Test <- SigInteractionFixture()
  result <- PlotInteractionEffectsMatrix(df_Test, "moderator", "null", "x", "z", Pthresh = 1e-10)
  plot <- result$FDRCorrected$plot
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_equal(plot$coordinates$clip, "off")
  expect_false(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLabel"), logical(1))))
  for (dimensions in list(c(5, 6), c(9, 5))) {
    path <- tempfile(fileext = ".png")
    grDevices::png(path, width = dimensions[1], height = dimensions[2], units = "in", res = 100)
    expect_no_error(print(plot))
    grDevices::dev.off()
    unlink(path)
  }
})
