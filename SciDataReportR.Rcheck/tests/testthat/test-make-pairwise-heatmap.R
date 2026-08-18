test_that("MakePairwiseHeatmap requires a Referent", {
  df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 5)),
    y = rnorm(10)
  )

  expect_error(
    MakePairwiseHeatmap(df, group_var = "Group", variables = "y"),
    "Referent"
  )
})

test_that("MakePairwiseHeatmap returns group-minus-referent Z-score contrasts", {
  skip_if_not_installed("emmeans")

  df <- data.frame(
    Group = factor(
      rep(c("Control", "A", "B"), each = 12),
      levels = c("Control", "A", "B")
    ),
    y1 = c(seq(-1, 1, length.out = 12),
           seq(1.5, 3.5, length.out = 12),
           seq(-3.5, -1.5, length.out = 12)),
    y2 = c(seq(0, 2, length.out = 12),
           seq(2, 4, length.out = 12),
           seq(4, 6, length.out = 12))
  )

  df$y1 <- sjlabelled::set_label(df$y1, "Marker one")

  res <- MakePairwiseHeatmap(
    data = df,
    group_var = "Group",
    variables = c("y1", "y2"),
    Referent = "Control",
    Parametric = TRUE,
    adjust_scope = "matrix",
    p_adjust_method = "fdr"
  )

  expect_s3_class(res, "SciDataReportRPairwiseHeatmap")
  expect_s3_class(res$Plot, "ggplot")
  expect_equal(unique(res$Results$Referent), "Control")
  expect_false("Control" %in% as.character(res$Results$ColumnLabel))
  expect_true(any(res$Results$VariableLabel == "Marker one"))
  expect_equal(levels(res$Results$ColumnLabel), c("A", "B"))
  expect_equal(levels(res$Results$RowLabel), rev(c("Marker one", "y2")))
  expect_equal(unique(res$Results$ScoreType), "ZScore")
  expect_equal(
    unique(res$Results$Test),
    "Linear model + emmeans referent contrast"
  )
  expect_equal(unique(res$Results$Contrast), "Group - Referent")
  expect_true(all(grepl("Referent-scaled Z-score ~ Group", res$Results$ModelFormula, fixed = TRUE)))
  expect_true(all(grepl("FDR correction across all displayed cells", res$Results$Adjustment, fixed = TRUE)))
  expect_true("ModelFormulaInternal" %in% names(res$Results))

  z_ref <- CreateZScoreObject(
    df[df$Group == "Control", , drop = FALSE],
    variables = c("y1", "y2"),
    names_prefix = ".PairwiseHeatmap_Z_",
    RetainLabels = FALSE
  )
  z_obj <- ProjectZScore(
    df,
    variables = c("y1", "y2"),
    parameters = z_ref,
    ParameterInputType = "ZScoreObj",
    names_prefix = ".PairwiseHeatmap_Z_",
    RetainLabels = FALSE
  )
  fit <- stats::lm(`.PairwiseHeatmap_Z_y1` ~ Group, data = z_obj$DataWithZ)
  emm <- emmeans::emmeans(fit, ~ Group)
  ctr <- as.data.frame(summary(
    emmeans::contrast(
      emm,
      method = list(A = c(-1, 1, 0), B = c(-1, 0, 1))
    ),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ))

  got_a <- res$Results$EstimatedMeanDifference[
    res$Results$Variable == "y1" & res$Results$Group == "A"
  ]
  got_b <- res$Results$EstimatedMeanDifference[
    res$Results$Variable == "y1" & res$Results$Group == "B"
  ]

  expect_equal(got_a, ctr$estimate[ctr$contrast == "A"], tolerance = 1e-8)
  expect_equal(got_b, ctr$estimate[ctr$contrast == "B"], tolerance = 1e-8)
  expect_gt(got_a, 0)
  expect_lt(got_b, 0)
})

test_that("MakePairwiseHeatmap uses M-scores for non-parametric mode", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("sandwich")

  df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 8), levels = c("Control", "A")),
    y = c(1:8, 6:13)
  )

  res <- MakePairwiseHeatmap(
    data = df,
    group_var = "Group",
    variables = "y",
    Referent = "Control",
    Parametric = FALSE,
    p_adjust_method = "none"
  )

  m_obj <- CreateMScoreObject(
    df[df$Group == "Control", , drop = FALSE],
    variables = "y",
    names_prefix = ".PairwiseHeatmap_M_",
    RetainLabels = FALSE
  )

  expect_equal(unique(res$Results$ScoreType), "MScore")
  expect_equal(
    unique(res$Results$Test),
    "Robust linear model (HC3) + emmeans referent contrast"
  )
  expect_equal(
    res$ScalingParameters$Median,
    m_obj$Parameters$Median,
    tolerance = 1e-10
  )
  expect_equal(
    res$ScalingParameters$MAD,
    m_obj$Parameters$MAD,
    tolerance = 1e-10
  )
})

test_that("MakePairwiseHeatmap adjustment scopes match p.adjust expectations", {
  skip_if_not_installed("emmeans")

  set.seed(101)
  df <- data.frame(
    Group = factor(rep(c("Control", "A", "B"), each = 10), levels = c("Control", "A", "B")),
    y1 = c(rnorm(10, 0), rnorm(10, 1.5), rnorm(10, -1.5)),
    y2 = c(rnorm(10, 0), rnorm(10, 1.0), rnorm(10, -1.0)),
    y3 = c(rnorm(10, 0), rnorm(10, 0.5), rnorm(10, -0.5))
  )

  res_matrix <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2", "y3"),
    Referent = "Control",
    adjust_scope = "matrix",
    p_adjust_method = "holm"
  )

  expect_equal(
    res_matrix$Results$AdjustedPValue,
    stats::p.adjust(res_matrix$Results$PValue, method = "holm"),
    tolerance = 1e-12
  )

  res_group <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2", "y3"),
    Referent = "Control",
    adjust_scope = "per_group",
    p_adjust_method = "bonferroni"
  )

  expected <- rep(NA_real_, nrow(res_group$Results))
  for (idx in split(seq_len(nrow(res_group$Results)), res_group$Results$Group)) {
    expected[idx] <- stats::p.adjust(res_group$Results$PValue[idx], method = "bonferroni")
  }

  expect_equal(res_group$Results$AdjustedPValue, expected, tolerance = 1e-12)
})

test_that("MakePairwiseHeatmap plot metadata and annotation layers are stable", {
  skip_if_not_installed("emmeans")

  df <- data.frame(
    Group = factor(rep(c("Control", "A"), each = 15), levels = c("Control", "A")),
    y1 = c(seq(-0.2, 0.2, length.out = 15), seq(4.8, 5.2, length.out = 15)),
    y2 = c(seq(-0.2, 0.2, length.out = 15), seq(-5.2, -4.8, length.out = 15))
  )

  res <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2"),
    Referent = "Control",
    star_p = "adjusted",
    adjusted_outline = TRUE
  )

  expect_equal(abs(res$Settings$fill_limits_used[1]), res$Settings$fill_limits_used[2])
  expect_equal(res$Settings$ScalingReference, "Referent")
  expect_equal(res$Settings$x_axis_text_angle, 0)
  expect_null(res$Plot$labels$caption)
  expect_match(res$Plot$scales$scales[[1]]$name, "Referent-scaled")
  expect_equal(
    unname(res$Settings$colors),
    c("#52BCA3FF", "white", "#E58606FF")
  )
  expect_true(any(res$Results$SignificanceLabel != ""))
  expect_true(any(res$Results$IsAdjustedSignificant))

  geoms <- vapply(res$Plot$layers, function(layer) class(layer$geom)[1], character(1))
  expect_gte(sum(geoms == "GeomTile"), 2)
  expect_true(any(geoms == "GeomText"))

  res_limited <- MakePairwiseHeatmap(
    df,
    group_var = "Group",
    variables = c("y1", "y2"),
    Referent = "Control",
    fill_limits = c(-2, 2)
  )
  expect_equal(res_limited$Settings$fill_limits_used, c(-2, 2))
})
