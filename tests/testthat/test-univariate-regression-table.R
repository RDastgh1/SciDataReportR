test_that("MakeUnivariateRegressionTable fits linear and logistic outcomes", {
  set.seed(909)
  df <- data.frame(
    y = rnorm(90),
    ybin = factor(
      sample(c("Control", "Case"), 90, replace = TRUE),
      levels = c("Control", "Case")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )
  attr(df$x1, "label") <- "Marker one"

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("y", "ybin"),
    predictor_vars = c("x1", "x2"),
    ReturnModels = TRUE
  )

  expect_named(res, c("FormattedTable", "LargeTable", "Results", "ModelSummaries", "Metadata"))
  expect_s3_class(res$FormattedTable, "gt_tbl")
  expect_s3_class(res$LargeTable, "gt_tbl")
  expect_s3_class(res$ModelSummaries$y$x1, "lm")
  expect_s3_class(res$ModelSummaries$ybin$x1, "glm")
  expect_equal(res$Metadata$Outcomes$OutcomeFamily, c("linear", "logistic"))
  expect_equal(res$Metadata$Outcomes$ReferenceLevel[[2]], "Control")
  expect_equal(res$Metadata$Outcomes$EventLevel[[2]], "Case")

  logistic_results <- res$Results[res$Results$OutcomeFamily == "logistic", ]
  expect_true(all(logistic_results$EffectType == "Odds ratio"))
  expect_true(all(logistic_results$Estimate > 0, na.rm = TRUE))
})

test_that("MakeUnivariateRegressionTable skips model objects by default", {
  set.seed(908)
  df <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = "x1"
  )

  expect_null(res$ModelSummaries)
  expect_false(res$Metadata$AnalysisSettings$ReturnModels)
})

test_that("MakeUnivariateRegressionTable returns a tidy Results dataframe", {
  set.seed(917)
  df <- data.frame(
    y = rnorm(90),
    ybin = factor(
      sample(c("Control", "Case"), 90, replace = TRUE),
      levels = c("Control", "Case")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )
  attr(df$x1, "label") <- "Marker one"
  attr(df$y, "label") <- "Outcome one"

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("y", "ybin"),
    predictor_vars = c("x1", "x2"),
    ReturnModels = TRUE
  )

  results <- res$Results
  expect_s3_class(results, "data.frame")
  expect_true(all(
    c(
      "Outcome", "OutcomeLabel", "OutcomeFamily", "EffectType", "Predictor",
      "PredictorLabel", "Term", "Level", "TermLabel", "N", "Estimate",
      "StdError", "ConfLow", "ConfHigh", "PValue", "Significant",
      "ReferenceValue"
    ) %in% names(results)
  ))
  # one row per outcome-predictor pair for continuous predictors
  expect_equal(nrow(results), 4)
  expect_equal(results$OutcomeLabel[results$Outcome == "y"][1], "Outcome one")
  expect_equal(results$TermLabel[results$Outcome == "y"], c("Marker one", "x2"))

  # linear estimates match the fitted models
  expect_equal(
    results$Estimate[results$Outcome == "y" & results$Predictor == "x1"],
    unname(stats::coef(res$ModelSummaries$y$x1)[["x1"]]),
    tolerance = 1e-8
  )
  # logistic estimates are exponentiated odds ratios
  expect_equal(
    results$Estimate[results$Outcome == "ybin" & results$Predictor == "x1"],
    exp(unname(stats::coef(res$ModelSummaries$ybin$x1)[["x1"]])),
    tolerance = 1e-8
  )
  expect_true(all(results$ConfLow < results$ConfHigh))
  expect_equal(results$ReferenceValue, c(0, 0, 1, 1))
  expect_equal(results$EffectType, c("Estimate", "Estimate", "Odds ratio", "Odds ratio"))
  expect_equal(results$Significant, results$PValue < 0.05)
})

test_that("UnivariateRegressionTable is a backwards-compatible synonym", {
  set.seed(918)
  df <- data.frame(
    y = rnorm(60),
    x1 = rnorm(60)
  )

  lifecycle::expect_deprecated(
    res_old <- UnivariateRegressionTable(
      data = df,
      outcome_vars = "y",
      predictor_vars = "x1",
      ReturnModels = TRUE
    )
  )
  res_new <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = "x1",
    ReturnModels = TRUE
  )

  expect_named(res_old, names(res_new))
  expect_equal(res_old$Results, res_new$Results)
  expect_equal(
    stats::coef(res_old$ModelSummaries$y$x1),
    stats::coef(res_new$ModelSummaries$y$x1)
  )
})

test_that("MakeUnivariateRegressionTable handles binary character outcomes", {
  set.seed(910)
  df <- data.frame(
    ychar = sample(c("No", "Yes"), 70, replace = TRUE),
    x1 = rnorm(70)
  )

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "ychar",
    predictor_vars = "x1",
    ReturnModels = TRUE
  )

  expect_s3_class(res$ModelSummaries$ychar$x1, "glm")
  expect_equal(res$Metadata$Outcomes$OutcomeFamily, "logistic")
  expect_equal(res$Metadata$Outcomes$ReferenceLevel, "No")
  expect_equal(res$Metadata$Outcomes$EventLevel, "Yes")
})

test_that("MakeUnivariateRegressionTable rejects unsupported categorical outcomes", {
  df <- data.frame(
    ymulti = factor(rep(c("A", "B", "C"), each = 10)),
    x1 = rnorm(30)
  )

  expect_error(
    MakeUnivariateRegressionTable(
      data = df,
      outcome_vars = "ymulti",
      predictor_vars = "x1"
    ),
    "Multinomial regression is not yet supported"
  )
})

test_that("MakeUnivariateRegressionTable does not standardize binary logistic outcomes", {
  set.seed(911)
  x1 <- rnorm(80)
  probability <- stats::plogis(-0.3 + 0.8 * x1)
  df <- data.frame(
    ybin = factor(
      ifelse(stats::runif(80) < probability, "Case", "Control"),
      levels = c("Control", "Case")
    ),
    x1 = x1
  )

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "ybin",
    predictor_vars = "x1",
    Standardize = TRUE,
    ReturnModels = TRUE
  )

  model_frame <- stats::model.frame(res$ModelSummaries$ybin$x1)
  expect_s3_class(res$ModelSummaries$ybin$x1, "glm")
  expect_equal(levels(model_frame$ybin), c("Control", "Case"))
})

test_that("PlotForestFromTable keeps outcome labels aligned with table order", {
  set.seed(912)
  df <- data.frame(
    CBF = rnorm(80, mean = 0.1),
    E = rnorm(80, mean = -0.1),
    VTT = rnorm(80, mean = 0.2),
    PS = rnorm(80, mean = -0.2),
    cohort = factor(
      rep(c("control", "Hiv_pos"), each = 40),
      levels = c("control", "Hiv_pos")
    )
  )
  df$CBF[df$cohort == "Hiv_pos"] <- df$CBF[df$cohort == "Hiv_pos"] + 0.10
  df$E[df$cohort == "Hiv_pos"] <- df$E[df$cohort == "Hiv_pos"] + 0.25
  df$VTT[df$cohort == "Hiv_pos"] <- df$VTT[df$cohort == "Hiv_pos"] - 0.35
  df$PS[df$cohort == "Hiv_pos"] <- df$PS[df$cohort == "Hiv_pos"] + 0.50

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("CBF", "E", "VTT", "PS"),
    predictor_vars = "cohort",
    Standardize = TRUE
  )
  p <- PlotForestFromTable(res)
  built_plot <- ggplot2::ggplot_build(p)
  facet_layout <- built_plot$layout$layout[, c("PANEL", "PlotFacet")]
  point_data <- built_plot$data[[1]][, c("PANEL", "x")]
  point_data$PlotFacet <- as.character(facet_layout$PlotFacet[match(point_data$PANEL, facet_layout$PANEL)])

  expected <- stats::setNames(res$Results$Estimate, res$Results$OutcomeLabel)
  observed <- stats::setNames(point_data$x, point_data$PlotFacet)

  expect_equal(names(observed), names(expected))
  expect_equal(unname(observed[names(expected)]), unname(expected), tolerance = 1e-8)
})

test_that("PlotForestFromTable can flip outcomes and terms", {
  set.seed(913)
  df <- data.frame(
    CBF = rnorm(80),
    E = rnorm(80),
    VTT = rnorm(80),
    PS = rnorm(80),
    cohort = factor(
      rep(c("control", "Hiv_pos"), each = 40),
      levels = c("control", "Hiv_pos")
    )
  )
  df$CBF[df$cohort == "Hiv_pos"] <- df$CBF[df$cohort == "Hiv_pos"] + 0.10
  df$E[df$cohort == "Hiv_pos"] <- df$E[df$cohort == "Hiv_pos"] + 0.25
  df$VTT[df$cohort == "Hiv_pos"] <- df$VTT[df$cohort == "Hiv_pos"] - 0.35
  df$PS[df$cohort == "Hiv_pos"] <- df$PS[df$cohort == "Hiv_pos"] + 0.50

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("CBF", "E", "VTT", "PS"),
    predictor_vars = "cohort",
    Standardize = TRUE
  )
  p <- PlotForestFromTable(res, Flip = TRUE)
  built_plot <- ggplot2::ggplot_build(p)
  facet_layout <- built_plot$layout$layout[, c("PANEL", "PlotFacet")]
  point_data <- built_plot$data[[1]][, c("PANEL", "x", "y")]
  point_data$PlotFacet <- as.character(facet_layout$PlotFacet[match(point_data$PANEL, facet_layout$PANEL)])

  expected <- stats::setNames(res$Results$Estimate, res$Results$OutcomeLabel)

  expect_equal(unique(point_data$PlotFacet), "cohort : Hiv_pos")
  expect_equal(sort(unique(point_data$y)), seq_along(expected))
  expect_equal(sort(point_data$x), sort(unname(expected)), tolerance = 1e-8)
})

test_that("PlotForestFromTable validates Flip", {
  df <- data.frame(
    y = rnorm(30),
    x1 = rnorm(30)
  )
  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = "x1"
  )

  expect_error(
    PlotForestFromTable(res, Flip = NA),
    "Flip must be TRUE or FALSE"
  )
})

test_that("PlotForestFromTable accepts a filtered Results dataframe", {
  set.seed(919)
  df <- data.frame(
    y1 = rnorm(80),
    y2 = rnorm(80),
    x1 = rnorm(80),
    x2 = rnorm(80)
  )
  df$y1 <- df$y1 + 0.9 * df$x1

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("y1", "y2"),
    predictor_vars = c("x1", "x2")
  )
  filtered <- res$Results[res$Results$Outcome == "y1", ]

  p <- PlotForestFromTable(filtered)
  built_plot <- ggplot2::ggplot_build(p)
  plot_layout <- built_plot$layout$layout
  point_data <- built_plot$data[[1]]

  expect_equal(as.character(plot_layout$PlotFacet), "y1")
  expect_equal(
    sort(point_data$x),
    sort(filtered$Estimate),
    tolerance = 1e-8
  )
})

test_that("PlotForestFromTable rejects dataframes missing required columns", {
  expect_error(
    PlotForestFromTable(data.frame(Estimate = 1)),
    "missing required columns"
  )
})

test_that("plotForestFromTable is a backwards-compatible synonym", {
  set.seed(921)
  df <- data.frame(
    y = rnorm(60),
    x1 = rnorm(60)
  )
  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = "x1"
  )

  lifecycle::expect_deprecated(p_old <- plotForestFromTable(res))
  p_new <- PlotForestFromTable(res)

  expect_equal(
    ggplot2::ggplot_build(p_old)$data[[1]]$x,
    ggplot2::ggplot_build(p_new)$data[[1]]$x,
    tolerance = 1e-8
  )
})

test_that("PlotForestFromTable uses outcome labels from metadata", {
  set.seed(914)
  df <- data.frame(
    cbf_m_l_100_g_min = rnorm(60),
    e = rnorm(60),
    cohort = factor(
      rep(c("control", "Hiv_pos"), each = 30),
      levels = c("control", "Hiv_pos")
    )
  )
  attr(df$cbf_m_l_100_g_min, "label") <- "CBF"
  attr(df$e, "label") <- "Extraction fraction"

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("cbf_m_l_100_g_min", "e"),
    predictor_vars = "cohort",
    Standardize = TRUE
  )

  p_default <- PlotForestFromTable(res)
  default_layout <- ggplot2::ggplot_build(p_default)$layout$layout
  expect_equal(
    as.character(default_layout$PlotFacet),
    c("CBF", "Extraction fraction")
  )

  p_flipped <- PlotForestFromTable(res, Flip = TRUE)
  flipped_data <- ggplot2::ggplot_build(p_flipped)$data[[1]]
  expect_equal(sort(unique(flipped_data$y)), c(1, 2))
})

test_that("PlotForestFromTable uses odds ratio null line for logistic models", {
  set.seed(916)
  df <- data.frame(
    ybin = factor(
      sample(c("Typical", "High"), 90, replace = TRUE),
      levels = c("Typical", "High")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )

  res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "ybin",
    predictor_vars = c("x1", "x2"),
    Standardize = TRUE
  )
  p <- PlotForestFromTable(res)
  built_plot <- ggplot2::ggplot_build(p)
  vline_data <- built_plot$data[[3]]

  expect_equal(unique(vline_data$xintercept), 1)
  expect_equal(p$labels$x, "Odds ratio")
})

test_that("MakeUnivariateRegressionTable still accepts deprecated argument names", {
  set.seed(922)
  df <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  deprecation_warnings <- testthat::capture_warnings(
    res <- MakeUnivariateRegressionTable(
      Data = df,
      OutcomeVars = "y",
      PredictorVars = "x1",
      ReturnModels = TRUE
    )
  )
  expect_length(deprecation_warnings, 3)
  expect_match(deprecation_warnings, "deprecated", all = TRUE)
  expect_s3_class(res$ModelSummaries$y$x1, "lm")
  expect_s3_class(res$Results, "data.frame")
})
