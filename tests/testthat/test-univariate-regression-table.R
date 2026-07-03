test_that("UnivariateRegressionTable fits linear and logistic outcomes", {
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = c("y", "ybin"),
    PredictorVars = c("x1", "x2")
  )

  expect_named(res, c("FormattedTable", "LargeTable", "ModelSummaries", "Metadata"))
  expect_s3_class(res$ModelSummaries$y$x1, "lm")
  expect_s3_class(res$ModelSummaries$ybin$x1, "glm")
  expect_equal(res$Metadata$Outcomes$OutcomeFamily, c("linear", "logistic"))
  expect_equal(res$Metadata$Outcomes$ReferenceLevel[[2]], "Control")
  expect_equal(res$Metadata$Outcomes$EventLevel[[2]], "Case")

  logistic_body <- res$LargeTable$tbls[[2]]$table_body
  expect_true(all(logistic_body$outcome_family == "logistic"))
  expect_true(all(logistic_body$effect_type == "Odds ratio"))
  expect_true(all(logistic_body$estimate > 0, na.rm = TRUE))
})

test_that("UnivariateRegressionTable handles binary character outcomes", {
  set.seed(910)
  df <- data.frame(
    ychar = sample(c("No", "Yes"), 70, replace = TRUE),
    x1 = rnorm(70)
  )

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = "ychar",
    PredictorVars = "x1"
  )

  expect_s3_class(res$ModelSummaries$ychar$x1, "glm")
  expect_equal(res$Metadata$Outcomes$OutcomeFamily, "logistic")
  expect_equal(res$Metadata$Outcomes$ReferenceLevel, "No")
  expect_equal(res$Metadata$Outcomes$EventLevel, "Yes")
})

test_that("UnivariateRegressionTable rejects unsupported categorical outcomes", {
  df <- data.frame(
    ymulti = factor(rep(c("A", "B", "C"), each = 10)),
    x1 = rnorm(30)
  )

  expect_error(
    UnivariateRegressionTable(
      Data = df,
      OutcomeVars = "ymulti",
      PredictorVars = "x1"
    ),
    "Multinomial regression is not yet supported"
  )
})

test_that("UnivariateRegressionTable does not standardize binary logistic outcomes", {
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = "ybin",
    PredictorVars = "x1",
    Standardize = TRUE
  )

  model_frame <- stats::model.frame(res$ModelSummaries$ybin$x1)
  expect_s3_class(res$ModelSummaries$ybin$x1, "glm")
  expect_equal(levels(model_frame$ybin), c("Control", "Case"))
})

test_that("plotForestFromTable keeps outcome labels aligned with table order", {
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = c("CBF", "E", "VTT", "PS"),
    PredictorVars = "cohort",
    Standardize = TRUE
  )
  p <- plotForestFromTable(res)
  built_plot <- ggplot2::ggplot_build(p)
  facet_layout <- built_plot$layout$layout[, c("PANEL", "PlotFacet")]
  point_data <- built_plot$data[[1]][, c("PANEL", "x")]
  point_data$PlotFacet <- as.character(facet_layout$PlotFacet[match(point_data$PANEL, facet_layout$PANEL)])

  expected <- vapply(
    names(res$FormattedTable$tbls),
    function(outcome) {
      table_body <- res$FormattedTable$tbls[[outcome]]$table_body
      table_body$estimate[!is.na(table_body$estimate)][[1]]
    },
    numeric(1)
  )
  observed <- stats::setNames(point_data$x, point_data$PlotFacet)

  expect_equal(names(observed), names(expected))
  expect_equal(unname(observed[names(expected)]), unname(expected), tolerance = 1e-8)
})

test_that("plotForestFromTable can flip outcomes and terms", {
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = c("CBF", "E", "VTT", "PS"),
    PredictorVars = "cohort",
    Standardize = TRUE
  )
  p <- plotForestFromTable(res, Flip = TRUE)
  built_plot <- ggplot2::ggplot_build(p)
  facet_layout <- built_plot$layout$layout[, c("PANEL", "PlotFacet")]
  point_data <- built_plot$data[[1]][, c("PANEL", "x", "y")]
  point_data$PlotFacet <- as.character(facet_layout$PlotFacet[match(point_data$PANEL, facet_layout$PANEL)])

  expected <- vapply(
    names(res$FormattedTable$tbls),
    function(outcome) {
      table_body <- res$FormattedTable$tbls[[outcome]]$table_body
      table_body$estimate[!is.na(table_body$estimate)][[1]]
    },
    numeric(1)
  )

  expect_equal(unique(point_data$PlotFacet), "cohort : Hiv_pos")
  expect_equal(sort(unique(point_data$y)), seq_along(expected))
  expect_equal(sort(point_data$x), sort(unname(expected)), tolerance = 1e-8)
})

test_that("plotForestFromTable validates Flip", {
  df <- data.frame(
    y = rnorm(30),
    x1 = rnorm(30)
  )
  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = "y",
    PredictorVars = "x1"
  )

  expect_error(
    plotForestFromTable(res, Flip = NA),
    "Flip must be TRUE or FALSE"
  )
})

test_that("plotForestFromTable uses outcome labels from metadata", {
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = c("cbf_m_l_100_g_min", "e"),
    PredictorVars = "cohort",
    Standardize = TRUE
  )

  p_default <- plotForestFromTable(res)
  default_layout <- ggplot2::ggplot_build(p_default)$layout$layout
  expect_equal(
    as.character(default_layout$PlotFacet),
    c("CBF", "Extraction fraction")
  )

  p_flipped <- plotForestFromTable(res, Flip = TRUE)
  flipped_data <- ggplot2::ggplot_build(p_flipped)$data[[1]]
  expect_equal(sort(unique(flipped_data$y)), c(1, 2))
})

test_that("plotForestFromTable preserves table spanner labels when metadata has raw names", {
  set.seed(915)
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

  res <- UnivariateRegressionTable(
    Data = df,
    OutcomeVars = c("cbf_m_l_100_g_min", "e"),
    PredictorVars = "cohort",
    Standardize = TRUE
  )
  res$Metadata$Outcomes$OutcomeLabel <- res$Metadata$Outcomes$Outcome

  p <- plotForestFromTable(res)
  plot_layout <- ggplot2::ggplot_build(p)$layout$layout

  expect_equal(
    as.character(plot_layout$PlotFacet),
    c("CBF", "Extraction fraction")
  )
})
