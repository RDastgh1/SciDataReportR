MakeVolcanoTestData <- function() {
  set.seed(2468)
  n <- 80
  outcome <- as.numeric(scale(stats::rnorm(n)))

  data.frame(
    outcome = outcome,
    strong = outcome + stats::rnorm(n, sd = 0.25),
    moderate = 0.25 * outcome + stats::rnorm(n),
    marginal = 0.22 * outcome + stats::rnorm(n),
    null1 = stats::rnorm(n),
    null2 = stats::rnorm(n),
    null3 = stats::rnorm(n),
    null4 = stats::rnorm(n),
    null5 = stats::rnorm(n)
  )
}

VolcanoTestVariables <- function() {
  c("strong", "moderate", "marginal", "null1", "null2", "null3", "null4", "null5")
}

test_that("PlotVolcanoEffects default coloring remains unchanged", {
  df_test <- MakeVolcanoTestData()
  vars_test <- VolcanoTestVariables()

  omitted <- PlotVolcanoEffects(
    df_test,
    vars_test,
    "outcome",
    InteractiveLabels = FALSE
  )
  explicit_null <- PlotVolcanoEffects(
    df_test,
    vars_test,
    "outcome",
    InteractiveLabels = FALSE,
    ColorBy = NULL
  )

  expect_named(omitted, c("RawPPlot", "FDRPlot", "ResultsTable"))
  expect_identical(omitted$ResultsTable, explicit_null$ResultsTable)
  expect_false("ColorCategory" %in% names(omitted$ResultsTable))
  expect_true(all(is.na(omitted$ResultsTable$PlotLabel)))
  expect_identical(
    rlang::get_expr(omitted$RawPPlot$mapping$colour),
    rlang::get_expr(explicit_null$RawPPlot$mapping$colour)
  )
  expect_identical(
    rlang::get_expr(omitted$FDRPlot$mapping$colour),
    rlang::get_expr(explicit_null$FDRPlot$mapping$colour)
  )
  expect_equal(
    omitted$ResultsTable$FDR,
    stats::p.adjust(omitted$ResultsTable$PValue, method = "fdr")
  )
})

test_that("PlotVolcanoEffects colors both plots from a category data frame", {
  df_test <- MakeVolcanoTestData()
  vars_test <- VolcanoTestVariables()
  df_categories <- data.frame(
    Variable = c("strong", "moderate", "marginal", "not_tested"),
    Category = c("Immune", "Immune", "Metabolic", "Ignored")
  )

  result <- PlotVolcanoEffects(
    df_test,
    vars_test,
    "outcome",
    InteractiveLabels = FALSE,
    ColorBy = df_categories
  )

  expected_categories <- c(
    "Immune", "Immune", "Metabolic",
    rep("Unmapped", length(vars_test) - 3)
  )

  expect_identical(result$ResultsTable$ColorCategory, expected_categories)
  expect_identical(result$RawPPlot$data$ColorCategory, expected_categories)
  expect_identical(result$FDRPlot$data$ColorCategory, expected_categories)
  expect_identical(
    rlang::get_expr(result$RawPPlot$mapping$colour),
    quote(.data$ColorCategory)
  )
  expect_identical(
    rlang::get_expr(result$FDRPlot$mapping$colour),
    quote(.data$ColorCategory)
  )

  raw_colors <- ggplot2::ggplot_build(result$RawPPlot)$data[[1]]$colour
  fdr_colors <- ggplot2::ggplot_build(result$FDRPlot)$data[[1]]$colour
  expect_identical(raw_colors, fdr_colors)
  expect_true(all(raw_colors[expected_categories == "Unmapped"] == "grey70"))
})

test_that("PlotVolcanoEffects accepts equivalent named-vector categories", {
  df_test <- MakeVolcanoTestData()
  vars_test <- VolcanoTestVariables()
  df_categories <- data.frame(
    Variable = c("strong", "moderate", "marginal"),
    Category = c("Immune", "Immune", "Metabolic")
  )
  vector_categories <- stats::setNames(
    df_categories$Category,
    df_categories$Variable
  )

  from_df <- PlotVolcanoEffects(
    df_test, vars_test, "outcome",
    InteractiveLabels = FALSE,
    ColorBy = df_categories
  )
  from_vector <- PlotVolcanoEffects(
    df_test, vars_test, "outcome",
    InteractiveLabels = FALSE,
    ColorBy = vector_categories
  )

  expect_identical(
    from_vector$ResultsTable$ColorCategory,
    from_df$ResultsTable$ColorCategory
  )
  expect_identical(
    ggplot2::ggplot_build(from_vector$RawPPlot)$data[[1]]$colour,
    ggplot2::ggplot_build(from_df$RawPPlot)$data[[1]]$colour
  )
})

test_that("PlotVolcanoEffects validates category mappings", {
  df_test <- MakeVolcanoTestData()
  vars_test <- VolcanoTestVariables()

  expect_error(
    PlotVolcanoEffects(df_test, vars_test, "outcome", ColorBy = c("A", "B")),
    "named atomic vector"
  )
  expect_error(
    PlotVolcanoEffects(
      df_test,
      vars_test,
      "outcome",
      ColorBy = data.frame(
        Variable = c("strong", "strong"),
        Category = c("A", "B")
      )
    ),
    "conflicting category mappings for: strong"
  )
})

test_that("raw and FDR label modes use different significance families", {
  df_test <- MakeVolcanoTestData()
  vars_test <- VolcanoTestVariables()

  raw <- PlotVolcanoEffects(
    df_test, vars_test, "outcome",
    LabelMode = "raw",
    InteractiveLabels = FALSE
  )
  significant_alias <- PlotVolcanoEffects(
    df_test, vars_test, "outcome",
    LabelMode = "significant",
    InteractiveLabels = FALSE
  )
  fdr <- PlotVolcanoEffects(
    df_test, vars_test, "outcome",
    LabelMode = "fdr",
    InteractiveLabels = FALSE
  )

  vars_raw <- raw$ResultsTable$Variable[!is.na(raw$ResultsTable$PlotLabel)]
  vars_alias <- significant_alias$ResultsTable$Variable[
    !is.na(significant_alias$ResultsTable$PlotLabel)
  ]
  vars_fdr <- fdr$ResultsTable$Variable[!is.na(fdr$ResultsTable$PlotLabel)]

  expect_identical(vars_raw, c("strong", "marginal"))
  expect_identical(vars_alias, vars_raw)
  expect_identical(vars_fdr, "strong")
  expect_true(all(raw$ResultsTable$PValue[raw$ResultsTable$Variable %in% vars_raw] < 0.05))
  expect_true(all(fdr$ResultsTable$FDR[fdr$ResultsTable$Variable %in% vars_fdr] < 0.05))
})
