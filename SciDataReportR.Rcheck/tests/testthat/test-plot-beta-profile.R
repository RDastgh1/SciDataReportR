MakeBetaProfileTestData <- function() {
  set.seed(1357)
  n <- 90
  age <- stats::rnorm(n, 50, 8)
  sex <- factor(rep(c("Female", "Male"), length.out = n))
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- 0.65 * x1 - 0.3 * x2 + 0.02 * age + as.numeric(sex) * 0.2 + stats::rnorm(n)
  x2[1:12] <- NA_real_

  data.frame(y = y, x1 = x1, x2 = x2, age = age, sex = sex)
}

test_that("PlotBetaProfile returns direct standardized model statistics", {
  df_test <- MakeBetaProfileTestData()
  result <- PlotBetaProfile(
    df_test,
    c("x1", "x2"),
    "y",
    covariates = c("age", "sex"),
    InteractiveLabels = FALSE
  )

  expect_named(result, c("Plot", "ResultsTable", "Metadata"))
  expect_named(
    result$ResultsTable,
    c(
      "Variable", "Label", "Category", "Beta", "SE", "CILow", "CIHigh",
      "PValue", "FDR", "N", "R", "AdjustedR", "Tooltip"
    )
  )
  df_model <- df_test[, c("y", "x1", "age", "sex")]
  df_model <- df_model[stats::complete.cases(df_model), ]
  df_model$y_scaled <- as.numeric(scale(df_model$y))
  df_model$x_scaled <- as.numeric(scale(df_model$x1))
  fit <- stats::lm(y_scaled ~ x_scaled + age + sex, data = df_model)
  fit_table <- summary(fit)$coefficients
  fit_ci <- stats::confint(fit, "x_scaled", level = 0.95)
  x1_result <- result$ResultsTable[result$ResultsTable$Variable == "x1", ]

  expect_equal(x1_result$Beta, unname(fit_table["x_scaled", "Estimate"]))
  expect_equal(x1_result$SE, unname(fit_table["x_scaled", "Std. Error"]))
  expect_equal(x1_result$PValue, unname(fit_table["x_scaled", "Pr(>|t|)"]))
  expect_equal(x1_result$CILow, unname(fit_ci[1]))
  expect_equal(x1_result$CIHigh, unname(fit_ci[2]))
  expect_equal(
    result$ResultsTable$FDR,
    stats::p.adjust(result$ResultsTable$PValue, method = "fdr")
  )
  expect_equal(result$ResultsTable$N, c(90L, 78L))
  expect_s3_class(result$Plot, "ggplot")
  expect_identical(result$Metadata$ConfidenceLevel, 0.95)
  expect_identical(result$Metadata$NModels, 2L)
})

test_that("unadjusted beta and R are the zero-order correlation", {
  df_test <- MakeBetaProfileTestData()
  result <- PlotBetaProfile(df_test, "x1", "y", InteractiveLabels = FALSE)

  expected_r <- stats::cor(df_test$x1, df_test$y)
  expect_equal(result$ResultsTable$Beta, expected_r)
  expect_equal(result$ResultsTable$R, expected_r)
  expect_true(is.na(result$ResultsTable$AdjustedR))
})

test_that("PlotBetaProfile resolves labels in the documented priority", {
  df_test <- MakeBetaProfileTestData()
  attr(df_test$x1, "label") <- "Attribute label"
  codebook <- data.frame(Variable = "x1", Label = "Codebook label")

  result <- PlotBetaProfile(
    df_test,
    c("x1", "x2"),
    "y",
    codebook = codebook,
    InteractiveLabels = FALSE
  )
  raw <- PlotBetaProfile(
    df_test,
    c("x1", "x2"),
    "y",
    codebook = codebook,
    Relabel = FALSE,
    InteractiveLabels = FALSE
  )

  expect_identical(result$ResultsTable$Label, c("Codebook label", "x2"))
  expect_identical(raw$ResultsTable$Label, c("x1", "x2"))
})

test_that("category inputs and within-category sorting preserve supplied order", {
  df_test <- MakeBetaProfileTestData()
  df_test$x3 <- stats::rnorm(nrow(df_test))
  predictors <- c("x1", "x2", "x3")
  categories <- c("Second", "First", "Second")
  df_categories <- data.frame(Variable = predictors, Category = categories)
  named_categories <- stats::setNames(categories, predictors)

  positional <- PlotBetaProfile(
    df_test, predictors, "y",
    VariableCategories = categories,
    Sort = "within_category_effect",
    InteractiveLabels = FALSE
  )
  mapped <- PlotBetaProfile(
    df_test, predictors, "y",
    VariableCategories = df_categories,
    Sort = "within_category_effect",
    InteractiveLabels = FALSE
  )
  named <- PlotBetaProfile(
    df_test, predictors, "y",
    VariableCategories = named_categories,
    Sort = "within_category_effect",
    InteractiveLabels = FALSE
  )

  expect_identical(positional$ResultsTable$Category, c("Second", "Second", "First"))
  expect_identical(mapped$ResultsTable$Variable, positional$ResultsTable$Variable)
  expect_identical(named$ResultsTable$Variable, positional$ResultsTable$Variable)
  expect_identical(mapped$ResultsTable$Category, positional$ResultsTable$Category)
})

test_that("all PlotBetaProfile sorting modes implement their stated order", {
  df_test <- MakeBetaProfileTestData()
  df_test$x3 <- stats::rnorm(nrow(df_test))
  predictors <- c("x3", "x1", "x2")

  original <- PlotBetaProfile(df_test, predictors, "y", InteractiveLabels = FALSE)
  pvalue <- PlotBetaProfile(df_test, predictors, "y", Sort = "pvalue", InteractiveLabels = FALSE)
  fdr <- PlotBetaProfile(df_test, predictors, "y", Sort = "fdr", InteractiveLabels = FALSE)
  effect <- PlotBetaProfile(df_test, predictors, "y", Sort = "effect", InteractiveLabels = FALSE)

  expect_identical(original$ResultsTable$Variable, predictors)
  expect_true(all(diff(pvalue$ResultsTable$PValue) >= 0))
  expect_true(all(diff(fdr$ResultsTable$FDR) >= 0))
  expect_true(all(diff(abs(effect$ResultsTable$Beta)) <= 0))
})

test_that("profile plot follows PlotZScore visual conventions", {
  df_test <- MakeBetaProfileTestData()
  categorized <- PlotBetaProfile(
    df_test,
    c("x1", "x2"),
    "y",
    VariableCategories = c("A", "B"),
    RemoveXAxisLabels = FALSE
  )
  single_color <- PlotBetaProfile(df_test, c("x1", "x2"), "y")
  built <- ggplot2::ggplot_build(categorized$Plot)

  expect_identical(categorized$Plot$labels$y, "Standardized Beta")
  expect_equal(categorized$Plot$theme$axis.text.x$angle, 45)
  expect_true(any(vapply(built$plot$layers, function(x) inherits(x$geom, "GeomHline"), logical(1))))
  expect_true("text" %in% names(categorized$Plot$layers[[2]]$mapping))
  expect_identical(single_color$Plot$theme$legend.position, "none")
  expect_s3_class(categorized$Plot$scales$get_scales("colour"), "ScaleDiscrete")
})

test_that("invalid predictors are removed with informative warnings", {
  df_test <- MakeBetaProfileTestData()
  df_test$constant <- 1
  df_test$character <- "value"

  expect_warning(
    PlotBetaProfile(df_test, c("x1", "constant"), "y", InteractiveLabels = FALSE),
    "constant.*Zero variance"
  )
  result <- suppressWarnings(PlotBetaProfile(
      df_test,
      c("x1", "missing", "character", "constant", "y", "age", "x1"),
      "y",
      covariates = "age",
      InteractiveLabels = FALSE
    ))
  expect_identical(result$ResultsTable$Variable, "x1")
  expect_error(
    suppressWarnings(PlotBetaProfile(df_test, c("constant", "missing"), "y")),
    "No predictors"
  )
})

test_that("PlotBetaProfile validates public inputs", {
  df_test <- MakeBetaProfileTestData()

  expect_error(PlotBetaProfile(1, "x1", "y"), "data frame")
  expect_error(PlotBetaProfile(df_test, character(0), "y"), "nonempty")
  expect_error(PlotBetaProfile(df_test, "x1", "missing"), "not found")
  expect_error(PlotBetaProfile(transform(df_test, y = factor(y)), "x1", "y"), "numeric")
  expect_error(PlotBetaProfile(df_test, "x1", "y", covariates = "missing"), "covariates")
  expect_error(PlotBetaProfile(df_test, "x1", "y", AdjustMethod = "invalid"), "AdjustMethod")
  expect_error(
    PlotBetaProfile(df_test, "x1", "y", codebook = data.frame(Variable = "x1")),
    "Variable and Label"
  )
  expect_error(PlotBetaProfile(df_test, "x1", "y", Sort = "invalid"), "arg")
})

test_that("PlotVolcanoEffects continuous public structure remains unchanged", {
  df_test <- MakeBetaProfileTestData()
  result <- PlotVolcanoEffects(
    df_test,
    c("x1", "x2"),
    "y",
    covariates = c("age", "sex"),
    OutcomeType = "continuous",
    InteractiveLabels = FALSE
  )

  expect_named(result, c("RawPPlot", "FDRPlot", "ResultsTable"))
  expect_false(any(c("Beta", "SE", "CILow", "CIHigh") %in% names(result$ResultsTable)))
  profile <- PlotBetaProfile(
    df_test,
    c("x1", "x2"),
    "y",
    covariates = c("age", "sex"),
    InteractiveLabels = FALSE
  )
  expect_equal(result$ResultsTable$Effect, profile$ResultsTable$Beta)
  expect_equal(result$ResultsTable$PValue, profile$ResultsTable$PValue)
  expect_equal(result$ResultsTable$R, profile$ResultsTable$R)
  expect_equal(result$ResultsTable$AdjustedR, profile$ResultsTable$AdjustedR)
})
