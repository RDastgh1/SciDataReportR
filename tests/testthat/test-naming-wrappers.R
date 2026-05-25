test_that("standardization compatibility wrappers return identical objects", {
  df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))

  expect_identical(
    CalcZScore(df, variables = c("a", "b"), RetainLabels = FALSE),
    CreateZScoreObject(df, variables = c("a", "b"), RetainLabels = FALSE)
  )

  expect_identical(
    CalcMScore(df, variables = c("a", "b"), RetainLabels = FALSE),
    CreateMScoreObject(df, variables = c("a", "b"), RetainLabels = FALSE)
  )
})

test_that("Project_ZScore remains a compatibility wrapper", {
  df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  z_obj <- CreateZScoreObject(df, variables = c("a", "b"), RetainLabels = FALSE)

  expect_identical(
    Project_ZScore(df, parameters = z_obj, ParameterInputType = "ZScoreObj"),
    ProjectZScore(df, parameters = z_obj, ParameterInputType = "ZScoreObj")
  )
})

test_that("data dictionary compatibility wrapper returns identical output", {
  skip_if_not_installed("codebook")

  df <- data.frame(
    group = factor(c("A", "B", "A")),
    value = c(1, 2, 3)
  )

  expect_identical(
    Make_DataDictionary(df),
    MakeDataDictionary(df)
  )
})

test_that("plot compatibility wrappers return ggplot objects", {
  skip_if_not_installed("plotrix")
  skip_if_not_installed("rstatix")

  df <- data.frame(
    group = rep(c("A", "B"), each = 4),
    x = c(1, 2, 3, 4, 2, 3, 4, 5),
    y = c(2, 3, 4, 5, 3, 4, 5, 6)
  )

  expect_s3_class(
    CreateZScorePlot(df, TargetVar = "group", Variables = c("x", "y")),
    "ggplot"
  )
  expect_s3_class(
    PlotZScore(df, TargetVar = "group", Variables = c("x", "y")),
    "ggplot"
  )
})

test_that("PCA and MCA compatibility wrappers expose the same structures", {
  skip_if_not_installed("psych")
  skip_if_not_installed("xtable")
  skip_if_not_installed("FactoMineR")

  pca_df <- data.frame(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(2, 3, 4, 5, 6, 7),
    c = c(6, 5, 4, 3, 2, 1)
  )

  pca_old <- CreatePCATable(pca_df, VarsToReduce = c("a", "b", "c"), numComponents = 2)
  pca_new <- CreatePCAObject(pca_df, VarsToReduce = c("a", "b", "c"), numComponents = 2)

  expect_named(pca_old, names(pca_new))

  mca_df <- data.frame(
    a = factor(c("low", "low", "mid", "mid", "high", "high")),
    b = factor(c("yes", "no", "yes", "no", "yes", "no")),
    c = factor(c("x", "x", "y", "y", "z", "z"))
  )

  mca_old <- CreateMCATable(mca_df, VarsToReduce = c("a", "b", "c"), numComponents = 2)
  mca_new <- CreateMCAObject(mca_df, VarsToReduce = c("a", "b", "c"), numComponents = 2)

  expect_named(mca_old, names(mca_new))
})
