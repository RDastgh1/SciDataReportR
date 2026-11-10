GroupedRegressionFixture <- function() {
  set.seed(91)
  df_Test <- data.frame(x = rnorm(160), z = rnorm(160), g = factor(rep(c("Men", "Women"), each = 80)))
  df_Test$y <- df_Test$x * rep(c(0.5, 1.5), each = 80) + df_Test$z + rnorm(160)
  df_Test$y2 <- -df_Test$x + rnorm(160)
  df_Test
}

test_that("stratified models share pooled complete-case SDs and retain pooled adjustment", {
  df_Test <- GroupedRegressionFixture()
  df_Test$z[1:6] <- NA_real_
  result <- MakeUnivariateRegressionTable(df_Test, "y", "x", "z", group_var = "g", Standardize = TRUE)
  df_Pooled <- na.omit(df_Test)
  ratio <- sd(df_Pooled$x) / sd(df_Pooled$y)
  expect_setequal(result$Results$Group, c("Men", "Women", "Pooled"))
  for (group in c("Men", "Women")) {
    fit <- lm(y ~ x + z, df_Pooled[df_Pooled$g == group, ])
    expect_equal(result$Results$Estimate[result$Results$Group == group], coef(fit)[["x"]] * ratio)
  }
  expect_equal(result$Results$Estimate[result$Results$Group == "Pooled"], coef(lm(y ~ x + z + g, df_Pooled))[["x"]] * ratio)
  expect_match(result$Metadata$ScalingScope, "Pooled complete")
  expect_equal(result$Metadata$PooledScaling$SD[result$Metadata$PooledScaling$Variable == "x"], sd(df_Pooled$x))
  expect_setequal(result$Metadata$GroupN$Group, c("Men", "Women", "Pooled"))
  excluded <- MakeUnivariateRegressionTable(df_Test, "y", "x", group_var = "g", include_pooled = FALSE)
  expect_setequal(excluded$Results$Group, c("Men", "Women"))
})

test_that("ungrouped columns stay stable and FDR families stay within groups", {
  df_Test <- GroupedRegressionFixture()
  plain <- MakeUnivariateRegressionTable(df_Test, "y", "x")
  explicit <- MakeUnivariateRegressionTable(df_Test, "y", "x", group_var = NULL)
  expect_identical(names(plain$Results), names(explicit$Results))
  expect_equal(plain$Results, explicit$Results)
  expect_true(all(is.na(plain$Results$Group)))
  expect_true(all(is.na(plain$Results$PValueFDR)))
  for (scope in c("matrix", "per_predictor", "per_outcome")) {
    result <- MakeUnivariateRegressionTable(df_Test, c("y", "y2"), c("x", "z"), group_var = "g", fdr_scope = scope)
    df_Result <- result$Results
    family <- paste(df_Result$Group, if (scope == "per_predictor") df_Result$Predictor else if (scope == "per_outcome") df_Result$Outcome else "matrix")
    expect_equal(df_Result$PValueFDR, as.numeric(ave(df_Result$PValue, family, FUN = function(p) p.adjust(p, "BH"))))
  }
})

test_that("forest groups dodge, named lists bind, and plain plots remain available", {
  df_Test <- GroupedRegressionFixture()
  result <- MakeUnivariateRegressionTable(df_Test, "y", "x", group_var = "g", fdr_scope = "matrix")
  for (flip in c(FALSE, TRUE)) {
    plot <- PlotForestFromTable(result, Flip = flip, PVar = "PValueFDR", xlim = c(-2, 2), facet_scales = "free_y")
    built <- ggplot2::ggplot_build(plot)
    expect_equal(length(unique(built$data[[1]]$y)), 3)
    expect_equal(length(unique(built$data[[1]]$colour)), 3)
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
  tables <- split(result$Results[result$Results$Group != "Pooled", ], result$Results$Group[result$Results$Group != "Pooled"])
  plot <- PlotForestFromTable(tables)
  expect_equal(length(unique(plot$data$Group)), 2)
  plain <- MakeUnivariateRegressionTable(df_Test, "y", "x")
  expect_s3_class(PlotForestFromTable(plain$Results), "ggplot")
  interactions <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x")
  annotated <- PlotForestFromTable(interactions$Results, show_interaction_p = TRUE)
  expect_equal(length(unique(annotated$data$Group)), 2)
  expect_s3_class(ggplot2::ggplotGrob(annotated), "gtable")
})
