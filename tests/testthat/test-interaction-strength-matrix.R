test_that("interaction shading separates strength, ending sign, and reversal", {
  slopes <- list(c(0.2, 0.8), c(-0.2, -0.8), c(0.8, 0.2),
                 c(-0.8, -0.2), c(-0.2, 0.8), c(0.8, -0.2), c(-0.5, 0.5))
  for (slope_pair in slopes) {
    df_Test <- data.frame(
      x = rep(seq(-2, 2, length.out = 80), 2),
      group = factor(rep(c("Reference", "Comparison"), each = 80),
                     levels = c("Reference", "Comparison"))
    )
    design <- stats::model.matrix(~ x * group, df_Test)
    noise <- stats::lm.fit(design, sin(seq_len(160)))$residuals * 0.03
    df_Test$y <- df_Test$x * rep(slope_pair, each = 80) + noise
    result <- PlotInteractionEffectsMatrix(df_Test, interVar = "group",
      outcome_vars = "y", predictor_vars = "x", Relabel = FALSE)
    df_Result <- result$raw_data
    ratio <- stats::sd(df_Test$x) / stats::sd(df_Test$y)
    expect_equal(df_Result$Slope1Std, slope_pair[1] * ratio, ignore_attr = TRUE)
    expect_equal(df_Result$Slope2Std, slope_pair[2] * ratio, ignore_attr = TRUE)
    expect_equal(df_Result$StrengthChange,
                 (abs(slope_pair[2]) - abs(slope_pair[1])) * ratio,
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(df_Result$DirectionReversal, prod(slope_pair) < 0)
    if (abs(slope_pair[1]) == abs(slope_pair[2])) {
      expect_identical(df_Result$StrengthChange, 0)
      expect_equal(result$Unadjusted$plot$data$FillValue, df_Result$SlopeDeltaStd)
    }
    expect_lt(df_Result$P, 0.05)
    expect_equal(result$Unadjusted$plot$data$FillValue,
                 result$FDRCorrected$plot$data$FillValue)
    expect_s3_class(ggplot2::ggplotGrob(result$Unadjusted$plot), "gtable")
    expect_s3_class(ggplot2::ggplotGrob(result$FDRCorrected$plot), "gtable")
  }
})

test_that("standardized strength uses model complete cases and is unit invariant", {
  set.seed(42)
  df_Test <- data.frame(x = rnorm(100), moderator = rnorm(100), z = rnorm(100))
  df_Test$y <- with(df_Test, 0.2 * x + 1.5 * x * moderator + z + rnorm(100, sd = 0.1))
  df_Test$z[1:10] <- NA_real_
  Fit <- function(df_Input) PlotInteractionEffectsMatrix(df_Input,
    interVar = "moderator", outcome_vars = "y", predictor_vars = "x",
    covariates = "z", Relabel = FALSE)
  result <- Fit(df_Test)
  df_Model <- stats::model.frame(y ~ x * moderator + z, data = df_Test)
  expect_equal(result$raw_data$Slope1Std,
    result$raw_data$Slope1 * sd(df_Model$x) / sd(df_Model$y))
  df_Test$x <- df_Test$x * 100
  df_Test$y <- df_Test$y * 10
  expect_equal(Fit(df_Test)$raw_data$StrengthChange, result$raw_data$StrengthChange)
})
