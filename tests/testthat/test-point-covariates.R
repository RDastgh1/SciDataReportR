test_that("partial point-biserial results share the continuous correlation engine", {
  set.seed(281)
  n <- 120
  df_Test <- data.frame(check.names = FALSE,
    `age years` = rnorm(n), `site group` = factor(rep(c("A", "B", "C"), 40)))
  df_Test$binary <- ifelse(df_Test$`age years` + rnorm(n) > 0, "yes", "no")
  df_Test$outcome <- 3 * df_Test$`age years` + rnorm(n)
  df_Test$other <- rnorm(n)
  df_Test$outcome[1:9] <- NA
  df_Test$other[10:22] <- NA
  df_Test$`age years`[23:26] <- NA
  df_Test$`site group`[27:30] <- NA
  attr(df_Test$outcome, "label") <- "Outcome label"
  vars_Covariates <- c("age years", "site group")
  for (positive in c("yes", "no")) {
    df_Map <- createBinaryMapping(df_Test, "binary")
    df_Map$PositiveLevel <- positive
    df_Coded <- df_Test
    df_Coded$binary <- as.integer(df_Test$binary == positive)
    expected <- PlotCorrelationsHeatmap(df_Coded, "binary", c("outcome", "other"),
      covariates = vars_Covariates)
    actual <- PlotPointCorrelationsHeatmap(df_Test, "binary", c("outcome", "other"),
      covariates = vars_Covariates, binary_map = df_Map)
    df_Result <- actual$Unadjusted$PvalTable
    expect_equal(df_Result$correlation, as.numeric(expected$Unadjusted$r))
    expect_equal(df_Result$p_value, as.numeric(expected$Unadjusted$p))
    expect_equal(df_Result$nPairs, as.numeric(expected$Unadjusted$npairs))
    expect_equal(actual$Unadjusted$plot$data$correlation, df_Result$correlation)
    expect_equal(unname(df_Result$YLabel[1]), "Outcome label")
    df_Complete <- na.omit(df_Coded[c("binary", "outcome", vars_Covariates)])
    design <- model.matrix(~ ., df_Complete[vars_Covariates])
    r <- cor(lm.fit(design, df_Complete$binary)$residuals,
      lm.fit(design, df_Complete$outcome)$residuals)
    residual_df <- nrow(df_Complete) - qr(design)$rank - 1
    expect_equal(df_Result$p_value[1], 2 * pt(-abs(r * sqrt(residual_df / (1 - r^2))), residual_df))
    raw <- cor(df_Complete$binary, df_Complete$outcome)
    expect_gt(abs(raw - r), 0.2)
  }
})

test_that("shared pair engine handles collapsed and singular pairs", {
  df_Test <- data.frame(x = rep(0:1, 10), y = seq_len(20), z = rep(1:5, 4))
  df_Test$duplicate <- df_Test$z
  singular <- SciDataReportR:::ScidrCorrelationPair(df_Test, "pearson", 3, 1e-12, "y", "x")
  regular <- SciDataReportR:::ScidrCorrelationPair(df_Test[c("x", "y", "z")], "pearson", 3, 1e-12, "y", "x")
  expect_equal(singular, regular)
  df_Test$x <- 1
  expect_true(is.na(SciDataReportR:::ScidrCorrelationPair(df_Test, "pearson", 3, 1e-12, "y", "x")$r))
  expect_true(is.na(SciDataReportR:::ScidrCorrelationPair(df_Test[1:2, ], "pearson", 3, 1e-12, "y", "x")$p))
})
