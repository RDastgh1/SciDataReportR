InteractionMatrixFixture <- function(levels = c("A", "B"), n = 180L) {
  set.seed(205)
  df_Test <- data.frame(x = rnorm(n * length(levels)),
    z = rnorm(n * length(levels)),
    g = factor(rep(levels, each = n), levels = levels))
  slopes <- c(-0.4, 0.9, 1.6)[seq_along(levels)]
  df_Test$y <- 0.7 + df_Test$x * rep(slopes, each = n) +
    0.6 * df_Test$z + rnorm(nrow(df_Test), sd = 0.25)
  df_Test
}

test_that("binary slopes and intervals use pooled complete-case scaling", {
  df_Test <- InteractionMatrixFixture()
  df_Test$z[1:12] <- NA_real_
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "z")
  model <- lm(y ~ x * g + z, df_Test)
  df_Model <- model.frame(model)
  scale_factor <- sd(df_Model$x) / sd(df_Model$y)
  df_Result <- result$Results
  expect_equal(df_Result$Term, c("Interaction", "Slope: A", "Slope: B"))
  expect_equal(df_Result$Estimate, unname(c(coef(model)["x:gB"], coef(model)["x"],
    coef(model)["x"] + coef(model)["x:gB"])))
  expect_equal(df_Result$EstimateStd, df_Result$Estimate * scale_factor)
  expect_equal(df_Result$StdErrorStd, df_Result$StdError * scale_factor)
  expect_equal(df_Result$ConfLowStd, df_Result$ConfLow * scale_factor)
  expect_equal(df_Result$ConfHighStd, df_Result$ConfHigh * scale_factor)
  expect_true(all(df_Result$ConfLow < c(1.3, -0.4, 0.9)))
  expect_true(all(df_Result$ConfHigh > c(1.3, -0.4, 0.9)))
  expect_equal(df_Result$N, rep(nobs(model), 3))
  expect_equal(df_Result$NGroup[-1], as.numeric(table(df_Model$g)))
  for (name in c("Slope1Std", "Slope2Std", "SlopeDeltaStd")) {
    expect_equal(result$Unadjusted[[name]][1, 1], result$raw_data[[name]])
    expect_equal(result$FDRCorrected[[name]], result$Unadjusted[[name]])
  }
  weights <- c(0, 1, 0, 0, 1)
  names(weights) <- names(coef(model))
  weights[] <- 0
  weights[c("x", "x:gB")] <- 1
  se <- sqrt(drop(t(weights) %*% vcov(model) %*% weights))
  expect_equal(df_Result$StdError[3], se)
  expect_equal(df_Result$PValue[3], 2 * pt(abs(df_Result$Estimate[3] / se), df.residual(model), lower.tail = FALSE))
  expect_equal(result$Unadjusted$C[1, 1], coef(model)[["x:gB"]])
})

test_that("full covariate interactions equal stratified fits", {
  df_Test <- InteractionMatrixFixture(c("A", "B", "C"))
  df_Test$y <- df_Test$y + df_Test$z * as.numeric(df_Test$g) * 2
  df_Test$z[1:5] <- NA_real_
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "z",
                                        covariate_interactions = TRUE)
  df_Slopes <- result$Results[result$Results$Term != "Interaction", ]
  for (level in levels(df_Test$g)) {
    stratified <- lm(y ~ x + z, df_Test[df_Test$g == level, ])
    expect_equal(df_Slopes$Estimate[df_Slopes$Level == level], coef(stratified)[["x"]])
  }
  expect_true(result$Metadata$Arguments$covariate_interactions)
  expect_match(result$Metadata$InteractionTest, "Moderator-specific")
})

test_that("multilevel inference uses a joint test and records the largest pair", {
  df_Test <- InteractionMatrixFixture(c("A", "B", "C"))
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "z")
  full <- lm(y ~ x * g + z, df_Test)
  reduced <- lm(y ~ x + g + z, df_Test)
  expect_equal(result$raw_data$P, anova(reduced, full)$`Pr(>F)`[2])
  df_Interaction <- result$Results[result$Results$Term == "Interaction", ]
  expect_equal(df_Interaction$JointDF, 2)
  expect_equal(df_Interaction$ReferenceLevel, "A")
  expect_equal(df_Interaction$ComparisonLevel, "C")
  expect_equal(df_Interaction$LargestLevelDifference, result$raw_data$SlopeDeltaStd)
  expect_equal(nrow(result$Results), 4)
  expect_equal(result$raw_data$Slope2, coef(full)[["x"]] + coef(full)[["x:gC"]])
  expect_equal(df_Interaction$PValue, result$raw_data$P)
})

test_that("continuous slopes use model-sample moderator values and covariance", {
  set.seed(412)
  df_Test <- data.frame(x = rnorm(500), moderator = rnorm(500, 2, 0.7), z = rnorm(500))
  df_Test$y <- 0.5 * df_Test$x + 1.2 * df_Test$x * df_Test$moderator +
    0.3 * df_Test$z + rnorm(500, sd = 0.2)
  df_Test$z[1:20] <- NA_real_
  df_Test$moderator[1:20] <- 100
  result <- PlotInteractionEffectsMatrix(df_Test, "moderator", "y", "x", "z")
  model <- lm(y ~ x * moderator + z, df_Test)
  df_Model <- model.frame(model)
  values <- mean(df_Model$moderator) + c(-1, 1) * sd(df_Model$moderator)
  df_Slopes <- result$Results[result$Results$Term != "Interaction", ]
  expect_equal(df_Slopes$ModeratorValue, values)
  expect_equal(df_Slopes$Estimate, unname(coef(model)["x"] + coef(model)["x:moderator"] * values))
  expect_true(all(df_Slopes$ConfLow < 0.5 + 1.2 * values))
  expect_true(all(df_Slopes$ConfHigh > 0.5 + 1.2 * values))
  expect_equal(result$Results$Estimate[1], coef(model)[["x:moderator"]])
  expect_equal(result$raw_data$SlopeDeltaStd,
    diff(df_Slopes$EstimateStd), ignore_attr = TRUE)
  expect_true(all(is.na(df_Slopes$NGroup)))
})

test_that("redistributing group outcome variance retains pooled scaling", {
  df_Test <- InteractionMatrixFixture()
  df_Changed <- df_Test
  indices <- df_Changed$g == "B"
  df_Changed$y[indices] <- mean(df_Changed$y[indices]) +
    1.7 * (df_Changed$y[indices] - mean(df_Changed$y[indices]))
  # Restore the original full-sample SD, without equalizing the group SDs.
  df_Changed$y <- mean(df_Test$y) + (df_Changed$y - mean(df_Changed$y)) *
    sd(df_Test$y) / sd(df_Changed$y)
  expect_equal(sd(df_Changed$y), sd(df_Test$y))
  results <- lapply(list(df_Test, df_Changed), function(df_Input) {
    PlotInteractionEffectsMatrix(df_Input, "g", "y", "x")$raw_data
  })
  ratios <- vapply(results, function(df_Result) df_Result$Slope2Std / df_Result$Slope2, numeric(1))
  expect_equal(ratios[1], ratios[2])
  expect_equal(ratios[1], sd(df_Test$x) / sd(df_Test$y))
  # Scaling alone stays invariant: actual slopes may change with the outcomes.
  expect_false(isTRUE(all.equal(results[[1]]$Slope2Std, results[[2]]$Slope2Std)))
})

test_that("term assignment handles substrings, non-syntactic names and contrasts", {
  df_Test <- InteractionMatrixFixture()
  names(df_Test)[names(df_Test) == "z"] <- "x longer"
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "x longer",
                                        covariate_interactions = TRUE)
  model <- lm(y ~ x * g + `x longer` * g, df_Test)
  expect_equal(result$raw_data$C, coef(model)[["x:gB"]])
  contrasts(df_Test$g) <- contr.sum(2)
  sum_result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "x longer",
                                            covariate_interactions = TRUE)
  expect_equal(sum_result$raw_data$Slope1, result$raw_data$Slope1)
  expect_equal(sum_result$raw_data$Slope2, result$raw_data$Slope2)
  expect_equal(sum_result$raw_data$P, result$raw_data$P)
  expect_warning(PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", c("x longer", "g")),
                 "Dropping interVar")
})

test_that("failures warn through the shared contract and remain grey", {
  df_Test <- InteractionMatrixFixture()
  df_Test$x <- 1
  expect_warning(result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x"),
                 "outcome y, predictor x.*rank deficient")
  expect_equal(result$Results$Status, "failed")
  expect_equal(nrow(result$Metadata$Failed), 1)
  built <- ggplot2::ggplot_build(result$Unadjusted$plot)
  expect_equal(built$data[[2]]$fill, "#BDBDBD")
  expect_s3_class(ggplot2::ggplotGrob(result$FDRCorrected$plot), "gtable")
})

test_that("plot thresholds, no-hit notes, common limits and FDR families agree", {
  df_Test <- InteractionMatrixFixture()
  set.seed(999)
  df_Test$null <- rnorm(nrow(df_Test))
  df_Test$other <- rnorm(nrow(df_Test))
  for (scope in c("matrix", "per_outcome", "per_predictor")) {
    result <- PlotInteractionEffectsMatrix(df_Test, "g", c("y", "null"), c("x", "other"),
      fdr_scope = scope, Pthresh = 0.01, PVar = "PValueFDR")
    df_Raw <- result$raw_data
    expected <- if (scope == "matrix") p.adjust(df_Raw$P, "BH") else {
      groups <- if (scope == "per_outcome") df_Raw$Outcome else df_Raw$Predictor
      ave(df_Raw$P, groups, FUN = function(p) p.adjust(p, "BH"))
    }
    expect_equal(df_Raw$P_FDR, expected)
    df_Interaction <- result$Results[result$Results$Term == "Interaction", ]
    expect_equal(df_Interaction$PValueFDR, expected)
    for (adjusted in c(FALSE, TRUE)) {
      plot <- if (adjusted) result$FDRCorrected$plot else result$Unadjusted$plot
      p_values <- if (adjusted) df_Raw$P_FDR else df_Raw$P
      expect_equal(is.finite(plot$data$FillValue), p_values <= 0.01)
      expect_equal(plot$data$FillValue[p_values <= 0.01], df_Raw$SlopeDeltaStd[p_values <= 0.01])
      expect_equal(plot$scales$get_scales("fill")$limits, result$SlopeDifferenceLimits)
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    }
  }
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "null", "other", Pthresh = 1e-10)
  built <- ggplot2::ggplot_build(result$Unadjusted$plot)
  labels <- unlist(lapply(built$data, function(df_Layer) df_Layer$label))
  expect_true(any(grepl("No interactions with p <=", labels)))
  expect_true(all(is.na(result$Unadjusted$plot$data$FillValue)))
})

test_that("threshold boundary is inclusive and failures do not hide valid cells", {
  df_Test <- InteractionMatrixFixture()
  df_Test$constant <- 1
  baseline <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x")
  at_boundary <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x",
                                             Pthresh = baseline$raw_data$P)
  expect_true(is.finite(at_boundary$Unadjusted$plot$data$FillValue))
  expect_true(at_boundary$Results$Significant[at_boundary$Results$Term == "Interaction"])
  expect_warning(mixed <- PlotInteractionEffectsMatrix(df_Test, "g", "y", c("x", "constant")),
                 "outcome y, predictor constant")
  expect_equal(mixed$raw_data$Status, c("ok", "failed"))
  expect_equal(mixed$Metadata$ModelCounts$Fit, 1L)
  expect_equal(mixed$Metadata$ModelCounts$Failed, 1L)
  expect_equal(mixed$Unadjusted$Slope1[1, 1], baseline$Unadjusted$Slope1[1, 1])
  expect_true(is.na(mixed$Unadjusted$Slope1[1, 2]))
  expect_error(PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "typo"), "typo")
})

test_that("multilevel selected differences can exclude the baseline group", {
  df_Test <- InteractionMatrixFixture(c("A", "B", "C"))
  df_Test$y <- df_Test$y + df_Test$x * rep(c(0.4, -2.9, 0.4), each = 180)
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", "z",
                                        covariate_interactions = TRUE)
  df_Interaction <- result$Results[result$Results$Term == "Interaction", ]
  expect_equal(df_Interaction$ReferenceLevel, "B")
  expect_equal(df_Interaction$ComparisonLevel, "C")
  full <- lm(y ~ x * g + z * g, df_Test)
  reduced <- lm(y ~ x + g + z * g, df_Test)
  expect_equal(df_Interaction$PValue, anova(reduced, full)$`Pr(>F)`[2])
  expect_equal(df_Interaction$Estimate, coef(full)[["x:gC"]] - coef(full)[["x:gB"]])
  expect_true(result$raw_data$DirectionReversal)
})

test_that("all numeric and factor covariates can have group-specific effects", {
  df_Test <- InteractionMatrixFixture()
  df_Test$h <- factor(rep(c("Low", "High"), length.out = nrow(df_Test)))
  df_Test$y <- df_Test$y + 0.8 * (df_Test$h == "High") * (df_Test$g == "B")
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x", c("z", "h"),
                                        covariate_interactions = TRUE)
  df_Slopes <- result$Results[result$Results$Term != "Interaction", ]
  for (level in levels(df_Test$g)) {
    stratified <- lm(y ~ x + z + h, df_Test[df_Test$g == level, ])
    expect_equal(df_Slopes$Estimate[df_Slopes$Level == level], coef(stratified)[["x"]])
  }
})
