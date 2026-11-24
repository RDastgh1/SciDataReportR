HeatmapCovariateData <- function() {
  set.seed(718)
  df_Test <- data.frame(check.names = FALSE,
    `age years` = rnorm(120), site = factor(rep(c("A", "B", "C"), 40)))
  df_Test$g <- factor(ifelse(df_Test$`age years` + rnorm(120) > 0, "yes", "no"))
  df_Test$h <- factor(ifelse(df_Test$`age years` + rnorm(120) > 0, "high", "low"))
  df_Test$a <- 2 * df_Test$`age years` + rnorm(120)
  df_Test$b <- df_Test$`age years` + rnorm(120)
  df_Test$a[1:5] <- NA
  df_Test$h[6:9] <- NA
  df_Test$`age years`[10:13] <- NA
  df_Test$site[14:16] <- NA
  attr(df_Test$g, "label") <- "Group label"
  df_Test
}

test_that("adjusted binary correlations match independent residual inference", {
  df_Test <- HeatmapCovariateData()
  vars_Covariates <- c("age years", "site")
  df_Map <- createBinaryMapping(df_Test, c("g", "h"))
  df_Coded <- df_Test
  for (variable in c("g", "h")) {
    df_Coded[[variable]] <- as.integer(df_Test[[variable]] == df_Map$PositiveLevel[df_Map$Variable == variable])
  }
  df_Complete <- na.omit(df_Coded[c("g", "h", vars_Covariates)])
  design <- model.matrix(~ ., df_Complete[vars_Covariates])
  r <- cor(lm.fit(design, df_Complete$g)$residuals,
    lm.fit(design, df_Complete$h)$residuals)
  residual_df <- nrow(df_Complete) - qr(design)$rank - 1
  expected_p <- 2 * pt(-abs(r * sqrt(residual_df / (1 - r^2))), residual_df)
  result <- PlotPhiHeatmap(df_Test, c("g", "h"), covariates = vars_Covariates,
    binary_map = df_Map, triangle = "upper", cluster_rows = TRUE)
  df_Result <- result$Unadjusted$PvalTable
  df_Result <- df_Result[df_Result$XVar != df_Result$YVar, ]
  expect_equal(df_Result$Phi, rep(r, 2))
  expect_equal(df_Result$p_value, rep(expected_p, 2))
  expect_equal(df_Result$nPairs, rep(nrow(df_Complete), 2))
  expect_true(all(result$Results$EffectType == "Partial correlation"))
  expect_identical(result$Metadata$Arguments$covariates, vars_Covariates)
  expect_match(result$Metadata$FormulaTemplate, "rank")
  expect_true("Group label" %in% df_Result$XLabel)
  df_Map$PositiveLevel[df_Map$Variable == "g"] <- setdiff(levels(df_Test$g), df_Map$PositiveLevel[df_Map$Variable == "g"])
  reversed <- PlotPhiHeatmap(df_Test, c("g", "h"), covariates = vars_Covariates, binary_map = df_Map)
  expect_equal(reversed$Results$Estimate, -result$Results$Estimate)
  expect_equal(reversed$Results$PValue, result$Results$PValue)
  df_Test$duplicate <- df_Test$`age years`
  singular <- PlotPhiHeatmap(df_Test, c("g", "h"), covariates = c(vars_Covariates, "duplicate"))
  expect_equal(singular$Results$Estimate, result$Results$Estimate)
  expect_equal(singular$Results$PValue, result$Results$PValue)
})

test_that("directional adjusted blocks agree with standalone components", {
  df_Test <- HeatmapCovariateData()
  vars_Covariates <- c("age years", "site")
  for (scope in c("matrix", "per_outcome", "per_predictor")) {
    phi <- PlotPhiHeatmap(df_Test, c("g", "h"), covariates = vars_Covariates, fdr_scope = scope)
    continuous <- PlotCorrelationsHeatmap(df_Test, c("a", "b"), covariates = vars_Covariates, fdr_scope = scope)
    point <- PlotPointCorrelationsHeatmap(df_Test, c("g", "h"), c("a", "b"), covariates = vars_Covariates, fdr_scope = scope)
    square <- PlotDirectionalHeatmaps(df_Test, c("g", "h", "a", "b"), covariates = vars_Covariates, fdr_scope = scope)
    for (component in list(phi, continuous, point)) {
      df_Expected <- component$Results
      for (i in seq_len(nrow(df_Expected))) {
        df_Actual <- square$Results
        df_Actual <- df_Actual[df_Actual$Predictor == df_Expected$Predictor[i] & df_Actual$Outcome == df_Expected$Outcome[i], ]
        expect_equal(nrow(df_Actual), 1)
        for (column in c("Estimate", "PValue", "PValueFDR", "N")) {
          expect_equal(df_Actual[[column]], df_Expected[[column]][i], info = column)
        }
      }
    }
    rectangle <- PlotDirectionalHeatmaps(df_Test, predictor_vars = c("g", "a"),
      outcome_vars = c("h", "b"), covariates = vars_Covariates, fdr_scope = scope)
    df_Rectangle <- rectangle$Results
    df_Square <- square$Results
    for (i in seq_len(nrow(df_Rectangle))) {
      df_Expected <- df_Square[df_Square$Predictor == df_Rectangle$Predictor[i] & df_Square$Outcome == df_Rectangle$Outcome[i], ]
      expect_equal(df_Rectangle$Estimate[i], df_Expected$Estimate)
      expect_equal(df_Rectangle$PValue[i], df_Expected$PValue)
      expect_equal(df_Rectangle$N[i], df_Expected$N)
    }
    expect_true(all(rectangle$Results$EffectType == "Partial correlation"))
    expect_identical(square$Metadata$Arguments$covariates, vars_Covariates)
  }
})

test_that("covariate validation and defaults preserve legacy behavior", {
  df_Test <- HeatmapCovariateData()
  for (Compute in list(PlotPhiHeatmap, PlotDirectionalHeatmaps)) {
    expect_error(Compute(df_Test, c("g", "h"), covariates = "absent"), "absent.*covariates")
    expect_warning(result <- Compute(df_Test, c("g", "h", "site"), covariates = "site"), "Excluding covariates")
    expect_false(any(result$Results$Predictor == "site" | result$Results$Outcome == "site"))
    expect_false("site" %in% c(result$Metadata$Arguments$predictor_vars,
      result$Metadata$Arguments$outcome_vars))
    raw <- Compute(df_Test, c("g", "h"))
    explicit <- Compute(df_Test, c("g", "h"), covariates = NULL)
    expect_equal(raw$Unadjusted, explicit$Unadjusted)
    expect_equal(raw$Results, explicit$Results)
  }
  automatic <- PlotDirectionalHeatmaps(df_Test, covariates = c("age years", "site"))
  expect_false(any(c("age years", "site") %in% automatic$AxisOrder$columns))
  raw <- PlotPhiHeatmap(df_Test, c("g", "h"))
  df_Complete <- na.omit(df_Test[c("g", "h")])
  expected <- suppressWarnings(chisq.test(table(df_Complete$g, df_Complete$h), correct = FALSE)$p.value)
  expect_equal(raw$Results$PValue, rep(expected, nrow(raw$Results)))
})

test_that("failed adjustments propagate through safe names and mirrored cells", {
  df_Test <- HeatmapCovariateData()
  names(df_Test)[names(df_Test) == "g"] <- "binary group"
  df_Test$collapsed <- factor(rep("only", nrow(df_Test)))
  for (Compute in list(PlotPhiHeatmap, PlotDirectionalHeatmaps)) {
    expect_warning(result <- Compute(df_Test, c("binary group", "h"), covariates = "collapsed"), "adjust")
    expect_true(all(result$Results$Status == "failed"))
    expect_true(all(is.na(result$Results$PValue)))
    expect_true(all(result$Results$AdjustmentStatus == "unadjusted_fallback"))
    expect_true(all(result$Metadata$Failed$Stage == "adjustment"))
    expect_false(any(grepl("binary.group", result$Metadata$Failed$Predictor, fixed = TRUE)))
  }
  for (rectangular in c(FALSE, TRUE)) {
    warnings <- character()
    result <- withCallingHandlers(if (rectangular) {
      PlotDirectionalHeatmaps(df_Test, predictor_vars = c("binary group", "a"),
        outcome_vars = c("h", "b"), covariates = "collapsed")
    } else {
      PlotDirectionalHeatmaps(df_Test, c("binary group", "h", "a", "b"), covariates = "collapsed")
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_true(length(warnings) > 0)
    expect_true(all(grepl("adjust", warnings)))
    expect_true(all(result$Results$Status == "failed"))
    expect_true(all(is.na(result$Results$PValue)))
    expect_true(all(is.na(result$Unadjusted$data$P_adj)))
    expect_false(any(grepl("*", result$Unadjusted$data$stars, fixed = TRUE)))
  }
  df_Tiny <- df_Test[1:4, ]
  df_Tiny$h <- factor(c("low", "high", "low", "high"))
  df_Tiny$`binary group` <- factor(c("yes", "no", "no", "yes"))
  df_Tiny$z <- seq_len(4)
  df_Tiny$w <- c(1, 0, 0, 0)
  expect_warning(result <- PlotPhiHeatmap(df_Tiny, c("binary group", "h"), covariates = c("z", "w")), "Insufficient")
  expect_true(all(result$Results$Status == "failed"))
})

test_that("pair-specific adjustment failures leave healthy correction families intact", {
  df_Test <- HeatmapCovariateData()
  df_Test$h[df_Test$site != "A" | is.na(df_Test$site)] <- NA
  for (scope in c("matrix", "per_outcome", "per_predictor")) {
    result <- suppressWarnings(PlotDirectionalHeatmaps(df_Test, c("g", "h", "a", "b"),
      covariates = "site", fdr_scope = scope))
    df_Failed <- result$Results[result$Results$Predictor == "h" | result$Results$Outcome == "h", ]
    expect_true(all(df_Failed$Status == "failed"))
    expect_true(all(is.na(df_Failed$PValue)))
    expected <- PlotCorrelationsHeatmap(df_Test, c("a", "b"), covariates = "site", fdr_scope = scope)
    df_Actual <- result$Results[result$Results$Predictor %in% c("a", "b") &
      result$Results$Outcome %in% c("a", "b"), ]
    indices <- match(paste(expected$Results$Predictor, expected$Results$Outcome),
      paste(df_Actual$Predictor, df_Actual$Outcome))
    expect_equal(df_Actual$PValueFDR[indices], expected$Results$PValueFDR)
    expect_true(all(df_Actual$Status == "ok"))
    df_Point <- result$Unadjusted$data
    df_Point <- df_Point[df_Point$test == "Point Correlation" & df_Point$XVar == "g", ]
    expect_equal(df_Point$P_adj, ApplyFDRCorrection(df_Point$p_value,
      fdr_scope = scope, outcome_ids = df_Point$YVar, predictor_ids = df_Point$XVar))
  }
})


test_that("adjusted safe names and residual degeneracy remain inspectable", {
  df_Test <- HeatmapCovariateData()
  names(df_Test)[names(df_Test) == "g"] <- "binary group"
  result <- PlotDirectionalHeatmaps(df_Test, predictor_vars = "binary group",
    outcome_vars = c("h", "a"), covariates = c("age years", "site"))
  expect_true(all(result$Results$Status == "ok"))
  expect_true(all(result$Results$Predictor == "binary group"))
  expect_true(all(result$Results$PredictorLabel == "Group label"))
  expect_true("binary group" %in% result$BinaryMapping$Variable)
  df_Test$copy <- df_Test$`binary group`
  expect_warning(result <- PlotPhiHeatmap(df_Test, c("binary group", "h"),
    covariates = "copy"), "could not be estimated")
  expect_true(all(result$Results$Status == "failed"))
  expect_true(all(is.na(result$Results$Estimate)))
})
