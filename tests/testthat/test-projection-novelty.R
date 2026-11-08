test_that("calibrated combined LOO flags match the requested training rate", {
  model <- ProjectionTestModel()
  calibration <- model$ModelInfo_SOM$ProjectionReference$calibration
  df_Individual <- model$ProbFit$individual
  far <- df_Individual$SOM_Distance_LOO > calibration$overall_cutoff |
    df_Individual$SOM_Distance_LOO > calibration$cluster_cutoffs[as.character(df_Individual$Cluster)]
  expect_lte(abs(sum(far) - length(far) * .05), 1)
  expect_equal(mean(far), calibration$training_rate)
  expect_gte(calibration$q_star, .95)
})

test_that("the embedded projection cohort gives no qualifying group or excess", {
  model <- ProjectionTestModel()
  result <- suppressWarnings(suppressMessages(ProjectCluster(model, EmbeddedTestData("Projection"))))
  expect_false(any(result$NoveltyGroups$groups$Candidate))
  expect_false(result$NoveltyGroups$overall_excess$Excess)
  expect_equal(result$NoveltyGroups$overall_excess$ExpectedN,
    nrow(EmbeddedTestData("Projection")) * model$ModelInfo_SOM$ProjectionReference$calibration$training_rate)
})

test_that("compact embedded-data shifts are detected across seeds", {
  model <- ProjectionTestModel()
  vars_Features <- model$vars_used
  df_Training <- EmbeddedTestData(variables = c(vars_Features, "TruthCluster"))
  df_Projection <- EmbeddedTestData("Projection")
  mean_profile <- colMeans(df_Training[df_Training$TruthCluster == levels(df_Training$TruthCluster)[1], vars_Features])
  sd_profile <- vapply(df_Training[vars_Features], stats::sd, numeric(1))
  performance <- vapply(1:10, function(seed) {
    set.seed(seed)
    df_Mixed <- df_Projection[sample(nrow(df_Projection), 130), vars_Features]
    for (variable in vars_Features) {
      df_Mixed[101:130, variable] <- mean_profile[[variable]] + 1.5 * sd_profile[[variable]] +
        rnorm(30, sd = .1 * sd_profile[[variable]])
    }
    result <- suppressWarnings(suppressMessages(ProjectCluster(model, df_Mixed)))
    groups <- result$NoveltyGroups$groups
    members <- unique(unlist(groups$Members[groups$Candidate]))
    expect_true(all(vapply(seq_len(nrow(groups)), function(i) grepl(paste0("n=", groups$N[i]), groups$Label[i], fixed = TRUE), logical(1))))
    c(found = any(groups$Candidate), sensitivity = mean(101:130 %in% members), false_inclusion = mean(1:100 %in% members))
  }, numeric(3))
  expect_gt(mean(performance["found", ]), .5)
  expect_gt(mean(performance["sensitivity", ]), .5)
})
