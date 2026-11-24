test_that("conformal novelty p-values match their definition", {
  reference <- c(1, 2, 3, 4)
  expect_equal(.SOMConformalP(c(0, 2, 2.5, 5, NA), reference), c(5, 4, 3, 1, NA) / 5)
  expect_equal(.SOMTrainingConformalP(c(1, 2, 3, 4)), c(4, 3, 2, 1) / 4)
  columns <- .SOMNoveltyColumns(c(.001, .2, .9, NA), novelty_fdr = .10, group_screen_alpha = .05)
  expect_equal(columns$NoveltyQ[1:3], stats::p.adjust(c(.001, .2, .9), "BH"))
  expect_identical(columns$Flag_Outlier, c(TRUE, FALSE, FALSE, NA))
  expect_identical(columns$Flag_Unusual, c(TRUE, FALSE, FALSE, NA))
})

test_that("fit classes give group membership and outliers precedence over membership", {
  fit <- .SOMFitClass(c(TRUE, TRUE, TRUE, TRUE, FALSE), c(FALSE, FALSE, TRUE, TRUE, NA),
    c(FALSE, FALSE, FALSE, TRUE, NA), clear = c(TRUE, FALSE, FALSE, FALSE, NA))
  expect_identical(as.character(fit), c("Good fit", "Uncertain membership", "Outlier",
    "Member of candidate new group", NA))
  expect_identical(levels(.SOMFitClass(TRUE, FALSE)), .SOMHardFitClassLevels)
})

test_that("novelty permutations leave the caller's random stream unchanged", {
  group_z <- matrix(c(0, 1, 0, 1, 0, 1), 3)
  pool_z <- matrix(seq_len(20) / 7, 10)
  set.seed(42)
  expected <- stats::runif(3)
  set.seed(42)
  coherence <- .SOMGroupCoherence(group_z, pool_z, permutations = 19)
  expect_identical(stats::runif(3), expected)
  expect_true(coherence$P > 0 && coherence$P <= 1)
})

test_that("training people carry novelty columns on the projection scale", {
  model <- ProjectionTestModel()
  df_Individual <- model$ProbFit$individual
  expect_true(all(c("NoveltyP", "NoveltyQ", "Flag_Outlier", "MembershipProbability",
    "FitClass") %in% names(df_Individual)))
  expect_equal(df_Individual$NoveltyP, .SOMTrainingConformalP(df_Individual$SOM_Distance_LOO))
  expect_identical(levels(df_Individual$FitClass), .SOMFitClassLevels)
  expect_null(model$ModelInfo_SOM$ProjectionReference$calibration)
})

test_that("the embedded projection cohort gives no outliers, group or excess", {
  model <- ProjectionTestModel()
  result <- suppressWarnings(suppressMessages(ProjectCluster(model, EmbeddedTestData("Projection"))))
  df_Individual <- result$ProjectionFit$individual
  expect_equal(sum(df_Individual$Flag_Outlier), 0)
  expect_false(any(result$NoveltyGroups$groups$Candidate))
  expect_false(result$NoveltyGroups$overall_excess$Excess)
  expect_identical(levels(df_Individual$FitClass), .SOMFitClassLevels)
  expect_identical(result$DataWithClusters$FitClass, df_Individual$FitClass)
  expect_false(any(c("Projection_Fit_Class", "FitClass_Legacy", "Flag_FarFromMap") %in%
    names(result$DataWithClusters)))
  shift <- result$ProjectionFit$cohort_shift
  expect_identical(shift$Test, c("Distance shift", "Distance distribution",
    "Unusual-profile excess", "Node occupancy", "Cluster prevalence"))
  expect_gt(shift$P[shift$Test == "Node occupancy"], .05)
  expect_gt(shift$P[shift$Test == "Cluster prevalence"], .05)
  expect_s3_class(result$ProjectionFit$plots$novelty_p_hist, "ggplot")
})

test_that("a planted extra cluster is found as a group and as outliers across seeds", {
  model <- ProjectionTestModel()
  vars_Features <- model$vars_used
  df_Training <- EmbeddedTestData(variables = c(vars_Features, "TruthCluster"))
  df_Projection <- EmbeddedTestData("Projection")
  mean_profile <- colMeans(df_Training[df_Training$TruthCluster == levels(df_Training$TruthCluster)[1], vars_Features])
  sd_profile <- vapply(df_Training[vars_Features], stats::sd, numeric(1))
  performance <- vapply(1:5, function(seed) {
    set.seed(seed)
    df_Mixed <- df_Projection[sample(nrow(df_Projection), 130), vars_Features]
    for (variable in vars_Features) {
      df_Mixed[101:130, variable] <- mean_profile[[variable]] + 1.5 * sd_profile[[variable]] +
        rnorm(30, sd = .1 * sd_profile[[variable]])
    }
    result <- suppressWarnings(suppressMessages(ProjectCluster(model, df_Mixed)))
    groups <- result$NoveltyGroups$groups
    members <- unique(unlist(groups$Members[groups$Candidate]))
    df_Individual <- result$ProjectionFit$individual
    expect_true(all(vapply(seq_len(nrow(groups)), function(i)
      grepl(paste0("n=", groups$N[i]), groups$Label[i], fixed = TRUE), logical(1))))
    expect_identical(df_Individual$Flag_NewGroupMember, df_Individual$.row_id %in% members)
    c(found = any(groups$Candidate), sensitivity = mean(101:130 %in% members),
      false_inclusion = mean(1:100 %in% members),
      planted_outlier = mean(df_Individual$Flag_Outlier[101:130]),
      shift_p = result$ProjectionFit$cohort_shift$P[1])
  }, numeric(5))
  expect_gt(mean(performance["found", ]), .8)
  expect_gt(mean(performance["sensitivity", ]), .8)
  expect_lt(mean(performance["false_inclusion", ]), .10)
  expect_gt(mean(performance["planted_outlier", ]), .8)
  expect_true(all(performance["shift_p", ] < .05))
})
