test_that("soft and hard certainty use the correct probabilities and keep node sizes", {
  for (assignment in c("soft", "hard")) {
    model <- ProjectionTestModel(assignment)
    result <- suppressWarnings(suppressMessages(ProjectCluster(model, EmbeddedTestData("Projection"))))
    df_Individual <- result$ProjectionFit$individual
    expected <- if (assignment == "soft") apply(as.matrix(df_Individual[, grep("^SoftProb_", names(df_Individual)), drop = FALSE]), 1, max) else df_Individual$max_prob
    expect_equal(df_Individual$MembershipProbability, expected)
    expect_equal(df_Individual$NodeN, model$ProbFit$node$NodeN[df_Individual$SOM_Node])
    expect_equal(result$DataWithClusters$Flag_OnePersonNode, df_Individual$NodeN == 1)
    expect_identical(levels(df_Individual$FitClass), .SOMFitClassLevels)
    expect_false(result$NoveltyGroups$overall_excess$Excess)
    expect_false(any(result$NoveltyGroups$groups$Candidate))
  }
})

test_that("legacy LOO recovery and explicit training references work", {
  model <- ProjectionTestModel()
  reference <- suppressWarnings(suppressMessages(ProjectCluster(model, EmbeddedTestData("Projection"))))
  legacy <- model
  legacy$ProbFit$individual$SOM_Distance_LOO <- NULL
  recovered <- suppressWarnings(suppressMessages(ProjectCluster(legacy, EmbeddedTestData("Projection"))))
  expect_equal(recovered$ProjectionFit$individual$NoveltyP,
    reference$ProjectionFit$individual$NoveltyP)
  legacy$ZScoreObject$ZScores <- NULL
  legacy$DataWithClusters <- NULL
  legacy$ModelInfo_SOM$som_model$data <- NULL
  expect_message(suppressWarnings(ProjectCluster(legacy, EmbeddedTestData("Projection"))), "cannot be recovered")
  training <- suppressWarnings(suppressMessages(ProjectCluster(model, EmbeddedTestData("Projection"),
    distance_reference = "training")))
  df_Individual <- training$ProjectionFit$individual
  expect_equal(df_Individual$NoveltyP,
    .SOMConformalP(df_Individual$SOM_Distance, model$ProbFit$individual$SOM_Distance))
})

test_that("same-population unusual rates are closer to chance with LOO", {
  model <- ProjectionTestModel()
  set.seed(902)
  df_Test <- EmbeddedTestData("Projection")
  rates <- vapply(c("training", "loo"), function(reference) {
    projected <- suppressWarnings(suppressMessages(ProjectCluster(model, df_Test, distance_reference = reference)))
    mean(projected$ProjectionFit$individual$Flag_Unusual)
  }, numeric(1))
  expect_lt(abs(rates[["loo"]] - .05), abs(rates[["training"]] - .05))
})
