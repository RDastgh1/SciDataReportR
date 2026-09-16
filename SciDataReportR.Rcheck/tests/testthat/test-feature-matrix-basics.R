test_that("AssessFeatureMissingness returns an auditable feature summary", {
  df_features <- data.frame(
    FeatureA = c(1, 2, NA, 4),
    FeatureB = c(2, 3, 4, 5),
    Group = c("A", "A", "B", "B")
  )

  result <- AssessFeatureMissingness(df_features, variables = c("FeatureA", "FeatureB"))

  expect_named(result, c("summary", "plot"))
  expect_equal(result$summary$PctMissing[result$summary$Variable == "FeatureA"], 25)
  expect_false(result$summary$Flagged[result$summary$Variable == "FeatureA"])
})

test_that("TransformFeatureMatrix retains unselected metadata and supports no transformation", {
  df_features <- data.frame(
    FeatureA = c(1, 2, 3),
    FeatureB = c(4, 5, 6),
    Group = c("A", "A", "B")
  )

  result <- TransformFeatureMatrix(
    df_features,
    variables = c("FeatureA", "FeatureB"),
    method = "none"
  )

  expect_equal(result$data, df_features)
  expect_equal(result$summary$Method, rep("none", 2))
})
