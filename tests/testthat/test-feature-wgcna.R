test_that("BuildFeatureWGCNA validates selected wide features", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  set.seed(1)
  df_data <- data.frame(SampleID = paste0("S", seq_len(30)), matrix(rnorm(30 * 9), ncol = 9), check.names = FALSE)
  names(df_data)[-1] <- paste0("Protein", seq_len(9))
  obj <- BuildFeatureWGCNA(df_data, variables = names(df_data)[-1], sample_id = "SampleID", soft_power = 4, min_module_size = 3, seed = 1)
  expect_s3_class(obj, "feature_wgcna_obj")
  expect_equal(nrow(obj$modules), 9)
  expect_equal(attr(obj$datExpr, "sample_ids"), df_data$SampleID)
  expect_error(BuildFeatureWGCNA(df_data, c("Protein1", "Protein1", "Protein2")), "unique")
  df_data$Protein1 <- 1
  expect_error(BuildFeatureWGCNA(df_data, names(df_data)[-1], soft_power = 4, min_module_size = 3), "constant")
})

test_that("AnalyzeFeatureModuleTraits requires an exact sample-ID join", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  set.seed(2)
  df_data <- data.frame(SampleID = paste0("S", seq_len(30)), matrix(rnorm(30 * 9), ncol = 9), check.names = FALSE)
  names(df_data)[-1] <- paste0("Protein", seq_len(9))
  obj <- BuildFeatureWGCNA(df_data, names(df_data)[-1], sample_id = "SampleID", soft_power = 4, min_module_size = 3)
  df_traits <- data.frame(SampleID = rev(df_data$SampleID), Age = seq_len(30), Group = rep(c("A", "B"), 15))
  result <- AnalyzeFeatureModuleTraits(obj, df_traits, sample_id = "SampleID", outcome = c("Age", "Group"))
  expect_s3_class(result, "feature_module_trait_obj")
  expect_true(all(result$join_audit$Matched))
  df_traits$SampleID[1] <- "other"
  expect_error(AnalyzeFeatureModuleTraits(obj, df_traits, sample_id = "SampleID", outcome = "Age"), "incomplete")
})
