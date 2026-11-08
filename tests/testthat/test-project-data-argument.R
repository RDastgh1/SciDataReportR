CheckProjectionAlias <- function(model, df_New) {
  expected <- suppressWarnings(suppressMessages(ProjectCluster(model, data = df_New)))
  positional <- suppressWarnings(suppressMessages(ProjectCluster(model, df_New)))
  expect_equal(positional$DataWithClusters, expected$DataWithClusters)
  expect_warning(actual <- withCallingHandlers(suppressMessages(ProjectCluster(model, new_df = df_New)),
    warning = function(condition) {
      if (!inherits(condition, "lifecycle_warning_deprecated")) invokeRestart("muffleWarning")
    }), class = "lifecycle_warning_deprecated")
  expect_equal(actual$DataWithClusters, expected$DataWithClusters)
  expect_equal(actual$ProbFit$individual, expected$ProbFit$individual)
}

test_that("data and new_df match across numeric, PCA and distance methods", {
  rlang::local_options(lifecycle_verbosity = "warning")
  skip_if_not_installed("mclust")
  set.seed(122)
  df_Test <- EmbeddedTestData()
  df_Projection <- EmbeddedTestData("Projection")
  for (constructor in list(CreateClusterModel_MClust, CreateClusterModel_PCA_MClust)) {
    model <- suppressWarnings(constructor(df_Test, names(df_Test), method = "finalize", final_k = 2, final_model = 1))
    CheckProjectionAlias(model, df_Projection)
  }
  for (constructor in list(CreateClusterModel_KMeans, CreateClusterModel_PCA_KMeans, CreateClusterModel_Gower_PAM)) {
    model <- suppressWarnings(constructor(df_Test, names(df_Test), method = "finalize", final_k = 2))
    CheckProjectionAlias(model, df_Projection)
  }
  if (requireNamespace("dbscan", quietly = TRUE)) {
    model <- suppressWarnings(CreateClusterModel_HDBSCAN(df_Test, names(df_Test), method = "finalize", final_minPts = 5, final_cluster_selection_epsilon = 0))
    CheckProjectionAlias(model, df_Projection)
  }
})

test_that("data and new_df match across categorical methods", {
  rlang::local_options(lifecycle_verbosity = "warning")
  skip_if_not_installed("poLCA")
  skip_if_not_installed("FactoMineR")
  set.seed(123)
  df_Test <- EmbeddedTestData(variables = paste0("CatVar", 1:3))
  df_Projection <- EmbeddedTestData("Projection", paste0("CatVar", 1:3))
  model <- suppressWarnings(CreateClusterModel_LatentClass(df_Test, names(df_Test), method = "finalize", final_k = 2, nrep = 1))
  CheckProjectionAlias(model, df_Projection)
  model <- suppressWarnings(CreateClusterModel_MCA_MClust(df_Test, names(df_Test), method = "finalize", final_k = 2, final_model = 1))
  CheckProjectionAlias(model, df_Projection)
})

test_that("data and new_df match across SOM methods and rejection methods", {
  rlang::local_options(lifecycle_verbosity = "warning")
  model <- ProjectionTestModel()
  CheckProjectionAlias(model, EmbeddedTestData("Projection"))
  if (requireNamespace("dbscan", quietly = TRUE)) {
    df_Test <- model$DataWithClusters[, model$vars_used]
    density <- suppressWarnings(suppressMessages(CreateClusterModel_SOM_HDBSCAN(df_Test,
      names(df_Test), method = "finalize", final_minPts = 3,
      final_cluster_selection_epsilon = 0, som_xdim = 4, som_ydim = 4)))
    CheckProjectionAlias(density, EmbeddedTestData("Projection"))
  }
  for (class_name in c("Pipeline_SOM_HDBSCAN_Exploratory", "unsupported")) {
    object <- structure(list(), class = class_name)
    expect_error(ProjectCluster(object, data = data.frame()), "requires|not a supported")
    expect_warning(expect_error(ProjectCluster(object, new_df = data.frame()),
      "requires|not a supported"), class = "lifecycle_warning_deprecated")
  }
})
