test_that("every finalizable clustering reference documents Stability field by field", {
  root <- testthat::test_path("..", "..")
  # R CMD check runs the tests from the installed package, where the source
  # man/ directory is not on disk. This check reads the generated .Rd files,
  # so it only applies when the tests are run from the package source.
  skip_if_not(dir.exists(file.path(root, "man")), "source man/ is unavailable")
  pages <- c(
    "CreateClusterModel_MClust.Rd",
    "CreateClusterModel_KMeans.Rd",
    "CreateClusterModel_PCA_MClust.Rd",
    "CreateClusterModel_PCA_KMeans.Rd",
    "CreateClusterModel_HDBSCAN.Rd",
    "CreateClusterModel_Gower_PAM.Rd",
    "CreateClusterModel_LatentClass.Rd",
    "CreateClusterModel_MCA_MClust.Rd",
    "CreateClusterModel_SOM_MClust.Rd",
    "CreateClusterModel_SOM_HDBSCAN.Rd"
  )
  fields <- c(
    "Stability output", "coassignment_limit", "noise_policy", "refit_scope",
    "resample_type", "resample_fraction", "subsample_without_replacement",
    "FowlkesMallows", "StabilityARI_Mean", "StabilityJaccard_Mean",
    "ReproducibilityScore", "participant_inclusion", "cluster_inclusion",
    "row_ids", "Hubert and Arabie", "Monti et al."
  )

  for (page in pages) {
    path <- file.path(root, "man", page)
    expect_true(file.exists(path), info = page)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    for (field in fields) {
      expect_match(text, field, fixed = TRUE, info = paste(page, field))
    }
  }
})

test_that("SOM plus Mclust reference documents every returned stability field", {
  path <- testthat::test_path("..", "..", "man", "CreateClusterModel_SOM_MClust.Rd")
  skip_if_not(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  fields <- c(
    "resamples", "seed", "refit_scope", "coassignment_limit", "noise_policy",
    "Model", "Classes", "Replicate", "Status", "ARI", "VI", "NMI",
    "FowlkesMallows", "Error", "Cluster", "Jaccard", "StabilitySuccessRate",
    "StabilityARI_P05", "StabilityJaccard_Min", "failures", "RowIndex",
    "SuccessfulRefits", "InclusionProbability", "MeanInclusion", "P05Inclusion",
    "MinInclusion", "coassignment", "matrix", "row_ids", "cluster_recovery",
    "partition_metrics", "cluster_inclusion"
  )
  for (field in fields) {
    expect_match(text, field, fixed = TRUE, info = field)
  }
})
