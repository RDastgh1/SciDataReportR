test_that("clustering phenotyping vignette documents every implemented workflow", {
  vignette_path <- testthat::test_path("..", "..", "vignettes", "ClusteringPhenotyping.qmd")
  testthat::skip_if_not(file.exists(vignette_path))

  vignette_text <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")
  expect_match(vignette_text, "CreateSOMClusterModel")
  expect_match(vignette_text, "CreateMclustClusterModel")
  expect_match(vignette_text, "CreateKMeansClusterModel")
  expect_match(vignette_text, "CreateHDBSCANClusterModel")
  expect_match(vignette_text, "CreateGowerPAMClusterModel")
  expect_match(vignette_text, "CreateLatentClassClusterModel")
  expect_match(vignette_text, "CreateMCAMclustClusterModel")
})
