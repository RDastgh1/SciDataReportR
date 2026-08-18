test_that("ProjectCluster reference documents the SOM projection contract", {
  path <- testthat::test_path("..", "..", "man", "ProjectCluster.Rd")
  skip_if_not(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  fields <- c(
    "SOM + Mclust projection results", "ProbFit$individual",
    "ProjectionFit$individual", "ProjectionFit$summary",
    "ProjectionFit$by_cluster", "ProjectionFit$out_of_support",
    "node_occupancy_js_divergence", "Jensen--Shannon",
    "essentially constant"
  )
  for (field in fields) {
    expect_match(text, field, fixed = TRUE, info = field)
  }
})
