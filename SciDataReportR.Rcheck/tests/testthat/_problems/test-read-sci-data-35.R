# Extracted from test-read-sci-data.R:35

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
path <- tempfile(fileext = ".tsv")
writeLines(
    c(
      "id\tvalue",
      "1\t2",
      "2\t3"
    ),
    path
  )
df <- ReadSciData(path, inspect = FALSE, fast_delimited = FALSE)
