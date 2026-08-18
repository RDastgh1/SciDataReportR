# Extracted from test-clustering-pipelines.R:81

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
DefaultRange <- function(fun, argument) eval(formals(fun)[[argument]])
expect_identical(DefaultRange(CreateMclustClusterModel, "k_range"), 2:10)
expect_identical(DefaultRange(CreateKMeansClusterModel, "k_range"), 2:10)
expect_identical(DefaultRange(CreatePCAMclustClusterModel, "k_range"), 2:10)
expect_identical(DefaultRange(CreatePCAKMeansClusterModel, "k_range"), 2:10)
expect_identical(DefaultRange(CreateGowerPAMClusterModel, "k_range"), 2:10)
expect_identical(DefaultRange(CreateLatentClassClusterModel, "k_range"), 2:10)
