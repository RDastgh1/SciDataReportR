# Extracted from test-clustering-pipelines.R:89

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
assign("ProjectCluster.Pipeline_AliasTest",
    function(object, new_df, ...) "projected", envir = .GlobalEnv)
on.exit(rm(ProjectCluster.Pipeline_AliasTest, envir = .GlobalEnv), add = TRUE)
object <- structure(list(), class = "Pipeline_AliasTest")
expect_identical(Project_SOMClust(object, data.frame(x = 1)), "projected")
