# Extracted from test-clustering-pipelines.R:146

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    Sex = factor(c(rep("Male", 9), "Female", rep("Male", 2), rep("Female", 8))),
    Cluster = rep(c(1L, 2L), each = 10L))
plot <- PlotClusterComposition(df_Test, "Sex", df_Test$Cluster,
    style = "enrichment")
