# Extracted from test-clustering-pipelines.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("SimulatedPhenotypeData")
df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")
by_variable <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster)
by_cluster <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster,
    facet_by = "cluster")
enrichment <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster,
    style = "enrichment")
