test_that("simulated phenotype benchmark has stable training and projection structure", {
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  data("SimulatedPhenotypeVariableTypes", package = "SciDataReportR")
  expect_equal(nrow(SimulatedPhenotypeData), 480)
  expect_equal(table(SimulatedPhenotypeData$Cohort)[["Training"]], 320)
  expect_equal(table(SimulatedPhenotypeData$Cohort)[["Projection"]], 160)
  expect_equal(
    unname(table(
      SimulatedPhenotypeData$TruthCluster,
      SimulatedPhenotypeData$Cohort)[, "Training"]),
    rep(80L, 4))
  expect_equal(
    unname(table(
      SimulatedPhenotypeData$TruthCluster,
      SimulatedPhenotypeData$Cohort)[, "Projection"]),
    rep(40L, 4))
  expect_true(all(names(SimulatedPhenotypeData) %in% SimulatedPhenotypeVariableTypes$Variable))
  expect_equal(sum(is.na(SimulatedPhenotypeData$Var12)), 5)
  expect_true(all(
    SimulatedPhenotypeData$Cohort[is.na(SimulatedPhenotypeData$Var12)] ==
      "Projection"))
  expect_true(all(c("DenseA", "DenseB", "Noise") %in% levels(SimulatedPhenotypeData$TruthDensityGroup)))
  expect_equal(attr(SimulatedPhenotypeData$Var1, "label"), "Var 1")
  expect_equal(attr(SimulatedPhenotypeData$CatVar1, "label"), "CatVar 1")
})

test_that("neutral benchmark retains at least four PCA components", {
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  model <- CreateClusterModel_PCA_KMeans(
    df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4,
    nstart = 5)
  expect_gte(ncol(model$PCA$Scores), 4)
})

test_that("categorical benchmark supports truth-known LCA and MCA examples", {
  skip_if_not_installed("poLCA")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  df_Projection <- subset(SimulatedPhenotypeData, Cohort == "Projection")
  vars_Categorical <- paste0("CatVar", 1:3)

  model_LCA <- CreateClusterModel_LatentClass(
    df_Training, vars_Categorical, method = "finalize", final_k = 4,
    nrep = 5)
  projected_LCA <- ProjectCluster(model_LCA, df_Projection)
  model_MCA <- CreateClusterModel_MCA_MClust(
    df_Training, vars_Categorical, method = "finalize", final_k = 4,
    final_model = 1)
  projected_MCA <- ProjectCluster(model_MCA, df_Projection)

  expect_gte(mclust::adjustedRandIndex(
    df_Training$TruthCluster,
    model_LCA$ProbFit$individual$Cluster), 0.9)
  expect_gte(mclust::adjustedRandIndex(
    df_Projection$TruthCluster,
    projected_LCA$ProbFit$individual$Cluster), 0.8)
  expect_gte(mclust::adjustedRandIndex(
    df_Training$TruthCluster,
    model_MCA$ProbFit$individual$Cluster), 0.9)
  expect_gte(mclust::adjustedRandIndex(
    df_Projection$TruthCluster,
    projected_MCA$ProbFit$individual$Cluster), 0.9)
})
