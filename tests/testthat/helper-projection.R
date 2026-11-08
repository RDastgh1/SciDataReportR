EmbeddedTestData <- function(cohort = "Training", variables = paste0("Var", 1:3)) {
  data("SimulatedPhenotypeData", package = "SciDataReportR", envir = environment())
  SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == cohort, variables, drop = FALSE]
}

EmbeddedAB <- function(cohort = "Training") {
  df_Data <- EmbeddedTestData(cohort, c("Var1", "Var2"))
  names(df_Data) <- c("A", "B")
  df_Data
}

ProjectionTestModel <- function(assignment = "soft") {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  df_Training <- EmbeddedTestData()
  suppressWarnings(suppressMessages(CreateClusterModel_SOM_MClust(df_Training,
    names(df_Training), method = "finalize", final_model = 1, final_k = 4,
    som_xdim = 5, som_ydim = 5, assignment = assignment,
    mixture_bootstrap_n = 0, Relabel = FALSE)))
}
