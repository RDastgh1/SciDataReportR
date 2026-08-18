# Extracted from test-scidata-colors.R:375

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(SimulatedPhenotypeData, envir = environment())
training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
numeric_vars <- paste0("Var", 1:12)
categorical_vars <- paste0("CatVar", 1:3)
BuildEveryPlot <- function(x, path = "") {
    if (inherits(x, "ggplot")) {
      err <- tryCatch(
        {
          invisible(ggplot2::ggplot_build(x))
          NULL
        },
        error = function(e) paste0(path, ": ", conditionMessage(e))
      )
      return(err)
    }
    if (is.list(x)) {
      return(unlist(lapply(
        names(x), function(n) BuildEveryPlot(x[[n]], paste0(path, "$", n))
      )))
    }
    NULL
  }
models <- list(
    Mclust = function() {
      CreateMclustClusterModel(
        training, numeric_vars, method = "finalize",
        final_k = 4, final_model = "EEI")
    },
    KMeans = function() {
      CreateKMeansClusterModel(
        training, numeric_vars, method = "finalize", final_k = 4, nstart = 10)
    },
    HDBSCAN = function() {
      CreateHDBSCANClusterModel(
        training, c("DensityX", "DensityY"), method = "finalize",
        final_minPts = 10, final_cluster_selection_epsilon = 0)
    },
    LatentClass = function() {
      CreateLatentClassClusterModel(
        training, categorical_vars, method = "finalize", final_k = 3)
    }
  )
