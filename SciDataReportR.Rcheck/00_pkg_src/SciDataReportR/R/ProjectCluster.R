#' Project cases through a fitted clustering model
#'
#' @description Projects `new_df` through the frozen preprocessing, reduction,
#' and clustering layers stored in `object`. Dispatch is determined by the
#' fitted model's `Pipeline_*` class; callers do not select a method-specific
#' projector.
#'
#' @param object A finalized object returned by `CreateClusterModel_*()`.
#' @param new_df New data to project into the frozen cluster structure.
#' @param ... Method-specific projection options. SOM + Mclust accepts its
#'   frozen-distance and posterior-probability threshold overrides.
#' @return A method-specific projection result with the common `ProjectionFit`
#'   contract and `ProbFit` assignment table.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_New <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:6), method = "finalize", final_k = 3
#' )
#' projected <- ProjectCluster(model, df_New)
#' projected$ProjectionFit$summary
#' }
#' @export
ProjectCluster <- function(object, new_df, ...) {
  UseMethod("ProjectCluster")
}

#' @export
ProjectCluster.default <- function(object, new_df, ...) {
  stop(
    "object is not a supported finalized clustering model. Expected a Pipeline_* object returned by CreateClusterModel_*().",
    call. = FALSE
  )
}
