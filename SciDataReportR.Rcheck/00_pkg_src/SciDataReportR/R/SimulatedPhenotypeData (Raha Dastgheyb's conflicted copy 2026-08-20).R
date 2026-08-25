#' Neutral simulated clustering and phenotyping benchmark
#'
#' A fixed 480-participant benchmark for evaluating train-once/project-many
#' clustering workflows. Twelve numeric variables support at least four retained
#' principal components and four balanced truth groups. Three four-level
#' categorical variables provide strongly separated, class-dependent response
#' patterns for categorical-method demonstrations, while `DensityX` and
#' `DensityY` provide two unequal-density groups plus noise.
#' `TruthCluster` and `TruthDensityGroup` are retained solely for teaching and
#' test evaluation; they must not be used as clustering features.
#'
#' @format `SimulatedPhenotypeData` has 480 rows: 320 training participants and
#'   160 projection participants. It contains `Var1` through `Var12`, `CatVar1`
#'   through `CatVar3`, density coordinates, truth fields, and the fixed cohort
#'   split. `SimulatedPhenotypeVariableTypes` is its complete labelled codebook.
#' @section Why the benchmark works:
#' The twelve numeric variables are organized in blocks, and each cluster is
#' high on a different block. That is the structure a clustering method is
#' supposed to recover, and the reason a method that finds nothing here is
#' genuinely failing rather than facing an impossible problem.
#'
#' `DensityX` and `DensityY` are a separate and deliberately harder problem:
#' two groups of unequal density plus a noise group, for methods that are meant
#' to detect outliers rather than partition everything.
#'
#' `TruthCluster` and `TruthDensityGroup` exist only so results can be scored.
#' They must never be used as clustering features.
#'
#' @seealso [SampleData] for the labelling and reporting examples, and
#'   [CreateClusterModel_SOM_MClust()] or [CreateClusterModel_PCA_KMeans()] for the
#'   pipelines this benchmark was built for.
#'
#' @examples
#' data(SimulatedPhenotypeData)
#'
#' # 480 participants, with the truth clusters balanced across both cohorts
#' dim(SimulatedPhenotypeData)
#' table(SimulatedPhenotypeData$Cohort, SimulatedPhenotypeData$TruthCluster)
#'
#' \donttest{
#' data(SimulatedPhenotypeVariableTypes)
#'
#' ShowTable <- function(x, caption = NULL, height = NULL) {
#'   htmltools::browsable(htmltools::HTML(as.character(
#'     kableExtra::kable_styling(
#'       knitr::kable(x, format = "html", caption = caption, row.names = FALSE),
#'       bootstrap_options = c("striped", "hover", "condensed"),
#'       full_width = FALSE
#'     )
#'   )))
#' }
#'
#' # A slice: identifiers, clustering variables, and the truth labels
#' ShowTable(
#'   utils::head(
#'     SimulatedPhenotypeData[, c(
#'       "ParticipantID", "Cohort", "TruthCluster",
#'       paste0("Var", 1:6), "CatVar1", "DensityX", "DensityY"
#'     )],
#'     8
#'   ),
#'   "First eight participants"
#' )
#'
#' # Cluster means, showing the block structure
#' df_Training <- SimulatedPhenotypeData[
#'   SimulatedPhenotypeData$Cohort == "Training",
#' ]
#' df_ClusterMeans <- stats::aggregate(
#'   df_Training[paste0("Var", 1:12)],
#'   by = list(Cluster = df_Training$TruthCluster),
#'   FUN = mean
#' )
#' df_ClusterMeans[-1] <- round(df_ClusterMeans[-1], 2)
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_ClusterMeans, full_width = TRUE)
#' )))
#'
#' # The density variables: two unequal groups plus noise
#' table(SimulatedPhenotypeData$TruthDensityGroup)
#'
#' ggplot2::ggplot(
#'   SimulatedPhenotypeData,
#'   ggplot2::aes(
#'     x = DensityX, y = DensityY, colour = TruthDensityGroup
#'   )
#' ) +
#'   ggplot2::geom_point(alpha = 0.7, size = 1.6) +
#'   ggplot2::labs(
#'     title = "Density groups, with noise",
#'     colour = "Truth group"
#'   ) +
#'   ggplot2::theme_bw()
#'
#' # The companion codebook
#' ShowTable(SimulatedPhenotypeVariableTypes, "SimulatedPhenotypeVariableTypes")
#' }
#'
#' @name SimulatedPhenotypeData
#' @aliases SimulatedPhenotypeVariableTypes
#' @usage data(SimulatedPhenotypeData)
#' @docType data
NULL
