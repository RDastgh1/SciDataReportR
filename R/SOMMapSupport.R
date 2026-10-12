#' Review SOM map support before fitting a SOM + Mclust model
#'
#' @description Calculates the competing codebook-node support and participant
#'   occupancy constraints for a planned SOM + Mclust search. Mclust is fit to
#'   SOM codebook vectors, so its effective sample size is the number of nodes;
#'   participant count instead limits how sparse a map should become.
#'
#' @param n Number of complete training participants.
#' @param p Number of post-preprocessing clustering variables entering Mclust.
#' @param k_range Requested profile counts.
#' @param models TidyLPA Mclust model IDs (1, 2, 3, or 6).
#' @param min_nodes_per_parameter Target codebook nodes per estimated mixture
#'   parameter. Default `3`.
#' @param min_participants_per_node Target complete training participants per
#'   SOM node. Default `10`; this is an engineering map-health target, not a
#'   universal statistical threshold.
#'
#' @return A tibble with one row per model/k candidate and common map-sizing
#'   diagnostics. `SupportStatus` is `"supported"`, `"under-supported"`, or
#'   `"over-parameterized"` at the recommended automatic map size.
#' @examples
#' SOMMapSupport(n = 4330, p = 11, k_range = 2:8, models = 1:3)
#' @export
SOMMapSupport <- function(n, p, k_range, models,
    min_nodes_per_parameter = 3,
    min_participants_per_node = 10) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("n must be one positive finite number.")
  }
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 1) {
    stop("p must be one positive finite number.")
  }
  if (!is.numeric(k_range) || !length(k_range) || anyNA(k_range) ||
      any(k_range < 2) || any(k_range != as.integer(k_range))) {
    stop("k_range must contain whole-number profile counts of at least 2.")
  }
  models <- .ResolveMclustModels(models)
  if (!is.numeric(min_nodes_per_parameter) || length(min_nodes_per_parameter) != 1L ||
      !is.finite(min_nodes_per_parameter) || min_nodes_per_parameter <= 0) {
    stop("min_nodes_per_parameter must be one positive finite number.")
  }
  if (!is.numeric(min_participants_per_node) || length(min_participants_per_node) != 1L ||
      !is.finite(min_participants_per_node) || min_participants_per_node <= 0) {
    stop("min_participants_per_node must be one positive finite number.")
  }
  df_Support <- tidyr::expand_grid(Model = as.integer(models), Classes =
    sort(unique(as.integer(k_range)))) %>%
    dplyr::mutate(
      MixtureParameters = mapply(.SOMMclustParameterCount,
        .data$Model, as.integer(p), .data$Classes),
      ParameterFloorSide = ceiling(sqrt(min_nodes_per_parameter *
        .data$MixtureParameters)))
  base_side <- ceiling(sqrt(20 * max(df_Support$Classes)))
  floor_side <- max(df_Support$ParameterFloorSide)
  occupancy_ceiling_side <- floor(sqrt(n / min_participants_per_node))
  requested_side <- max(base_side, floor_side)
  recommended_side <- max(2L, min(requested_side,
    max(2L, occupancy_ceiling_side)))
  binding <- if (requested_side > occupancy_ceiling_side) {
    "participant_occupancy"
  } else if (floor_side > base_side) {
    "parameter_floor"
  } else {
    "k_informed"
  }
  .SOMSupportAtSide(df_Support, n = n, side = recommended_side,
    min_nodes_per_parameter = min_nodes_per_parameter,
    min_participants_per_node = min_participants_per_node) %>%
    dplyr::mutate(
      BaseSide = as.integer(base_side),
      FloorSide = as.integer(floor_side),
      OccupancyCeilingSide = as.integer(occupancy_ceiling_side),
      RecommendedSide = as.integer(recommended_side),
      BindingConstraint = binding,
      MinNodesPerParameter = min_nodes_per_parameter,
      MinParticipantsPerNode = min_participants_per_node)
}

.SOMSupportAtSide <- function(df_Support, n, side,
    min_nodes_per_parameter, min_participants_per_node) {
  .SOMSupportAtNodes(df_Support, n = n, nodes = as.integer(side)^2,
    min_nodes_per_parameter = min_nodes_per_parameter,
    min_participants_per_node = min_participants_per_node) %>%
    dplyr::mutate(SOMSide = as.integer(side))
}

.SOMSupportAtNodes <- function(df_Support, n, nodes,
    min_nodes_per_parameter, min_participants_per_node) {
  df_Support %>%
    dplyr::mutate(
      MixtureNodes = as.integer(nodes),
      ParticipantsPerNode = n / nodes,
      OccupancyCompliant = .data$ParticipantsPerNode >= min_participants_per_node,
      NodesPerParameter = .data$MixtureNodes / .data$MixtureParameters,
      Overparameterized = .data$NodesPerParameter < 1,
      SupportStatus = dplyr::case_when(
        .data$Overparameterized ~ "over-parameterized",
        .data$NodesPerParameter < min_nodes_per_parameter ~ "under-supported",
        TRUE ~ "supported"))
}
