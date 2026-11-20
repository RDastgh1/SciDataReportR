# SOM grid geometry ------------------------------------------------------------

# kohonen places hexagonal and rectangular grids on unit spacing, so two nodes
# are neighbours when their grid centres are one unit apart.
.SOMGridAdjacency <- function(pts) {
  adjacent <- as.matrix(stats::dist(pts)) <= 1.01
  diag(adjacent) <- FALSE
  adjacent
}

# A node is on a cluster border when any neighbour belongs to another cluster.
.SOMBorderNodes <- function(clusters, adjacent) {
  vapply(seq_along(clusters), function(i) {
    neighbours <- clusters[adjacent[i, ]]
    any(!is.na(clusters[i]) & !is.na(neighbours) & neighbours != clusters[i])
  }, logical(1))
}

# Polygon vertices for every node: pointy-top hexagons (kohonen offsets
# alternate rows by half a node) or unit squares.
.SOMNodeShapes <- function(pts, topo = "hexagonal") {
  if (identical(topo, "hexagonal")) {
    angles <- pi / 6 + (0:5) * pi / 3
    radius <- 1 / sqrt(3)
  } else {
    angles <- pi / 4 + (0:3) * pi / 2
    radius <- sqrt(2) / 2
  }
  data.frame(
    NodeID = rep(seq_len(nrow(pts)), each = length(angles)),
    x = rep(pts[, 1], each = length(angles)) + radius * cos(angles),
    y = rep(pts[, 2], each = length(angles)) + radius * sin(angles))
}

# The shared edge between each pair of neighbouring nodes in different
# clusters: perpendicular to the line joining their centres, through its
# midpoint, one hexagon side (or one square side) long.
.SOMBorderSegments <- function(pts, clusters, adjacent, topo = "hexagonal") {
  pairs <- which(adjacent & upper.tri(adjacent), arr.ind = TRUE)
  differs <- !is.na(clusters[pairs[, 1]]) & !is.na(clusters[pairs[, 2]]) &
    clusters[pairs[, 1]] != clusters[pairs[, 2]]
  pairs <- pairs[differs, , drop = FALSE]
  if (!nrow(pairs)) {
    return(data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
  }
  from <- pts[pairs[, 1], , drop = FALSE]
  to <- pts[pairs[, 2], , drop = FALSE]
  middle <- (from + to) / 2
  direction <- (to - from) / sqrt(rowSums((to - from)^2))
  half <- if (identical(topo, "hexagonal")) 1 / (2 * sqrt(3)) else 0.5
  data.frame(
    x = middle[, 1] - direction[, 2] * half, y = middle[, 2] + direction[, 1] * half,
    xend = middle[, 1] + direction[, 2] * half, yend = middle[, 2] - direction[, 1] * half)
}

.SOMMapModelParts <- function(object, caller) {
  if (!is.list(object) || is.data.frame(object)) {
    stop(caller, "() needs a model from CreateClusterModel_SOM_MClust().", call. = FALSE)
  }
  som_grid <- object$ModelInfo_SOM$som_grid
  node_df <- object$ProbFit$node
  if (is.null(som_grid$pts) || !is.data.frame(node_df) ||
      !all(c("NodeID", "Cluster") %in% names(node_df))) {
    stop(caller, "() needs a model from CreateClusterModel_SOM_MClust().", call. = FALSE)
  }
  pts <- as.matrix(som_grid$pts)
  nodes <- seq_len(nrow(pts))
  node_df <- node_df[match(nodes, node_df$NodeID), , drop = FALSE]
  node_n <- if ("NodeN" %in% names(node_df)) node_df$NodeN else
    tabulate(object$ProbFit$individual$SOM_Node, length(nodes))
  list(pts = pts, topo = if (is.null(som_grid$topo)) "hexagonal" else som_grid$topo,
    node_df = node_df, clusters = node_df$Cluster, node_n = node_n,
    codes = as.matrix(object$ModelInfo_SOM$som_codes))
}

#' Plot the SOM of a SOM + Mclust model
#'
#' @description Draws the trained self-organizing map as its grid of nodes,
#' coloured by cluster, by the number of participants on each node, by the
#' share of boundary participants, by the U-matrix, or by one clustering
#' variable. Black lines mark cluster borders: the edges between neighbouring
#' nodes that belong to different clusters.
#' @param object A model from [CreateClusterModel_SOM_MClust()].
#' @param fill What to colour nodes by: `"Cluster"` (each node's Mclust
#'   cluster), `"People"` (training participants whose best-matching node it
#'   is), `"Boundary"` (the share of those participants flagged as boundary
#'   participants under soft assignment), `"UMatrix"` (the mean distance from
#'   the node's prototype to its neighbours' prototypes), or the name of a
#'   clustering variable (the node prototypes' values for that variable, on the
#'   scale the SOM was trained on).
#' @param borders Logical; draw cluster borders.
#' @param highlight_nodes Optional node IDs to outline, for example
#'   `ExplainSOMAssignment(model, row)$nodes$NodeID`.
#' @param title Plot title. `NULL` uses a title that names `fill`.
#' @details Node IDs follow kohonen's numbering, row by row from the bottom
#'   left of the map; they match `ProbFit$node$NodeID` and the participant
#'   column `SOM_Node`. Nodes with no training participants are grey in the
#'   boundary map. High U-matrix values mark places where neighbouring
#'   prototypes differ most; when they line up with cluster borders, the
#'   clusters follow natural gaps in the data rather than cutting a continuum
#'   (see `ModelInfo_SOM$diagnostics$border_alignment`).
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_SOM_MClust(
#'   df_Training, paste0("Var", 1:4), method = "finalize",
#'   final_model = "EEI", final_k = 3, som_xdim = 6, som_ydim = 6
#' )
#' PlotSOMMap(model)
#' PlotSOMMap(model, fill = "Boundary")
#' PlotSOMMap(model, fill = "Var1")
#' }
#' @export
PlotSOMMap <- function(object, fill = "Cluster", borders = TRUE,
    highlight_nodes = NULL, title = NULL) {
  parts <- .SOMMapModelParts(object, "PlotSOMMap")
  if (!is.character(fill) || length(fill) != 1L || is.na(fill)) {
    stop("fill must be a single string.", call. = FALSE)
  }
  n_nodes <- nrow(parts$pts)
  codes <- parts$codes
  variable_column <- NA_integer_
  if (!fill %in% c("Cluster", "People", "Boundary", "UMatrix")) {
    variable_column <- match(fill, colnames(codes))
    if (is.na(variable_column)) variable_column <- match(paste0("Z_", fill), colnames(codes))
    if (is.na(variable_column) && length(object$vars_used) == ncol(codes)) {
      variable_column <- match(fill, object$vars_used)
    }
    if (is.na(variable_column)) {
      stop("fill must be \"Cluster\", \"People\", \"Boundary\", \"UMatrix\", or a clustering variable: ",
        paste(object$vars_used, collapse = ", "), ".", call. = FALSE)
    }
  }

  adjacent <- .SOMGridAdjacency(parts$pts)
  df_Shapes <- .SOMNodeShapes(parts$pts, parts$topo)
  caption <- NULL
  if (fill == "Cluster") {
    values <- .ClusterFactor(parts$clusters)
    fill_scale <- ggplot2::scale_fill_manual(
      values = .SciDataNamedValues(levels(values), fixed = c(Noise = "grey70")),
      name = "Cluster", drop = FALSE)
    default_title <- "SOM nodes by cluster"
  } else if (fill == "People") {
    values <- parts$node_n
    fill_scale <- ggplot2::scale_fill_viridis_c(name = "Participants")
    default_title <- "Participants per SOM node"
  } else if (fill == "Boundary") {
    individual <- object$ProbFit$individual
    flagged <- !is.na(individual$SOM_Node) & !is.na(individual$Flag_BoundaryParticipant)
    values <- rep(NA_real_, n_nodes)
    shares <- tapply(as.numeric(individual$Flag_BoundaryParticipant[flagged]),
      individual$SOM_Node[flagged], mean)
    values[as.integer(names(shares))] <- as.numeric(shares)
    fill_scale <- ggplot2::scale_fill_viridis_c(name = "Boundary\nparticipants",
      limits = c(0, 1), labels = scales::label_percent(), na.value = "grey85")
    default_title <- "Boundary participants per SOM node"
    if (anyNA(values)) caption <- "Grey: no training participants on the node"
  } else if (fill == "UMatrix") {
    u_matrix <- object$ModelInfo_SOM$diagnostics$u_matrix
    values <- if (is.data.frame(u_matrix) && nrow(u_matrix) == n_nodes) {
      u_matrix$U[match(seq_len(n_nodes), u_matrix$NodeID)]
    } else {
      rowSums(as.matrix(stats::dist(codes)) * adjacent) / rowSums(adjacent)
    }
    fill_scale <- ggplot2::scale_fill_viridis_c(name = "Distance to\nneighbours")
    default_title <- "U-matrix: distance between neighbouring node prototypes"
  } else {
    values <- codes[, variable_column]
    limit <- max(abs(values), na.rm = TRUE)
    fill_scale <- ggplot2::scale_fill_gradient2(
      low = "#1864AB", mid = "white", high = "#C92A2A", midpoint = 0,
      limits = c(-limit, limit), name = "Node\nprototype")
    variable <- if (variable_column <= length(object$vars_used)) {
      object$vars_used[[variable_column]]
    } else fill
    default_title <- .ClusterVariableLabel(object$DataWithClusters, variable)
  }
  df_Shapes$Value <- values[df_Shapes$NodeID]

  plot <- ggplot2::ggplot(df_Shapes, ggplot2::aes(
      x = .data$x, y = .data$y, group = .data$NodeID, fill = .data$Value)) +
    ggplot2::geom_polygon(colour = "white", linewidth = 0.3) +
    fill_scale
  if (isTRUE(borders)) {
    df_Borders <- .SOMBorderSegments(parts$pts, parts$clusters, adjacent, parts$topo)
    if (nrow(df_Borders)) {
      plot <- plot + ggplot2::geom_segment(data = df_Borders,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        inherit.aes = FALSE, colour = "black", linewidth = 0.8, lineend = "round")
    }
  }
  if (length(highlight_nodes)) {
    highlight_nodes <- unique(as.integer(highlight_nodes))
    if (anyNA(highlight_nodes) || any(highlight_nodes < 1L | highlight_nodes > n_nodes)) {
      stop("highlight_nodes must be node IDs between 1 and ", n_nodes, ".", call. = FALSE)
    }
    # A black outline on a white halo shows on every fill and cannot be
    # mistaken for a cluster border or a cluster colour.
    df_Highlight <- df_Shapes[df_Shapes$NodeID %in% highlight_nodes, , drop = FALSE]
    highlight_aes <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$NodeID)
    plot <- plot +
      ggplot2::geom_polygon(data = df_Highlight, highlight_aes, inherit.aes = FALSE,
        fill = NA, colour = "white", linewidth = 2.4) +
      ggplot2::geom_polygon(data = df_Highlight, highlight_aes, inherit.aes = FALSE,
        fill = NA, colour = "black", linewidth = 1)
  }
  plot +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::labs(title = if (is.null(title)) default_title else title,
      caption = caption)
}

# Scale rows with the model's frozen training preprocessing, as ProjectCluster()
# does, and return them in the SOM's column order.
.SOMScaledRows <- function(object, data) {
  z_vars <- object$ZScoreVars
  if (is.null(z_vars)) z_vars <- paste0("Z_", object$vars_used)
  if (identical(object$ZScoreType, "PreZScored")) {
    missing_vars <- setdiff(z_vars, names(data))
    if (length(missing_vars)) {
      stop("data is missing required pre-zscored columns: ",
        paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    z <- data[, z_vars, drop = FALSE]
  } else {
    missing_vars <- setdiff(object$vars_used, names(data))
    if (length(missing_vars)) {
      stop("data is missing clustering variable(s): ",
        paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    z <- ProjectZScore(data = data, variables = object$vars_used,
      parameters = object$ZScoreObject, ParameterInputType = "ZScoreObj",
      names_prefix = "Z_")$ZScores[, z_vars, drop = FALSE]
  }
  z <- as.matrix(z)
  colnames(z) <- colnames(object$ModelInfo_SOM$som_codes)
  z
}

#' Explain participants' SOM + Mclust assignments
#'
#' @description Shows how a participant's cluster is reached: their
#' best-matching SOM node (hard assignment), the nearest nodes used by soft
#' assignment with their distances and weights, each of those nodes' Mclust
#' cluster probabilities, and the blended soft probabilities. Everything is
#' computed with the model's own settings, so the results reproduce the
#' model's `Cluster_Hard`, `Cluster_Soft`, `SoftProb_*` and boundary flags.
#' @param object A model from [CreateClusterModel_SOM_MClust()].
#' @param row Row number(s) to explain: rows of `data`, or of the training
#'   data (`object$DataWithClusters`) when `data` is `NULL`.
#' @param data Optional new data. Rows are scaled with the frozen training
#'   preprocessing, as [ProjectCluster()] does. `NULL` explains training
#'   participants.
#' @details Soft probabilities are a weighted average of the nearest nodes'
#'   Mclust probabilities, \eqn{p_c = \sum_i w_i q_{ic}}, where the weights
#'   \eqn{w_i} sum to one and come from the model's `soft_kernel`: inverse
#'   squared distance to the node prototype, or a Gaussian kernel with the
#'   bandwidth stored at training. This is fuzzy k-nearest-neighbour
#'   classification with SOM node prototypes as the neighbours. A participant
#'   is a boundary participant when the margin between their two largest
#'   probabilities is below `boundary_margin_threshold` or their normalized
#'   entropy is above `boundary_entropy_threshold`. Distances are Euclidean, on
#'   the scale the SOM was trained on.
#' @return An object of class `SOMAssignmentExplanation`: a list with
#'   \describe{
#'     \item{person}{One row per explained row: `Row`, best-matching node
#'       (`SOM_Node`, `SOM_Distance`, `NodeN` training participants on it),
#'       `Cluster_Hard`, `Cluster_Soft`, `SoftMargin`, `SoftEntropy`,
#'       `Flag_BoundaryParticipant`, `BoundaryReason`, `SoftProb_*`, and
#'       `Cluster` under the model's `assignment` rule.}
#'     \item{nodes}{One row per explained row and nearest node: `Row`, `Rank`,
#'       `NodeID`, `Distance`, `Weight`, the node's cluster (`NodeCluster`),
#'       `NodeN`, and the node's Mclust probabilities `Prob_*`.}
#'     \item{settings}{The assignment settings used.}
#'   }
#'   Printing it writes the explanation out. Pass `nodes$NodeID` to
#'   [PlotSOMMap()]'s `highlight_nodes` to show the nodes on the map.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_SOM_MClust(
#'   df_Training, paste0("Var", 1:4), method = "finalize",
#'   final_model = "EEI", final_k = 3, som_xdim = 6, som_ydim = 6
#' )
#' explanation <- ExplainSOMAssignment(model, row = 1)
#' explanation
#' PlotSOMMap(model, highlight_nodes = explanation$nodes$NodeID)
#' }
#' @export
ExplainSOMAssignment <- function(object, row, data = NULL) {
  parts <- .SOMMapModelParts(object, "ExplainSOMAssignment")
  prob_cols <- grep("^prob_[0-9]+$", names(parts$node_df), value = TRUE)
  if (!length(prob_cols) || is.null(object$ModelInfo_SOM$som_model)) {
    stop("ExplainSOMAssignment() needs a model from CreateClusterModel_SOM_MClust().",
      call. = FALSE)
  }
  source_data <- if (is.null(data)) object$DataWithClusters else data
  if (!is.data.frame(source_data)) stop("data must be a data frame.", call. = FALSE)
  if (missing(row) || !is.numeric(row) || !length(row) || any(!is.finite(row)) ||
      any(row != round(row)) || any(row < 1 | row > nrow(source_data))) {
    stop("row must contain row numbers between 1 and ", nrow(source_data), ".",
      call. = FALSE)
  }
  row <- as.integer(row)
  zmat <- .SOMScaledRows(object, source_data[row, , drop = FALSE])
  incomplete <- !stats::complete.cases(zmat)
  if (any(incomplete)) {
    stop("Row(s) ", paste(row[incomplete], collapse = ", "),
      " have missing clustering variables, so they are not assigned to a SOM node.",
      call. = FALSE)
  }

  spec <- object$Specification$projection
  Setting <- function(name, default) if (is.null(spec[[name]])) default else spec[[name]]
  settings <- list(
    assignment = Setting("assignment", "hard"),
    soft_n_nodes = as.integer(Setting("soft_n_nodes", 3L)),
    soft_kernel = Setting("soft_kernel", "inverse"),
    soft_bandwidth = spec$soft_bandwidth,
    boundary_margin_threshold = Setting("boundary_margin_threshold", 0.20),
    boundary_entropy_threshold = Setting("boundary_entropy_threshold", 0.60))
  soft <- .SOMSoftProjection(zmat, parts$codes, parts$node_df,
    n_nodes = settings$soft_n_nodes,
    margin_threshold = settings$boundary_margin_threshold,
    entropy_threshold = settings$boundary_entropy_threshold,
    soft_kernel = settings$soft_kernel, soft_bandwidth = settings$soft_bandwidth)
  mapping <- kohonen::map(object$ModelInfo_SOM$som_model, newdata = zmat)
  bmu <- mapping$unit.classif
  node_prob <- as.matrix(parts$node_df[, prob_cols, drop = FALSE])

  person <- tibble::tibble(
    Row = row, SOM_Node = bmu, SOM_Distance = sqrt(mapping$distances),
    NodeN = parts$node_n[bmu], Cluster_Hard = parts$clusters[bmu],
    Cluster_Soft = soft$cluster, SoftMargin = soft$margin,
    SoftEntropy = soft$entropy, Flag_BoundaryParticipant = soft$boundary,
    BoundaryReason = soft$reason)
  for (j in seq_along(prob_cols)) person[[paste0("SoftProb_", j)]] <- soft$probability[, j]
  person$Cluster <- if (identical(settings$assignment, "soft")) {
    person$Cluster_Soft
  } else person$Cluster_Hard

  nodes <- dplyr::bind_rows(lapply(seq_along(row), function(i) {
    ids <- soft$nearest[i, ]
    out <- tibble::tibble(Row = row[[i]], Rank = seq_along(ids), NodeID = ids,
      Distance = soft$distance[i, ], Weight = soft$weights[i, ],
      NodeCluster = parts$clusters[ids], NodeN = parts$node_n[ids])
    for (j in seq_along(prob_cols)) out[[paste0("Prob_", j)]] <- node_prob[ids, j]
    out
  }))
  structure(list(person = person, nodes = nodes, settings = settings),
    class = "SOMAssignmentExplanation")
}

#' @export
print.SOMAssignmentExplanation <- function(x, digits = 2, ...) {
  s <- x$settings
  weights <- if (identical(s$soft_kernel, "gaussian")) {
    sprintf("Gaussian kernel, bandwidth %.2f", s$soft_bandwidth)
  } else "inverse squared distance"
  cat(sprintf("SOM + Mclust assignment (%s assignment). Soft weights: %s, %d nearest nodes.\n",
    s$assignment, weights, s$soft_n_nodes))
  soft_cols <- grep("^SoftProb_", names(x$person), value = TRUE)
  for (i in seq_len(nrow(x$person))) {
    p <- x$person[i, ]
    cat(sprintf("\nRow %d: best-matching node %d (Cluster %d, %d training participant%s), distance %.2f\n",
      p$Row, p$SOM_Node, p$Cluster_Hard, p$NodeN, if (isTRUE(p$NodeN == 1)) "" else "s",
      p$SOM_Distance))
    nodes <- x$nodes[x$nodes$Row == p$Row, , drop = FALSE]
    prob_cols <- grep("^Prob_", names(nodes), value = TRUE)
    shown <- data.frame(Rank = nodes$Rank, Node = nodes$NodeID,
      Distance = round(nodes$Distance, digits), Weight = round(nodes$Weight, digits),
      Cluster = nodes$NodeCluster, Participants = nodes$NodeN)
    for (column in prob_cols) {
      shown[[paste0("P(", sub("^Prob_", "", column), ")")]] <- round(nodes[[column]], digits)
    }
    print(shown, row.names = FALSE)
    probabilities <- unlist(p[soft_cols])
    cat("Soft probabilities: ", paste(sprintf("Cluster %d %.2f", seq_along(probabilities),
      probabilities), collapse = ", "), "\n", sep = "")
    cat(sprintf("Hard cluster %d, soft cluster %d. Margin %.2f, entropy %.2f: %s.\n",
      p$Cluster_Hard, p$Cluster_Soft, p$SoftMargin, p$SoftEntropy,
      if (isTRUE(p$Flag_BoundaryParticipant)) {
        paste0("boundary participant (", gsub("_", " ", p$BoundaryReason), ")")
      } else "not a boundary participant"))
  }
  invisible(x)
}
