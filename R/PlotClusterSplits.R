#' Show how participants move between clusters as k increases
#'
#' @description
#' Answers "as clusters are added, is an existing cluster being split, or are
#' participants reshuffled?" For each model family, participants are followed
#' from the solution with k clusters to the solution with k + 1 clusters and
#' shown as an alluvial plot, with a transition table summarizing each step.
#'
#' @param object An exploratory model from [CreateClusterModel_SOM_MClust()],
#'   or a data frame whose columns are cluster assignments for the same
#'   participants at increasing k (column names such as `"k2"`, `"k3"`; one
#'   row per participant; treated as a single model family).
#' @param models Optional integer vector of model families to show (SOM +
#'   Mclust objects only). Defaults to every family with at least two
#'   consecutive k.
#' @param k_range Optional integer vector restricting the k shown.
#' @param title Plot title.
#'
#' @details
#' Cluster labels are arbitrary at every k, so each solution's clusters are
#' reordered to follow their largest parent in the previous solution; this
#' keeps a cluster that is carried forward in the same vertical position.
#'
#' For each step from k to k + 1 the transition table reports:
#' \itemize{
#'   \item `ParentPurity`: the share of participants whose cluster at k + 1
#'     came mostly from their own cluster at k (each cluster at k + 1 is
#'     credited to its largest parent). Values near 1 mean new clusters are
#'     carved out of existing ones (a hierarchical split); lower values mean
#'     participants are reshuffled between clusters.
#'   \item `SplitCluster`: the cluster at k whose members are divided most
#'     evenly among clusters at k + 1 (its largest child holds the smallest
#'     share of it), with `SplitRetained` giving that largest share.
#'   \item `SmallestNewCluster`: the smallest cluster at k + 1 as a share of
#'     participants.
#' }
#' With SOM + Mclust, assignments are the phenotypes of participants'
#' best-matching SOM nodes.
#'
#' @return A list with `plot` (a ggplot) and `transitions` (a tibble with one
#'   row per model and step).
#'
#' @examples
#' \donttest{
#' data("SampleData")
#' data("SampleVariableTypes")
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' model <- CreateClusterModel_SOM_MClust(
#'   df_Labelled, c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
#'   method = "exploratory", k_range = 2:4, models = 1,
#'   som_xdim = 5, som_ydim = 5, min_nodes_per_cluster = NULL,
#'   Relabel = FALSE)
#' splits <- PlotClusterSplits(model)
#' splits$plot
#' splits$transitions
#' }
#' @export
PlotClusterSplits <- function(object, models = NULL, k_range = NULL,
    title = "Where participants go as clusters are added") {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for PlotClusterSplits().")
  }
  assignments <- .ClusterSplitAssignments(object, models, k_range)
  if (!length(assignments)) {
    stop("No model family has cluster solutions at two or more consecutive k.")
  }

  transitions <- dplyr::bind_rows(lapply(names(assignments), function(model_name) {
    .ClusterSplitTransitions(assignments[[model_name]]) %>%
      dplyr::mutate(Model = model_name, .before = 1)
  }))

  flows <- dplyr::bind_rows(lapply(names(assignments), function(model_name) {
    aligned <- .AlignClusterLabels(assignments[[model_name]])
    path <- do.call(paste, c(aligned, sep = "-"))
    counts <- as.data.frame(table(path), stringsAsFactors = FALSE)
    lodes <- dplyr::bind_rows(lapply(seq_along(aligned), function(j) {
      dplyr::tibble(
        path = counts$path,
        k = as.integer(sub("^k", "", names(aligned)[j])),
        Cluster = vapply(strsplit(counts$path, "-", fixed = TRUE),
          function(parts) parts[[j]], character(1)),
        n = counts$Freq)
    }))
    lodes %>%
      dplyr::mutate(Model = model_name)
  }))
  flows$Cluster <- factor(flows$Cluster,
    levels = as.character(sort(unique(as.integer(flows$Cluster)))))
  flows$k <- factor(flows$k, levels = sort(unique(flows$k)))

  palette <- stats::setNames(
    .SciDataColorValues(nlevels(flows$Cluster)), levels(flows$Cluster))

  plot <- ggplot2::ggplot(flows, ggplot2::aes(
      x = .data$k, stratum = .data$Cluster, alluvium = .data$path, y = .data$n)) +
    ggalluvial::geom_flow(ggplot2::aes(fill = .data$Cluster),
      stat = "alluvium", lode.guidance = "frontback", aes.flow = "forward",
      alpha = 0.45) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = .data$Cluster), width = 0.35,
      colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::facet_wrap(~Model, ncol = 1, scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "Number of clusters (k)", y = "Participants",
      fill = "Cluster (aligned to previous k)")

  list(plot = plot, transitions = transitions)
}

# Participant-by-k assignment frames, one per model family.
.ClusterSplitAssignments <- function(object, models = NULL, k_range = NULL) {
  if (is.data.frame(object)) {
    out <- list(Clusters = as.data.frame(lapply(object, as.integer)))
    names(out$Clusters) <- names(object)
  } else {
    mclust_info <- object$ModelInfo_MClust
    if (is.null(mclust_info) || is.null(mclust_info$node_assignments)) {
      stop("object must be an exploratory CreateClusterModel_SOM_MClust() model or a data frame of assignments.")
    }
    individual <- identical(object$assignment, "individual") &&
      !is.null(mclust_info$individual_assignments)
    per_fit <- if (individual) {
      mapped <- !is.na(object$ProbFit$individual$SOM_Node)
      lapply(mclust_info$individual_assignments, function(x) x[mapped])
    } else {
      nodes <- object$ProbFit$individual$SOM_Node
      nodes <- nodes[!is.na(nodes)]
      lapply(mclust_info$node_assignments, function(x) x[nodes])
    }
    fit_model <- as.integer(sub("^model_([0-9]+)_class_.*$", "\\1", names(per_fit)))
    fit_k <- as.integer(sub("^.*_class_([0-9]+)$", "\\1", names(per_fit)))
    families <- sort(unique(fit_model))
    if (!is.null(models)) families <- intersect(families, models)
    out <- lapply(stats::setNames(families, families), function(m) {
      ks <- sort(fit_k[fit_model == m])
      if (!is.null(k_range)) ks <- intersect(ks, k_range)
      frame <- as.data.frame(lapply(ks, function(k) per_fit[[paste0("model_", m, "_class_", k)]]))
      names(frame) <- paste0("k", ks)
      frame
    })
    labels <- .MclustModelLabels[as.character(families)]
    names(out) <- ifelse(is.na(labels), paste("Model", families), labels)
  }
  # keep the longest run of consecutive k in each family
  out <- lapply(out, function(frame) {
    ks <- as.integer(sub("^k", "", names(frame)))
    breaks <- c(0, which(diff(ks) != 1), length(ks))
    runs <- lapply(seq_len(length(breaks) - 1), function(i) (breaks[i] + 1):breaks[i + 1])
    frame[, runs[[which.max(lengths(runs))]], drop = FALSE]
  })
  Filter(function(frame) ncol(frame) >= 2, out)
}

# Relabel each solution so a cluster carried forward keeps the label of its
# largest parent, and new clusters take the next free labels.
.AlignClusterLabels <- function(frame) {
  frame <- frame[stats::complete.cases(frame), , drop = FALSE]
  aligned <- frame
  first <- as.integer(factor(frame[[1]], levels = names(sort(table(frame[[1]]), decreasing = TRUE))))
  aligned[[1]] <- first
  for (j in seq_len(ncol(frame))[-1]) {
    tab <- table(child = frame[[j]], parent = aligned[[j - 1]])
    children <- rownames(tab)
    parent_of <- colnames(tab)[apply(tab, 1, which.max)]
    child_size <- rowSums(tab)
    label <- stats::setNames(rep(NA_integer_, length(children)), children)
    # each parent label goes to its largest child; remaining children get new labels
    for (p in unique(parent_of)) {
      candidates <- children[parent_of == p]
      label[candidates[which.max(child_size[candidates])]] <- as.integer(p)
    }
    unused <- setdiff(seq_along(children), label)
    label[is.na(label)] <- unused[seq_len(sum(is.na(label)))]
    aligned[[j]] <- unname(label[as.character(frame[[j]])])
  }
  aligned
}

.ClusterSplitTransitions <- function(frame) {
  frame <- frame[stats::complete.cases(frame), , drop = FALSE]
  aligned <- .AlignClusterLabels(frame)
  dplyr::bind_rows(lapply(seq_len(ncol(frame) - 1), function(j) {
    parent <- aligned[[j]]
    child <- aligned[[j + 1]]
    tab <- table(parent, child)
    child_n <- colSums(tab)
    purity <- sum(apply(tab, 2, max)) / sum(tab)
    retained <- apply(tab, 1, max) / rowSums(tab)
    split <- names(which.min(retained))
    dplyr::tibble(
      From = as.integer(sub("^k", "", names(frame)[j])),
      To = as.integer(sub("^k", "", names(frame)[j + 1])),
      ParentPurity = purity,
      SplitCluster = as.integer(split),
      SplitRetained = unname(retained[split]),
      SmallestNewCluster = min(child_n) / sum(child_n)
    )
  }))
}
