# Novelty scoring for SOM projections ---------------------------------------
#
# Three levels, kept separate on purpose:
#   * person: a conformal p-value for "is this person farther from the map than
#     training people are?", with Benjamini-Hochberg control across the cohort;
#   * cohort: does the projected cohort sit farther from the map, or on
#     different parts of it, than training did?
#   * group: do unusual people concentrate on one part of the map and resemble
#     each other more than training's own unusual people do?
# Membership probability is deliberately not used here: a low probability means
# "between known clusters", which is a different question from "outside them".

.SOMFitClassLevels <- c("Good fit", "Uncertain membership", "Outlier",
  "Member of candidate new group")
.SOMHardFitClassLevels <- c("Good fit", "Outlier", "Member of candidate new group")

# Conformal p-value of each distance against a reference sample:
# (1 + #reference >= d) / (n + 1). Valid without distributional assumptions
# when the new person is exchangeable with the reference people.
.SOMConformalP <- function(distance, reference) {
  reference <- sort(reference[is.finite(reference)])
  n <- length(reference)
  p <- rep(NA_real_, length(distance))
  ok <- is.finite(distance)
  if (!n) return(p)
  p[ok] <- (1 + n - findInterval(distance[ok], reference, left.open = TRUE)) / (n + 1)
  p
}

# The same p-value for each training person against the other training people.
.SOMTrainingConformalP <- function(distance) {
  ok <- is.finite(distance)
  reference <- sort(distance[ok])
  n <- length(reference)
  p <- rep(NA_real_, length(distance))
  if (n < 2L) return(p)
  p[ok] <- (n - findInterval(distance[ok], reference, left.open = TRUE)) / n
  p
}

# Person-level novelty columns from conformal p-values.
.SOMNoveltyColumns <- function(p, novelty_fdr, group_screen_alpha) {
  q <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  q[ok] <- stats::p.adjust(p[ok], method = "BH")
  list(NoveltyP = p, NoveltyQ = q,
    Flag_Outlier = ifelse(ok, q <= novelty_fdr, NA),
    Flag_Unusual = ifelse(ok, p <= group_screen_alpha, NA))
}

# Attainable chance rate of p <= alpha for a reference of size n.
.SOMChanceRate <- function(alpha, n) {
  if (!is.finite(n) || n < 1) return(alpha)
  floor(alpha * (n + 1)) / (n + 1)
}

.SOMFitClass <- function(has_distance, outlier, group_member = FALSE, clear = NULL) {
  n <- length(has_distance)
  outlier <- rep_len(if (is.null(outlier)) FALSE else outlier, n) %in% TRUE
  group_member <- rep_len(if (is.null(group_member)) FALSE else group_member, n) %in% TRUE
  if (is.null(clear)) {
    out <- ifelse(group_member, "Member of candidate new group",
      ifelse(outlier, "Outlier", "Good fit"))
    levels <- .SOMHardFitClassLevels
  } else {
    out <- ifelse(group_member, "Member of candidate new group",
      ifelse(outlier, "Outlier", ifelse(!is.na(clear) & clear, "Good fit", "Uncertain membership")))
    levels <- .SOMFitClassLevels
  }
  out[!has_distance] <- NA_character_
  factor(out, levels = levels)
}

.SOMMeanPairwiseDistance <- function(x) {
  if (nrow(x) < 2L) return(NA_real_)
  mean(stats::dist(x))
}

# Is a group of unusual people more alike than the same number of training
# people who were unusual in the same sense? Permutation p-value, one-sided
# for a tighter group.
.SOMGroupCoherence <- function(group_z, pool_z, permutations = 999L, seed = 1L) {
  observed <- .SOMMeanPairwiseDistance(group_z)
  size <- nrow(group_z)
  if (is.na(observed) || is.null(pool_z) || nrow(pool_z) < 2L) {
    return(list(Spread = observed, ReferenceSpread = NA_real_, P = NA_real_))
  }
  replace <- nrow(pool_z) < size
  reference <- .WithClusterSeed(seed, vapply(seq_len(permutations), function(i) {
    .SOMMeanPairwiseDistance(pool_z[sample.int(nrow(pool_z), size, replace = replace), , drop = FALSE])
  }, numeric(1)))
  list(Spread = observed, ReferenceSpread = stats::median(reference),
    P = (1 + sum(reference <= observed)) / (permutations + 1))
}

.WithClusterSeed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# Chi-square goodness of fit of projected counts to training proportions,
# Monte Carlo p-value because node counts are sparse. A half-count prior on
# training keeps nodes empty in training from having probability zero.
.SOMOccupancyTest <- function(training_counts, projected_counts, B = 2000L, seed = 1L) {
  if (sum(projected_counts) < 1 || length(training_counts) < 2L) {
    return(list(statistic = NA_real_, p = NA_real_))
  }
  probability <- (training_counts + 0.5) / sum(training_counts + 0.5)
  test <- .WithClusterSeed(seed, suppressWarnings(stats::chisq.test(projected_counts,
    p = probability, simulate.p.value = TRUE, B = B)))
  list(statistic = unname(test$statistic), p = unname(test$p.value))
}

# Cohort-level shift between the training reference and the projected cohort.
.SOMCohortShift <- function(training_distance, projected_distance, projected_p,
    group_screen_alpha, training_nodes, projected_nodes, n_nodes,
    training_clusters, projected_clusters, seed = 1L) {
  training_distance <- training_distance[is.finite(training_distance)]
  projected_distance <- projected_distance[is.finite(projected_distance)]
  projected_p <- projected_p[!is.na(projected_p)]
  n_train <- length(training_distance)
  n_proj <- length(projected_distance)
  rank_test <- if (n_train && n_proj) suppressWarnings(stats::wilcox.test(projected_distance,
    training_distance, alternative = "greater", exact = FALSE)) else NULL
  ks <- if (n_train && n_proj) suppressWarnings(stats::ks.test(projected_distance, training_distance)) else NULL
  chance <- .SOMChanceRate(group_screen_alpha, n_train)
  unusual <- sum(projected_p <= group_screen_alpha)
  excess_p <- if (length(projected_p)) stats::pbinom(unusual - 1L, length(projected_p), chance,
    lower.tail = FALSE) else NA_real_
  occupancy <- .SOMOccupancyTest(tabulate(training_nodes, n_nodes),
    tabulate(projected_nodes, n_nodes), seed = seed)
  cluster_levels <- sort(unique(stats::na.omit(c(training_clusters, projected_clusters))))
  prevalence <- .SOMOccupancyTest(
    tabulate(match(training_clusters, cluster_levels), length(cluster_levels)),
    tabulate(match(projected_clusters, cluster_levels), length(cluster_levels)), seed = seed)
  dplyr::tibble(
    Test = c("Distance shift", "Distance distribution", "Unusual-profile excess",
      "Node occupancy", "Cluster prevalence"),
    Question = c(
      "Are projected people farther from the map than training people?",
      "Does the distance distribution differ in shape or location?",
      paste0("Are more people unusual (p <= ", group_screen_alpha, ") than expected by chance?"),
      "Do projected people land on the map in the training proportions?",
      "Are clusters as common as in training?"),
    Statistic = c(
      if (is.null(rank_test)) NA_real_ else unname(rank_test$statistic) / (n_train * n_proj),
      if (is.null(ks)) NA_real_ else unname(ks$statistic),
      if (length(projected_p)) unusual / length(projected_p) else NA_real_,
      occupancy$statistic, prevalence$statistic),
    StatisticName = c("P(projected > training)", "Kolmogorov-Smirnov D",
      "Observed share", "Chi-square", "Chi-square"),
    Expected = c(0.5, 0, chance, NA_real_, NA_real_),
    P = c(if (is.null(rank_test)) NA_real_ else rank_test$p.value,
      if (is.null(ks)) NA_real_ else ks$p.value, excess_p, occupancy$p, prevalence$p))
}

# Group level: neighbourhoods of the map where unusual projected people pile up
# beyond that neighbourhood's own training rate, merged into groups, each
# checked for internal similarity against training's unusual people.
.SOMNoveltyGroups <- function(df_Individual, zmat, complete_rows, object, min_n,
    df_Training, chance_rate, training_z = NULL, permutations = 999L, seed = 1L) {
  if (length(min_n) != 1L || !is.finite(min_n) || min_n < 1 || min_n != as.integer(min_n)) {
    stop("novel_group_min_n must be a positive integer.")
  }
  nodes <- df_Individual$SOM_Node[complete_rows]
  unusual <- df_Individual$Flag_Unusual[complete_rows] %in% TRUE
  codes <- object$ModelInfo_SOM$som_codes
  adjacent <- .SOMGridAdjacency(object$ModelInfo_SOM$som_model$grid$pts)
  neighbourhoods <- lapply(seq_len(nrow(codes)), function(i) which(adjacent[i, ]))
  training_unusual <- df_Training$Flag_Unusual
  tests <- dplyr::bind_rows(lapply(seq_along(neighbourhoods), function(i) {
    inside <- nodes %in% neighbourhoods[[i]]
    observed <- sum(inside & unusual)
    training_inside <- df_Training$SOM_Node %in% neighbourhoods[[i]] & !is.na(training_unusual)
    training_n <- sum(training_inside)
    training_hits <- sum(training_unusual[training_inside])
    # Shrink toward the chance rate with ten pseudo-people, so a neighbourhood
    # with no unusual training people does not get an expected rate of zero.
    probability <- (training_hits + 10 * chance_rate) / (training_n + 10)
    tibble::tibble(NodeID = i, ProjectedN = sum(inside), TrainingN = training_n,
      TrainingUnusualN = training_hits, TrainingUnusualRate = probability,
      UnusualN = observed, ExpectedN = probability * sum(inside),
      P = if (observed > 0) stats::pbinom(observed - 1L, sum(inside), probability, lower.tail = FALSE) else 1)
  }))
  tests$HolmP <- stats::p.adjust(tests$P, "holm")
  sets <- neighbourhoods[tests$HolmP < .05]
  merged <- list()
  while (length(sets)) {
    group <- sets[[1]]
    sets <- sets[-1]
    repeat {
      overlaps <- vapply(sets, function(x) length(intersect(x, group)) > 0, logical(1))
      if (!any(overlaps)) break
      group <- unique(c(group, unlist(sets[overlaps])))
      sets <- sets[!overlaps]
    }
    merged[[length(merged) + 1L]] <- sort(group)
  }
  node_clusters <- object$ProbFit$node$Cluster
  cluster_profiles <- object$ModelInfo_SOM$ProjectionReference$cluster_profiles
  if (is.null(cluster_profiles)) {
    counts <- tabulate(df_Training$SOM_Node, nrow(codes))
    cluster_profiles <- lapply(split(seq_len(nrow(codes)), node_clusters), function(idx) {
      if (sum(counts[idx]) == 0) return(colMeans(codes[idx, , drop = FALSE]))
      colSums(codes[idx, , drop = FALSE] * counts[idx]) / sum(counts[idx])
    })
  }
  pool_z <- NULL
  if (!is.null(training_z)) {
    training_complete <- which(!is.na(df_Training$SOM_Node))
    if (nrow(training_z) == length(training_complete)) {
      pool_z <- training_z[training_unusual[training_complete] %in% TRUE, , drop = FALSE]
    }
  }
  groups <- dplyr::bind_rows(lapply(seq_along(merged), function(i) {
    selected <- unusual & nodes %in% merged[[i]]
    member_z <- zmat[selected, , drop = FALSE]
    profile <- colMeans(member_z)
    distances <- vapply(cluster_profiles, function(x) sum((x - profile)^2), numeric(1))
    nearest <- which.min(distances)
    coherence <- .SOMGroupCoherence(member_z, pool_z, permutations, seed + i)
    candidate <- sum(selected) >= min_n && !is.na(coherence$P) && coherence$P < .05
    tibble::tibble(Group = i, Nodes = list(merged[[i]]), N = sum(selected),
      Members = list(df_Individual$.row_id[complete_rows][selected]),
      ShareOfUnusual = sum(selected) / sum(unusual),
      Clusters = list(sort(unique(node_clusters[merged[[i]]]))),
      NearestCluster = as.integer(names(cluster_profiles)[nearest]),
      Spread = coherence$Spread, ReferenceSpread = coherence$ReferenceSpread,
      CoherenceP = coherence$P,
      MeanZ = list(profile), TrainingClusterMeanZ = list(cluster_profiles[[nearest]]),
      DifferenceZ = list(profile - cluster_profiles[[nearest]]), Candidate = candidate,
      Label = paste0(if (candidate) "candidate new phenotype" else "unusual-profile cluster on the map",
        " (n=", sum(selected), ")"))
  }))
  if (!length(merged)) groups <- tibble::tibble(Group = integer(), Nodes = list(), N = integer(),
    Members = list(), ShareOfUnusual = numeric(), Clusters = list(), NearestCluster = integer(),
    Spread = numeric(), ReferenceSpread = numeric(), CoherenceP = numeric(),
    MeanZ = list(), TrainingClusterMeanZ = list(), DifferenceZ = list(),
    Candidate = logical(), Label = character())
  candidate <- nrow(groups) > 0 && any(groups$Candidate)
  excess_p <- if (length(unusual)) stats::pbinom(sum(unusual) - 1L, length(unusual), chance_rate,
    lower.tail = FALSE) else NA_real_
  overall <- tibble::tibble(N = length(unusual), UnusualN = sum(unusual),
    ExpectedN = length(unusual) * chance_rate, P = excess_p, Excess = !is.na(excess_p) & excess_p < .05)
  list(groups = groups, neighbourhood_tests = tests, overall_excess = overall,
    verdict = paste(if (candidate) paste(groups$Label[groups$Candidate], collapse = ", ")
      else "no candidate new phenotype",
      if (overall$Excess) "more unusual profiles than expected by chance"
      else "no excess of unusual profiles", sep = "; "))
}
