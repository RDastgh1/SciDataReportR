#' Screen SOM grid sizes from the map alone
#'
#' @description Trains self-organizing maps over a range of square grid sizes
#'   and measures each map before any clustering model is fitted. No value of
#'   `k` and no Mclust covariance model enters the screen, so the chosen map is
#'   fixed before candidate cluster solutions are compared on it.
#'   [CreateClusterModel_SOM_MClust()] runs this screen automatically when no
#'   grid size is supplied.
#'
#'   For every side the screen reports:
#'   \itemize{
#'     \item resolution: variance captured and quantization error (mean distance
#'       from each person to their node);
#'     \item map error: topographic error (share of people whose two closest
#'       nodes are not adjacent);
#'     \item neighbourhood preservation: trustworthiness and continuity on a
#'       sample of people (Venna and Kaski, 2001);
#'     \item sparse nodes: empty, one-person and low-density nodes;
#'     \item stability: Cohen's kappa for "same or adjacent node" between pairs
#'       of maps, both across seeds and across seeds plus jitter (small normal
#'       noise added to every z-score; Hennig, 2007). Every participant is kept
#'       in every refit: subsamples overstate sparse nodes because the final
#'       map is trained on everyone.
#'   }
#'
#'   The grid is then chosen by [SelectSOMGrid()].
#'
#' @param data A data frame.
#' @param variables Numeric clustering variables. Defaults to all numeric
#'   columns.
#' @param ZScoreType,ZScoreObject,ZScoreVars Preprocessing, as in
#'   [CreateClusterModel_SOM_MClust()], so the screen sees the same matrix as
#'   the pipeline.
#' @param sides Optional integer vector of square grid sides to screen. By
#'   default the screen covers two sides below the smaller of the variable floor
#'   and the 5 * sqrt(n) target, up to two sides above the larger of the two,
#'   and never more than one side past the cap. Wide ranges are screened every
#'   other side, always keeping the floor, target and cap.
#' @param cores Workers for per-side refits, using the clustering stability
#'   backend. Default 1. Seeds are fixed per fit; serial and parallel metrics
#'   match. refit_status records failures and table reports seed/jitter success
#'   rates. If all seed fits fail for a side, screening stops with the cause.
#' @param refits Map refits per side for each perturbation (seed, and seed plus
#'   jitter). Default `5`, or `3` when there are more than 5,000 participants.
#' @param jitter_sd Noise standard deviation in z-score units, one value or one
#'   per variable. It sets both the jitter used for stability and the noise
#'   ceiling. When test-retest reliability is known, use
#'   `sqrt(1 - reliability)` per variable. Default `0.25`.
#' @param seed_som SOM seed. The first seed refit of each side is exactly the
#'   map [CreateClusterModel_SOM_MClust()] trains with the same seed.
#' @param som_topo,som_neigh,som_epochs SOM topology, neighbourhood function and
#'   batch epochs, matching the pipeline defaults.
#' @param min_nodes_per_variable,max_nodes_per_variable,max_nodes,max_one_person_prop,max_empty_prop,select
#'   Selection settings passed to [SelectSOMGrid()].
#' @param low_density_n Nodes holding fewer than this many people are counted
#'   as low-density (reported, not used for selection). Default `3`.
#' @param pair_sample Number of person pairs sampled for neighbour kappa.
#' @param trust_sample,trust_k Sample size and neighbourhood size for
#'   trustworthiness and continuity.
#' @param progress Print one line per screened side.
#'
#' @return A `SOMGridScreen` object with `table` (one row per side),
#'   `refit_status` (Side, Refit, Type, Status, Error for seed and jitter fits),
#'   `table$SuccessRate_seed` and `table$SuccessRate_jitter` (shares fitted),
#'   `refits` (one row per seed refit), `pairs` (one row per pair of refits),
#'   `selection` (chosen side, floor, ceilings, reason, warnings) and
#'   `settings`. Use `plot()` for the screening figure and
#'   [SOMGridDecisionTable()] for a formatted table.
#'
#' @references
#' Vesanto J, Himberg J, Alhoniemi E, Parhankangas J (2000). SOM Toolbox for
#' Matlab 5. Helsinki University of Technology, Report A57.
#'
#' Kiviluoto K (1996). Topology preservation in self-organizing maps.
#' Proceedings of the IEEE International Conference on Neural Networks, 294-299.
#'
#' Venna J, Kaski S (2001). Neighborhood preservation in nonlinear projection
#' methods: an experimental study. ICANN 2001, 485-491.
#'
#' Hennig C (2007). Cluster-wise assessment of cluster stability.
#' Computational Statistics & Data Analysis 52, 258-271.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- as.data.frame(matrix(rnorm(300 * 4), ncol = 4))
#' screen <- ScreenSOMGrid(df, names(df), refits = 2)
#' screen
#' plot(screen)
#' SOMGridDecisionTable(screen)
#' # Change a threshold without refitting any maps
#' SelectSOMGrid(screen, max_one_person_prop = 0.05)
#' }
#' @export
ScreenSOMGrid <- function(data, variables = NULL,
    ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj", "PreZScored"),
    ZScoreObject = NULL, ZScoreVars = NULL,
    sides = NULL, refits = NULL, jitter_sd = 0.25, seed_som = 934521L,
    som_topo = "hexagonal", som_neigh = "gaussian", som_epochs = 100L,
    min_nodes_per_variable = 10, max_nodes_per_variable = 30, max_nodes = NULL,
    max_one_person_prop = 0.10, max_empty_prop = 0.05, low_density_n = 3L,
    pair_sample = 20000L, trust_sample = 1000L, trust_k = 10L,
    select = c("target", "largest_feasible"), progress = FALSE, cores = 1L) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  ZScoreType <- match.arg(ZScoreType)
  select <- match.arg(select)
  .ValidateClusterStabilityCores(cores)
  if (is.null(variables)) variables <- names(data)[vapply(data, is.numeric, logical(1))]
  if (!is.character(variables) || !length(variables) || length(setdiff(variables, names(data)))) {
    stop("variables must name columns in data.")
  }
  zmat <- if (ZScoreType == "PreZScored") {
    as.matrix(data[, if (is.null(ZScoreVars)) variables else ZScoreVars, drop = FALSE])
  } else if (ZScoreType == "ZScoreObj") {
    if (is.null(ZScoreObject) || !"ZScoreObj" %in% class(ZScoreObject)) {
      stop("ZScoreType = 'ZScoreObj' requires a valid ZScoreObj.")
    }
    as.matrix(ProjectZScore(data = data, variables = variables, parameters = ZScoreObject,
      ParameterInputType = "ZScoreObj", names_prefix = "Z_")$ZScores)
  } else {
    as.matrix(CreateZScoreObject(data = data, variables = variables, names_prefix = "Z_",
      center = ZScoreType %in% c("Center and Scale", "Center Only"),
      scale = ZScoreType %in% c("Center and Scale", "Scale Only"))$ZScores)
  }
  zmat <- zmat[stats::complete.cases(zmat), , drop = FALSE]
  .ScreenSOMGridMatrix(zmat, sides = sides, refits = refits, jitter_sd = jitter_sd,
    seed_som = seed_som, som_topo = som_topo, som_neigh = som_neigh, som_epochs = som_epochs,
    min_nodes_per_variable = min_nodes_per_variable, max_nodes_per_variable = max_nodes_per_variable,
    max_nodes = max_nodes, max_one_person_prop = max_one_person_prop, max_empty_prop = max_empty_prop,
    low_density_n = low_density_n, pair_sample = pair_sample, trust_sample = trust_sample,
    trust_k = trust_k, select = select, progress = progress, cores = cores)
}

# Screen on an already preprocessed, complete matrix (used by the pipeline).
.ScreenSOMGridMatrix <- function(zmat, sides = NULL, refits = NULL, jitter_sd = 0.25,
    seed_som = 934521L, som_topo = "hexagonal", som_neigh = "gaussian", som_epochs = 100L,
    min_nodes_per_variable = 10, max_nodes_per_variable = 30, max_nodes = NULL,
    max_one_person_prop = 0.10, max_empty_prop = 0.05, low_density_n = 3L,
    pair_sample = 20000L, trust_sample = 1000L, trust_k = 10L,
    select = c("target", "largest_feasible"), progress = FALSE, cores = 1L) {
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required for ScreenSOMGrid().")
  }
  select <- match.arg(select)
  .ValidateClusterStabilityCores(cores)
  zmat <- as.matrix(zmat)
  n <- nrow(zmat)
  p <- ncol(zmat)
  if (n < 10L) stop("At least 10 complete participants are required to screen SOM grids.")
  if (!is.numeric(jitter_sd) || !length(jitter_sd) || any(!is.finite(jitter_sd)) || any(jitter_sd < 0) ||
      !length(jitter_sd) %in% c(1L, p)) {
    stop("jitter_sd must be one non-negative number or one per variable.")
  }
  jitter_sd <- rep_len(jitter_sd, p)
  if (is.null(refits)) refits <- if (n > 5000L) 3L else 5L
  if (!is.numeric(refits) || length(refits) != 1L || refits < 2L) stop("refits must be at least 2.")
  refits <- as.integer(refits)

  limits <- .SOMGridLimits(n, p, min_nodes_per_variable, max_nodes_per_variable, max_nodes)
  if (is.null(sides)) {
    lo <- max(2L, min(limits$floor, limits$target) - 2L)
    hi <- max(min(max(limits$floor, limits$target) + 2L, limits$cap + 1L), limits$floor + 1L)
    sides <- seq.int(lo, hi)
    # Wide ranges: every other side, always keeping the floor, target and cap
    if (length(sides) > 9L) {
      sides <- sort(unique(c(sides[seq(1L, length(sides), by = 2L)], limits$floor, limits$target,
        limits$cap, hi)))
      sides <- sides[sides >= lo & sides <= hi]
    }
  }
  if (!is.numeric(sides) || !length(sides) || any(sides < 2) || any(sides != as.integer(sides))) {
    stop("sides must be whole numbers of at least 2.")
  }
  sides <- sort(unique(as.integer(sides)))

  set.seed(seed_som + 1L)
  pairs <- cbind(sample.int(n, pair_sample, TRUE), sample.int(n, pair_sample, TRUE))
  pairs <- pairs[pairs[, 1] != pairs[, 2], , drop = FALSE]
  idx_Trust <- sort(sample.int(n, min(as.integer(trust_sample), n)))

  per_Refit <- vector("list", length(sides))
  per_Pair <- vector("list", length(sides))
  per_Trust <- vector("list", length(sides))
  per_Status <- vector("list", length(sides))
  fit_map <- .FitSOMGridScreen
  map_quality <- .SOMGridMapQuality
  for (i in seq_along(sides)) {
    side <- sides[[i]]
    t0 <- Sys.time()
    seed_Maps <- .ClusterParallelMap(seq_len(refits), function(r) tryCatch({
      m <- fit_map(zmat, side, seed_som + r - 1L, som_topo, som_neigh, som_epochs)
      list(quality = map_quality(m, zmat, low_density_n), pts = m$grid$pts)
    }, error = function(e) e), cores = cores)
    valid_seed <- !vapply(seed_Maps, inherits, logical(1), "error")
    if (!any(valid_seed)) stop("All SOM grid-screen refits failed for side ", side,
      ": ", conditionMessage(seed_Maps[[1]]))
    first_valid <- which(valid_seed)[1]
    pts <- seed_Maps[[first_valid]]$pts
    gd <- as.matrix(stats::dist(pts))
    set.seed(seed_som + 2L)
    per_Trust[[i]] <- dplyr::tibble(Side = side,
      !!!as.list(.SOMGridTrustContinuity(zmat, seed_Maps[[first_valid]]$quality$bmu, pts, idx_Trust, trust_k)))
    jitter_Bmu <- .ClusterParallelMap(seq_len(refits), function(r) tryCatch({
      set.seed(seed_som + 5000L + r)
      zj <- zmat + sweep(matrix(stats::rnorm(n * p), n, p), 2, jitter_sd, "*")
      m <- fit_map(zj, side, seed_som + 1000L + r, som_topo, som_neigh, som_epochs)
      map_quality(m, zmat, low_density_n)$bmu   # people placed by their real scores
    }, error = function(e) e), cores = cores)
    valid_jitter <- !vapply(jitter_Bmu, inherits, logical(1), "error")
    per_Status[[i]] <- dplyr::bind_rows(lapply(c("seed", "jitter"), function(type) {
      results <- if (type == "seed") seed_Maps else jitter_Bmu
      valid <- if (type == "seed") valid_seed else valid_jitter
      dplyr::tibble(Side = side, Refit = seq_len(refits), Type = type,
        Status = ifelse(valid, "success", "failed"),
        Error = vapply(results, function(x) if (inherits(x, "error")) conditionMessage(x) else NA_character_, character(1)))
    }))
    ij <- utils::combn(refits, 2)
    per_Refit[[i]] <- dplyr::bind_cols(dplyr::tibble(Side = side, Refit = seq_len(refits)),
      dplyr::as_tibble(do.call(rbind, lapply(seed_Maps, function(x) if (inherits(x, "error"))
        seed_Maps[[first_valid]]$quality$metrics * NA_real_ else x$quality$metrics))))
    # Seed-only refits often converge to the same arrangement (batch training starts with a map-wide
    # neighbourhood), so kappa across seeds is lumpy; jitter kappa is the smoother stability measure.
    per_Pair[[i]] <- dplyr::tibble(Side = side, PairA = ij[1, ], PairB = ij[2, ],
      KappaSeed = apply(ij, 2, function(x) if (!all(valid_seed[x])) NA_real_ else .SOMGridNeighbourKappa(seed_Maps[[x[1]]]$quality$bmu,
        seed_Maps[[x[2]]]$quality$bmu, gd, pairs)),
      KappaJitter = apply(ij, 2, function(x) if (!all(valid_jitter[x])) NA_real_ else .SOMGridNeighbourKappa(jitter_Bmu[[x[1]]],
        jitter_Bmu[[x[2]]], gd, pairs)))
    if (isTRUE(progress)) {
      message(sprintf("SOM grid screen: %d x %d done (%.1f s)", side, side,
        as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  df_Refits <- dplyr::bind_rows(per_Refit)
  df_Pairs <- dplyr::bind_rows(per_Pair)
  df_Status <- dplyr::bind_rows(per_Status)
  success_rates <- df_Status %>% dplyr::group_by(.data$Side, .data$Type) %>%
    dplyr::summarise(SuccessRate = mean(Status == "success"), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = "Type", values_from = "SuccessRate", names_prefix = "SuccessRate_")

  df_Table <- df_Refits %>%
    dplyr::group_by(.data$Side) %>%
    dplyr::summarise(dplyr::across(-dplyr::all_of("Refit"), list(Mean = function(x) mean(x, na.rm = TRUE), SD = function(x) stats::sd(x, na.rm = TRUE))), .groups = "drop") %>%
    dplyr::left_join(df_Pairs %>% dplyr::group_by(.data$Side) %>%
      dplyr::summarise(KappaSeed_SD = stats::sd(.data$KappaSeed, na.rm = TRUE), KappaSeed = mean(.data$KappaSeed, na.rm = TRUE),
        KappaJitter_SD = stats::sd(.data$KappaJitter, na.rm = TRUE), KappaJitter = mean(.data$KappaJitter, na.rm = TRUE),
        .groups = "drop"), by = "Side") %>%
    dplyr::left_join(dplyr::bind_rows(per_Trust), by = "Side") %>%
    dplyr::left_join(success_rates, by = "Side") %>%
    dplyr::mutate(Nodes = .data$Side^2, PeoplePerNode = n / .data$Nodes, .after = "Side") %>%
    dplyr::rename_with(~ sub("_Mean$", "", .x)) %>%
    dplyr::mutate(
      OnePersonSmooth = .SOMGridSmoothSparse(df_Refits, "OnePersonNodes", n, .data$Side),
      EmptySmooth = .SOMGridSmoothSparse(df_Refits, "EmptyNodes", n, .data$Side))

  out <- list(table = df_Table, refits = df_Refits, pairs = df_Pairs,
    refit_status = df_Status,
    settings = list(n = n, p = p, variables = colnames(zmat), refits = refits, cores = cores, jitter_sd = jitter_sd,
      seed_som = seed_som, som_topo = som_topo, som_neigh = som_neigh, som_epochs = som_epochs,
      low_density_n = low_density_n, pair_sample = nrow(pairs), trust_sample = length(idx_Trust),
      trust_k = trust_k))
  class(out) <- "SOMGridScreen"
  SelectSOMGrid(out, min_nodes_per_variable = min_nodes_per_variable,
    max_nodes_per_variable = max_nodes_per_variable, max_nodes = max_nodes,
    max_one_person_prop = max_one_person_prop, max_empty_prop = max_empty_prop, select = select)
}

#' Choose a SOM grid size from a grid screen
#'
#' @description Applies the map-only selection rule to a [ScreenSOMGrid()]
#'   result. Thresholds can be changed and the choice redone without refitting
#'   any maps.
#'
#'   \enumerate{
#'     \item Floor: `min_nodes_per_variable` nodes per clustering variable
#'       (default 10), so the mixture step has enough nodes to estimate spread in
#'       every dimension.
#'     \item Ceiling: the smallest of
#'       \itemize{
#'         \item the sparsity ceiling: the largest side before the smoothed
#'           one-person share exceeds `max_one_person_prop` or the smoothed empty
#'           share exceeds `max_empty_prop`;
#'         \item the noise ceiling: the first side whose quantization error is
#'           below the expected distance between a person and a retest,
#'           `sqrt(sum(noise_sd^2))`; larger maps only separate people by
#'           differences smaller than measurement noise;
#'         \item the cap: `max_nodes_per_variable` nodes per variable (default
#'           30), or `max_nodes` if smaller.
#'       }
#'     \item Choice: the 5 * sqrt(n) target (Vesanto et al., 2000), moved into
#'       the allowed range. If the floor is above the ceiling, the target is kept
#'       between the two and a warning is given.
#'   }
#'
#'   The chosen map is also checked: a warning is given if its topographic error
#'   is high compared with the other sides screened (a folded map), or if its
#'   trustworthiness is below 0.80.
#'
#' @param x A `SOMGridScreen` object.
#' @param min_nodes_per_variable Floor, in nodes per clustering variable.
#' @param max_nodes_per_variable Cap, in nodes per clustering variable.
#' @param max_nodes Optional absolute cap on the number of nodes.
#' @param max_one_person_prop Largest allowed smoothed share of one-person nodes.
#' @param max_empty_prop Largest allowed smoothed share of empty nodes.
#' @param noise_sd Measurement noise SD in z-score units (one value or one per
#'   variable) for the noise ceiling. Defaults to the jitter used in the screen.
#' @param select `"target"` (default) moves the 5 * sqrt(n) target into the
#'   allowed range; `"largest_feasible"` takes the top of the allowed range.
#'
#' @return The screen with updated `table` pass/fail columns and `selection`.
#' @export
SelectSOMGrid <- function(x, min_nodes_per_variable = 10, max_nodes_per_variable = 30,
    max_nodes = NULL, max_one_person_prop = 0.10, max_empty_prop = 0.05,
    noise_sd = x$settings$jitter_sd, select = c("target", "largest_feasible")) {
  if (!inherits(x, "SOMGridScreen")) stop("x must be a SOMGridScreen object from ScreenSOMGrid().")
  select <- match.arg(select)
  for (nm in c("min_nodes_per_variable", "max_nodes_per_variable", "max_one_person_prop", "max_empty_prop")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) stop(nm, " must be one positive number.")
  }
  if (max_nodes_per_variable < min_nodes_per_variable) {
    stop("max_nodes_per_variable must be at least min_nodes_per_variable.")
  }
  n <- x$settings$n
  p <- x$settings$p
  limits <- .SOMGridLimits(n, p, min_nodes_per_variable, max_nodes_per_variable, max_nodes)
  noise_Distance <- sqrt(sum(rep_len(noise_sd, p)^2))
  df_Table <- x$table %>%
    dplyr::mutate(
      PassFloor = .data$Side >= limits$floor,
      PassSparsity = .data$OnePersonSmooth <= max_one_person_prop & .data$EmptySmooth <= max_empty_prop,
      FinerThanNoise = .data$QuantizationError < noise_Distance,
      PassCap = .data$Side <= limits$cap)

  fail_Sparse <- df_Table$Side[!df_Table$PassSparsity]
  ceil_Sparse <- if (length(fail_Sparse)) min(fail_Sparse) - 1L else Inf
  ceil_Noise <- if (any(df_Table$FinerThanNoise)) min(df_Table$Side[df_Table$FinerThanNoise]) else Inf
  ceilings <- c(sparsity = ceil_Sparse, noise = ceil_Noise, cap = limits$cap,
    screened_range = max(df_Table$Side))
  side_Ceiling <- as.integer(min(ceilings))
  binding <- names(ceilings)[which.min(ceilings)]
  lab <- c(sparsity = "sparsity", noise = "noise",
    cap = sprintf("%g-nodes-per-variable cap", max_nodes_per_variable),
    screened_range = "screened-range")[[binding]]
  if (binding == "cap" && !is.null(max_nodes) && floor(sqrt(max_nodes)) < floor(sqrt(max_nodes_per_variable * p))) {
    lab <- sprintf("%g-node cap", max_nodes)
  }

  warn <- character()
  if (limits$floor > side_Ceiling) {
    lo <- side_Ceiling
    hi <- limits$floor
    warn <- c(warn, sprintf(paste0("The variable floor (%d x %d, %g nodes per variable) is above the %s ceiling ",
      "(%d x %d): the sample is small for %d variables. The 5 * sqrt(n) target is kept between the two."),
      limits$floor, limits$floor, min_nodes_per_variable, lab, side_Ceiling, side_Ceiling, p))
  } else {
    lo <- limits$floor
    hi <- side_Ceiling
  }
  side_Chosen <- as.integer(if (select == "target") min(max(limits$target, lo), hi) else hi)
  reason <- if (select == "largest_feasible") {
    sprintf("Largest side in the allowed range (%s ceiling).", lab)
  } else if (side_Chosen == limits$target) {
    sprintf("The 5 * sqrt(n) target (%d x %d) lies inside the allowed range %d x %d to %d x %d.",
      limits$target, limits$target, lo, lo, hi, hi)
  } else if (side_Chosen < limits$target) {
    sprintf("The 5 * sqrt(n) target (%d x %d) was moved down to the %s ceiling (%d x %d).",
      limits$target, limits$target, lab, side_Chosen, side_Chosen)
  } else {
    sprintf("The 5 * sqrt(n) target (%d x %d) was moved up to the variable floor (%d x %d).",
      limits$target, limits$target, side_Chosen, side_Chosen)
  }
  if (binding == "screened_range" && side_Chosen == side_Ceiling) {
    warn <- c(warn, "The ceiling is the largest side screened; screen larger sides to find the real ceiling.")
  }
  if (!side_Chosen %in% df_Table$Side) {
    warn <- c(warn, sprintf("The chosen side (%d x %d) was not screened; its map checks are unavailable.",
      side_Chosen, side_Chosen))
  } else {
    row <- df_Table[df_Table$Side == side_Chosen, , drop = FALSE]
    te <- df_Table$TopographicError
    if (nrow(df_Table) >= 3L && row$TopographicError > 0.05 &&
        row$TopographicError > stats::median(te) + 2 * stats::mad(te, constant = 1) + 0.02) {
      warn <- c(warn, sprintf("Topographic error on the chosen map (%.2f) is high compared with the other sides screened: the map may be folded.",
        row$TopographicError))
    }
    if (is.finite(row$Trustworthiness) && row$Trustworthiness < 0.80) {
      warn <- c(warn, sprintf("Trustworthiness on the chosen map is %.2f (< 0.80): neighbours on the map are often not neighbours in the data.",
        row$Trustworthiness))
    }
  }
  x$table <- df_Table %>% dplyr::mutate(Chosen = .data$Side == side_Chosen)
  x$selection <- list(side = side_Chosen, floor = limits$floor, ceiling = side_Ceiling,
    ceiling_binding = binding, ceiling_label = lab, ceilings = ceilings, target = limits$target,
    allowed = c(lo, hi), rule = select, reason = reason, warnings = warn)
  x$settings <- utils::modifyList(x$settings, list(noise_distance = noise_Distance,
    min_nodes_per_variable = min_nodes_per_variable, max_nodes_per_variable = max_nodes_per_variable,
    max_nodes = max_nodes, max_one_person_prop = max_one_person_prop, max_empty_prop = max_empty_prop))
  for (w in warn) warning(w, call. = FALSE)
  x
}

#' @export
print.SOMGridScreen <- function(x, ...) {
  s <- x$selection
  cat(sprintf("SOM grid screen: n = %d, %d variables, %d refits per side (seed, and seed + jitter SD %s)\n",
    x$settings$n, x$settings$p, x$settings$refits,
    paste(unique(round(x$settings$jitter_sd, 2)), collapse = "/")))
  cat(sprintf("Allowed %d x %d to %d x %d (floor: %g nodes per variable; ceiling: %s). Chosen: %d x %d.\n%s\n",
    s$allowed[[1]], s$allowed[[1]], s$allowed[[2]], s$allowed[[2]], x$settings$min_nodes_per_variable,
    s$ceiling_label, s$side, s$side, s$reason))
  if (length(s$warnings)) cat(paste0("! ", s$warnings, collapse = "\n"), "\n")
  cat("See plot() and SOMGridDecisionTable(); change thresholds with SelectSOMGrid().\n")
  invisible(x)
}

#' Formatted decision table for a SOM grid screen
#'
#' @param x A `SOMGridScreen` object.
#' @return A data frame with one row per screened side: map measures and
#'   pass/fail on each selection criterion.
#' @export
SOMGridDecisionTable <- function(x) {
  if (!inherits(x, "SOMGridScreen")) stop("x must be a SOMGridScreen object from ScreenSOMGrid().")
  pct <- function(v) sprintf("%.0f%%", 100 * v)
  tick <- function(v) ifelse(v, "\u2713", "\u2717")
  t <- x$table
  data.frame(
    Grid = paste0(t$Side, " x ", t$Side),
    `People per node` = round(t$PeoplePerNode, 1),
    `Variance captured` = round(t$VarianceCaptured, 3),
    `Quantization error` = round(t$QuantizationError, 2),
    `Topographic error` = round(t$TopographicError, 2),
    Trustworthiness = round(t$Trustworthiness, 3),
    Continuity = round(t$Continuity, 3),
    Empty = pct(t$EmptySmooth),
    `One-person` = pct(t$OnePersonSmooth),
    `Low-density` = pct(t$LowDensityNodes),
    `Kappa jitter` = round(t$KappaJitter, 2),
    `Kappa seeds` = round(t$KappaSeed, 2),
    Floor = tick(t$PassFloor),
    Sparsity = tick(t$PassSparsity),
    `Above noise` = tick(!t$FinerThanNoise),
    `Under cap` = tick(t$PassCap),
    Chosen = ifelse(t$Chosen, "\u25c0", ""),
    check.names = FALSE)
}

#' @export
plot.SOMGridScreen <- function(x, ...) {
  pal <- unname(SciDataPalette()[c("Navy", "Orange", "Emerald", "Cranberry", "Plum")])
  s <- x$selection
  lev <- c("Resolution: variance captured", "Map error (dashed = noise level)",
    "Neighbourhood preservation", "Sparse nodes (dashed = caps)", "Stability (neighbour kappa)")
  long <- function(df, cols) {
    dplyr::bind_rows(lapply(names(cols), function(nm)
      dplyr::tibble(Side = df$Side, Measure = nm, Value = df[[cols[[nm]]]])))
  }
  df <- dplyr::bind_rows(
    long(x$refits, c(`Variance captured` = "VarianceCaptured", `Quantization error` = "QuantizationError",
      `Topographic error` = "TopographicError", Empty = "EmptyNodes", `One person` = "OnePersonNodes",
      `Low density` = "LowDensityNodes")),
    long(x$table, c(Trustworthiness = "Trustworthiness", Continuity = "Continuity")),
    long(x$pairs, c(`Across seeds` = "KappaSeed", `With jitter` = "KappaJitter"))) %>%
    dplyr::mutate(Panel = factor(dplyr::case_when(
      .data$Measure == "Variance captured" ~ lev[1],
      .data$Measure %in% c("Quantization error", "Topographic error") ~ lev[2],
      .data$Measure %in% c("Trustworthiness", "Continuity") ~ lev[3],
      .data$Measure %in% c("Empty", "One person", "Low density") ~ lev[4],
      TRUE ~ lev[5]), levels = lev))
  df_Mean <- df %>%
    dplyr::group_by(.data$Panel, .data$Measure, .data$Side) %>%
    dplyr::summarise(Mean = mean(.data$Value, na.rm = TRUE),
      Lo = stats::quantile(.data$Value, .25, na.rm = TRUE),
      Hi = stats::quantile(.data$Value, .75, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$Panel) %>%
    dplyr::mutate(Colour = pal[as.integer(factor(.data$Measure))]) %>%
    dplyr::ungroup()
  df_Lab <- df_Mean %>% dplyr::group_by(.data$Panel, .data$Measure) %>%
    dplyr::filter(.data$Side == max(.data$Side)) %>% dplyr::ungroup()
  df_Ref <- dplyr::tibble(Panel = factor(lev[c(2, 4, 4)], levels = lev),
    y = c(x$settings$noise_distance, x$settings$max_one_person_prop, x$settings$max_empty_prop))
  rng <- range(df$Side)
  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = df_Lab, ggplot2::aes(label = .data$Measure, colour = I(.data$Colour)),
      hjust = 0, direction = "y", nudge_x = .04 * diff(rng) + .2, xlim = c(rng[2] + .2, NA),
      size = 2.7, segment.colour = NA, min.segment.length = Inf, box.padding = .1)
  } else {
    ggplot2::geom_text(data = df_Lab, ggplot2::aes(label = .data$Measure, colour = I(.data$Colour)),
      hjust = 0, nudge_x = .3, size = 2.7)
  }
  ggplot2::ggplot(df_Mean, ggplot2::aes(.data$Side, .data$Mean, group = .data$Measure)) +
    ggplot2::annotate("rect", xmin = s$allowed[[1]] - .5, xmax = s$allowed[[2]] + .5,
      ymin = -Inf, ymax = Inf, fill = "grey93") +
    ggplot2::geom_vline(xintercept = s$target, colour = "grey60", linetype = 3) +
    ggplot2::geom_vline(xintercept = s$side, colour = "grey25", linewidth = .6) +
    ggplot2::geom_hline(data = df_Ref, ggplot2::aes(yintercept = .data$y), linetype = 2, colour = "grey55") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Lo, ymax = .data$Hi, fill = I(.data$Colour)), alpha = .15) +
    ggplot2::geom_line(ggplot2::aes(colour = I(.data$Colour)), linewidth = .8) +
    ggplot2::geom_point(ggplot2::aes(colour = I(.data$Colour)), size = 1.2) +
    label_layer +
    ggplot2::facet_wrap(~Panel, scales = "free_y", ncol = 3) +
    ggplot2::scale_x_continuous(
      breaks = function(l) { b <- unique(round(pretty(rng, n = 5))); b[b >= rng[1] & b <= rng[2]] },
      labels = function(b) paste0(b, " x ", b),
      expand = ggplot2::expansion(add = c(.3, .3 + diff(rng) * .55))) +
    ggplot2::labs(x = sprintf(
      "Grid side. Shaded = allowed %d x %d to %d x %d (%s ceiling); dotted = 5 * sqrt(n) target %d x %d; solid = chosen %d x %d",
      s$allowed[[1]], s$allowed[[1]], s$allowed[[2]], s$allowed[[2]], s$ceiling_label,
      s$target, s$target, s$side, s$side), y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7), axis.title.x = ggplot2::element_text(size = 8))
}

# ---- Internals ---------------------------------------------------------------------------------

.SOMGridLimits <- function(n, p, min_nodes_per_variable, max_nodes_per_variable, max_nodes) {
  cap <- floor(sqrt(max_nodes_per_variable * p))
  if (!is.null(max_nodes)) {
    if (!is.numeric(max_nodes) || length(max_nodes) != 1L || !is.finite(max_nodes) || max_nodes < 4) {
      stop("max_nodes must be NULL or one number of at least 4.")
    }
    cap <- min(cap, floor(sqrt(max_nodes)))
  }
  list(target = as.integer(max(2L, round(sqrt(5 * sqrt(n))))),
    floor = as.integer(max(2L, ceiling(sqrt(min_nodes_per_variable * p)))),
    cap = as.integer(max(2L, cap)))
}

# Same training as FitSOMForGrid() in CreateClusterModel_SOM_MClust(): seeded sample of people as the
# initial codebook, batch supersom.
.FitSOMGridScreen <- function(zmat, side, seed, som_topo, som_neigh, som_epochs) {
  grid <- kohonen::somgrid(xdim = side, ydim = side, topo = som_topo, neighbourhood.fct = som_neigh)
  set.seed(seed)
  init <- zmat[sample.int(nrow(zmat), nrow(grid$pts), replace = nrow(zmat) < nrow(grid$pts)), , drop = FALSE]
  kohonen::supersom(list(zmat), grid = grid, init = list(init), mode = "batch",
    rlen = som_epochs, keep.data = FALSE, cores = 1L)
}

.SOMGridMapQuality <- function(m, zmat, low_density_n) {
  C <- as.matrix(if (is.list(m$codes)) m$codes[[1]] else m$codes)
  d2 <- pmax(outer(rowSums(zmat^2), rowSums(C^2), "+") - 2 * tcrossprod(zmat, C), 0)
  ord <- t(apply(d2, 1, function(v) order(v)[1:2]))
  bmu <- ord[, 1]
  gd <- as.matrix(stats::dist(m$grid$pts))
  hits <- tabulate(bmu, nrow(C))
  d_bmu <- d2[cbind(seq_len(nrow(zmat)), bmu)]
  list(bmu = bmu, metrics = c(
    QuantizationError = mean(sqrt(d_bmu)),
    VarianceCaptured = 1 - mean(d_bmu) / sum(apply(zmat, 2, stats::var)),
    TopographicError = mean(gd[cbind(bmu, ord[, 2])] > 1.01),
    EmptyNodes = mean(hits == 0),
    OnePersonNodes = mean(hits == 1),
    LowDensityNodes = mean(hits > 0 & hits < low_density_n),
    PeopleOnLowDensity = mean(hits[bmu] < low_density_n)))
}

# Cohen's kappa on "same or adjacent node" for sampled pairs of people; kappa removes the chance agreement
# that grows on small maps.
.SOMGridNeighbourKappa <- function(bmu_a, bmu_b, gd, pairs) {
  A <- gd[cbind(bmu_a[pairs[, 1]], bmu_a[pairs[, 2]])] < 1.01
  B <- gd[cbind(bmu_b[pairs[, 1]], bmu_b[pairs[, 2]])] < 1.01
  pe <- mean(A) * mean(B) + (1 - mean(A)) * (1 - mean(B))
  if (pe >= 1) return(NA_real_)
  (mean(A == B) - pe) / (1 - pe)
}

# Trustworthiness and continuity (Venna & Kaski 2001); map position = grid coordinates of the node,
# ties within a node broken by tiny noise.
.SOMGridTrustContinuity <- function(zmat, bmu, pts, idx, k) {
  m <- length(idx)
  k <- min(as.integer(k), floor((m - 1) / 3))
  if (k < 1L) return(c(Trustworthiness = NA_real_, Continuity = NA_real_))
  Dd <- as.matrix(stats::dist(zmat[idx, , drop = FALSE]))
  Dm <- as.matrix(stats::dist(pts[bmu[idx], , drop = FALSE] + matrix(stats::rnorm(2 * m, 0, 1e-3), m)))
  diag(Dd) <- Inf
  diag(Dm) <- Inf
  Rd <- t(apply(Dd, 1, rank))
  Rm <- t(apply(Dm, 1, rank))
  Nd <- Rd <= k
  Nm <- Rm <= k
  scale <- 2 / (m * k * (2 * m - 3 * k - 1))
  c(Trustworthiness = 1 - scale * sum((Rd - k)[Nm & !Nd]),
    Continuity = 1 - scale * sum((Rm - k)[Nd & !Nm]))
}

# Sparse-node shares are noisy per side, so a hard cap would flip between neighbouring sides. One smooth
# curve across all sides and refits: logit(share) ~ log(people per node), weighted by node count.
.SOMGridSmoothSparse <- function(df_Refits, measure, n, sides_out) {
  d <- data.frame(y = df_Refits[[measure]], w = df_Refits$Side^2, x = log(n / df_Refits$Side^2))
  d <- d[is.finite(d$y), , drop = FALSE]
  if (all(d$y == 0)) return(rep(0, length(sides_out)))
  if (length(unique(d$x)) < 2L) return(rep(mean(d$y), length(sides_out)))
  fit <- suppressWarnings(stats::glm(y ~ x, family = stats::quasibinomial(), weights = d$w, data = d))
  as.numeric(stats::predict(fit, newdata = data.frame(x = log(n / sides_out^2)), type = "response"))
}

# Warning for candidate k with fewer than candidate_min_nodes_per_cluster nodes each on a fixed map
.WarnSOMNodesPerCluster <- function(n_nodes, k_values, candidate_min_nodes_per_cluster) {
  if (is.null(candidate_min_nodes_per_cluster)) return(invisible(integer()))
  low <- sort(unique(k_values[n_nodes / k_values < candidate_min_nodes_per_cluster]))
  if (length(low)) {
    warning(sprintf(paste0("The %d-node map gives fewer than %g SOM nodes per cluster for k = %s. ",
      "Interpret these candidates cautiously, or use a larger grid."),
      n_nodes, candidate_min_nodes_per_cluster, paste(low, collapse = ", ")), call. = FALSE)
  }
  invisible(low)
}
