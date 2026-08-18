# Canonical significance-star formatting for the whole package.
#
# Thresholds are INCLUSIVE upper bounds: p <= 0.05 is "*". This matches
# rstatix::add_significance(), which PlotPhiHeatmap(),
# PlotPointCorrelationsHeatmap() and PlotAnovaRelationshipsMatrix() already use,
# so a p-value of exactly 0.05 or 0.001 now earns the same number of stars
# everywhere. Before 20.24.0 the package carried six separate implementations
# and they disagreed at the boundary.
#
# `tiers` is a named numeric vector of upper bounds; names are the labels. The
# default four-tier scale is rstatix's. Pass a three-tier vector for plots that
# do not distinguish p <= 0.0001.
ScidrPValueStars <- function(p,
    tiers = c("****" = 0.0001, "***" = 0.001, "**" = 0.01, "*" = 0.05),
    ns_label = "ns",
    na_label = NA_character_) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep(ns_label, length(p))

  # Assign least-stringent first so the most stringent tier wins the overwrite.
  for (i in order(tiers, decreasing = TRUE)) {
    out[!is.na(p) & p <= tiers[[i]]] <- names(tiers)[i]
  }

  out[is.na(p)] <- na_label
  out
}

# Three-tier scale used by the plotting functions.
ScidrStarTiersThree <- c("***" = 0.001, "**" = 0.01, "*" = 0.05)

# Single caption string describing the scale, so plots cannot drift from the
# thresholds the code actually applies. ASCII only (CRAN).
ScidrStarCaptionText <- function() {
  "* = p <= 0.05, ** = p <= 0.01, *** = p <= 0.001"
}
