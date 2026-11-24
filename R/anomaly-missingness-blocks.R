# Variables collected on the same form go blank together, so the missing/observed
# matrix carries the form structure without anyone naming it. Clustering columns
# by how similarly they go blank recovers those groups, and the group's member
# list is the explanation a reader needs: "you have 7 of these 8, and almost
# everyone who has any of them has all 8."
#
# Discovery needs co-missing events. A column whose blanks coincide with nothing
# cannot be placed, and an unplaced column yields no finding: there is no
# comparison group to judge it against, and "this column is mostly full" is a
# completeness fact, not an anomaly.
ScidrAnomalyMissingnessBlocks <- function(missing, args) {
  if (!is.matrix(missing) || !nrow(missing) || !ncol(missing)) return(list())
  counts <- colSums(missing)
  # A column blank in every record carries no contrast; one blank in too few
  # records cannot establish membership of anything.
  usable <- counts >= args$min_block_events & counts < nrow(missing)
  vars <- colnames(missing)[usable]
  if (length(vars) < max(2L, args$min_block_size)) return(list())
  m_Missing <- missing[, vars, drop = FALSE] * 1
  intersection <- crossprod(m_Missing)
  sizes <- counts[vars]
  union <- outer(sizes, sizes, "+") - intersection
  jaccard <- intersection / pmax(union, 1)
  # Complete linkage, so every pair inside a block meets the similarity bar.
  # Single linkage would chain two unrelated forms through one shared skipper.
  tree <- stats::hclust(stats::as.dist(1 - jaccard), method = "complete")
  blocks <- split(vars, stats::cutree(tree, h = 1 - args$block_similarity))
  blocks <- blocks[lengths(blocks) >= args$min_block_size]
  unname(lapply(blocks, sort))
}

# Coverage of one block per record: how many of its members the record carries.
ScidrAnomalyMissingnessCoverage <- function(missing, block) {
  observed <- rowSums(!missing[, block, drop = FALSE])
  list(Observed = observed, Size = length(block),
    Started = observed >= 1L, Complete = observed == length(block))
}

ScidrAnomalyMissingnessBlockText <- function(block, labels, limit = 6L) {
  names <- vapply(block, function(variable) {
    label <- labels[variable]
    if (!length(label) || is.na(label) || !nzchar(label)) variable else unname(label)
  }, character(1))
  if (length(names) <= limit) return(paste(names, collapse = ", "))
  paste0(paste(names[seq_len(limit)], collapse = ", "), ", and ", length(names) - limit, " more")
}

ScidrAnomalyMissingnessBlockReason <- function(observed, size, n_complete, n_started, block, labels) {
  paste0("Blank, although this record has ", observed, " of the ", size,
    " variables in a group that go missing together (", ScidrAnomalyMissingnessBlockText(block, labels),
    "), and ", n_complete, " of the ", n_started, " records carrying any of that group have all ", size, ".")
}
