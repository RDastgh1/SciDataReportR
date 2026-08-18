#' Apply multiple-comparison correction across a p-value matrix
#'
#' Central helper used by the SciDataReportR heatmap/matrix family to apply
#' multiple-comparison correction (FDR by default) with a selectable scope.
#' All non-finite p-values (`NA`, `NaN`) are left untouched and excluded from
#' the correction, matching the long-standing behavior of the plotting
#' functions that now delegate to this helper.
#'
#' @section Symmetric matrices:
#' When a variable set is correlated against itself, the resulting p-value
#' matrix is symmetric: every pair appears twice, above and below the diagonal,
#' and the diagonal holds self-comparisons that were never tested. Correcting
#' across every cell therefore describes a family of `n * (n - 1)` tests when
#' only `n * (n - 1) / 2` were actually run - 90 instead of 45 for ten
#' variables.
#'
#' How much that matters depends on the method, and on what sits on the
#' diagonal:
#'
#' * Benjamini-Hochberg (`method = "fdr"`, the default) is invariant to exact
#'   duplication, because doubling both a p-value's rank and the family size
#'   cancels out. Adjusted values are unchanged.
#' * `"bonferroni"`, `"holm"`, and `"hochberg"` are not invariant. Counting
#'   each pair twice made every adjusted p-value exactly twice as large as it
#'   should be.
#' * A diagonal carrying real values is the damaging case in every method. A
#'   self-correlation p-value of 0 enters the family as the most significant
#'   test there is, which drags the whole ranking down and makes the off-diagonal
#'   results look *stronger* than they are.
#'
#' `symmetric = "auto"` (the default) detects symmetry and corrects the unique
#' pairs once, then mirrors the adjusted values back into the lower triangle so
#' the matrix stays symmetric. The diagonal is excluded from the family and
#' returned as `NA` unless `include_diagonal = TRUE`. Set `symmetric = FALSE`
#' to force the old whole-matrix behavior, or `symmetric = TRUE` to require
#' symmetric handling and error out if `pmat` is not symmetric.
#'
#' Detection ignores the diagonal, since the heatmap functions routinely blank
#' it before correcting, and requires the two triangles to agree on both their
#' values and their missing cells. An asymmetric matrix is never affected.
#'
#' Under `"per_outcome"` or `"per_predictor"` scope each row or column is
#' already a family of unique tests, so symmetric input only has its diagonal
#' excluded. The result of a per-group correction on a symmetric matrix is not
#' itself symmetric.
#'
#' @param pmat A numeric matrix or data frame of p-values, or a plain numeric
#'   vector. Matrices and data frames keep their dimensions and dimnames.
#' @param fdr_scope Either `"matrix"` (default), `"per_outcome"`, or `"per_predictor"`.
#'   `"matrix"` corrects across all p-values at once (one family).
#'   `"per_outcome"` corrects separately within each outcome's p-values:
#'   for matrix input, groups run along `outcome_margin`; for vector input,
#'   groups are defined by `outcome_ids`.
#' @param outcome_margin For matrix/data-frame input with
#'   `fdr_scope = "per_outcome"`: `2` (default) if outcomes are columns, `1`
#'   if outcomes are rows. Ignored for `"matrix"` scope and for vector input.
#' @param method Correction method passed to [stats::p.adjust()]. Default
#'   `"fdr"` (Benjamini-Hochberg).
#' @param outcome_ids Optional vector (same length as `pmat`) identifying the
#'   outcome each p-value belongs to. Only used - and then required - when
#'   `pmat` is a vector and `fdr_scope = "per_outcome"`. This is how the
#'   long-format table functions (for example [PlotPhiHeatmap()] or
#'   [PlotChiSqCovar()]) group their p-values by outcome.
#' @param predictor_ids Optional vector (same length as `pmat`) identifying the
#'   predictor each p-value belongs to. Only used - and then required - when
#'   `pmat` is a vector and `fdr_scope = "per_predictor"`.
#' @param symmetric How to treat a square matrix whose two triangles hold the
#'   same p-values. `"auto"` (default) detects symmetry and corrects each pair
#'   once; `TRUE` requires it (and errors if `pmat` is not symmetric); `FALSE`
#'   restores the whole-matrix behavior that counts each pair twice. Ignored
#'   for vector input. See the Symmetric matrices section.
#' @param include_diagonal Logical; for symmetric input, whether the diagonal
#'   (self-comparisons) joins the family being corrected. Default `FALSE`,
#'   which excludes it and returns `NA` on the diagonal.
#'
#' @return An object of the same shape as `pmat` (matrix, data frame, or
#'   vector) containing adjusted p-values. Non-finite entries remain `NA`.
#'
#' @examples
#' pm <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
#'              nrow = 2,
#'              dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3")))
#'
#' # One family across the whole matrix (classic behavior)
#' ApplyFDRCorrection(pm)
#'
#' # Correct each outcome (column) separately
#' ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 2)
#'
#' # Vector input with explicit outcome grouping
#' ApplyFDRCorrection(c(0.01, 0.04, 0.02, 0.03),
#'                    fdr_scope = "per_outcome",
#'                    outcome_ids = c("y1", "y1", "y2", "y2"))
#'
#' # A symmetric matrix: six filled cells, but only three pairs tested
#' vars <- c("mmse", "trails", "digit_span")
#' pm_sym <- matrix(NA_real_, 3, 3, dimnames = list(vars, vars))
#' pm_sym[lower.tri(pm_sym)] <- c(0.001, 0.020, 0.040)
#' pm_sym[upper.tri(pm_sym)] <- t(pm_sym)[upper.tri(pm_sym)]
#' pm_sym
#'
#' # Detected automatically: three tests in the family, diagonal excluded
#' ApplyFDRCorrection(pm_sym)
#'
#' # Bonferroni, corrected once per pair and then against all six cells
#' ApplyFDRCorrection(pm_sym, method = "bonferroni")
#' ApplyFDRCorrection(pm_sym, method = "bonferroni", symmetric = FALSE)
#'
#' @export
ApplyFDRCorrection <- function(pmat,
                               fdr_scope = c("matrix", "per_outcome", "per_predictor"),
                               outcome_margin = 2,
                               method = "fdr",
                               outcome_ids = NULL,
                               predictor_ids = NULL,
                               symmetric = "auto",
                               include_diagonal = FALSE) {
  fdr_scope <- match.arg(fdr_scope)
  if (!outcome_margin %in% c(1, 2)) {
    stop("outcome_margin must be 1 (outcomes are rows) or 2 (outcomes are columns).")
  }
  if (!(identical(symmetric, "auto") || isTRUE(symmetric) || isFALSE(symmetric))) {
    stop("symmetric must be \"auto\", TRUE, or FALSE.")
  }
  if (!isTRUE(include_diagonal) && !isFALSE(include_diagonal)) {
    stop("include_diagonal must be TRUE or FALSE.")
  }

  adjust_vec <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method)
    out
  }

  # Detects the "same variables down both margins" case. The diagonal is
  # ignored: heatmap callers routinely blank it out before correcting.
  is_symmetric_pmat <- function(m) {
    if (nrow(m) != ncol(m) || nrow(m) < 2) return(FALSE)
    rn <- rownames(m)
    cn <- colnames(m)
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) return(FALSE)
    upper <- m[upper.tri(m)]
    lower <- t(m)[upper.tri(m)]
    missing_upper <- !is.finite(upper)
    missing_lower <- !is.finite(lower)
    if (any(missing_upper != missing_lower)) return(FALSE)
    paired <- !missing_upper
    # An entirely empty off-diagonal carries no evidence of symmetry.
    if (!any(paired)) return(FALSE)
    isTRUE(all.equal(upper[paired], lower[paired], tolerance = 1e-8))
  }

  # ---- vector input --------------------------------------------------------
  if (is.null(dim(pmat))) {
    p <- as.numeric(pmat)
    if (fdr_scope == "matrix") {
      res <- adjust_vec(p)
    } else {
      ids <- if (fdr_scope == "per_outcome") outcome_ids else predictor_ids
      id_name <- if (fdr_scope == "per_outcome") "outcome_ids" else "predictor_ids"
      if (is.null(ids)) {
        stop(id_name, " is required when pmat is a vector and fdr_scope = '", fdr_scope, "'.")
      }
      if (length(ids) != length(p)) {
        stop(id_name, " must have the same length as pmat.")
      }
      res <- rep(NA_real_, length(p))
      for (idx in split(seq_along(p), as.character(ids))) {
        res[idx] <- adjust_vec(p[idx])
      }
    }
    names(res) <- names(pmat)
    return(res)
  }

  # ---- matrix / data frame input -------------------------------------------
  is_df <- is.data.frame(pmat)
  m <- as.matrix(pmat)
  if (!is.numeric(m)) {
    stop("pmat must contain numeric p-values.")
  }
  is_symmetric <- if (identical(symmetric, "auto")) {
    is_symmetric_pmat(m)
  } else if (isTRUE(symmetric)) {
    if (!is_symmetric_pmat(m)) {
      stop("symmetric = TRUE, but pmat is not a symmetric p-value matrix.")
    }
    TRUE
  } else {
    FALSE
  }

  # A symmetric matrix holds each pair twice, plus untested self-comparisons on
  # the diagonal. Correct the unique pairs once, then restore the symmetry.
  if (is_symmetric) {
    if (!include_diagonal) diag(m) <- NA_real_
    if (fdr_scope == "matrix") {
      keep <- upper.tri(m, diag = include_diagonal)
      res <- m
      res[] <- NA_real_
      res[keep] <- adjust_vec(m[keep])
      lower <- lower.tri(res)
      res[lower] <- t(res)[lower]
      if (is_df) {
        res <- as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
      }
      return(res)
    }
  }

  res <- m
  if (fdr_scope == "matrix") {
    res[] <- adjust_vec(as.vector(m))
  } else if (fdr_scope == "per_outcome" && outcome_margin == 2) {
    for (j in seq_len(ncol(m))) res[, j] <- adjust_vec(m[, j])
  } else if (fdr_scope == "per_outcome") {
    for (i in seq_len(nrow(m))) res[i, ] <- adjust_vec(m[i, ])
  } else if (outcome_margin == 2) {
    for (i in seq_len(nrow(m))) res[i, ] <- adjust_vec(m[i, ])
  } else {
    for (j in seq_len(ncol(m))) res[, j] <- adjust_vec(m[, j])
  }
  if (is_df) {
    res <- as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
  }
  res
}
