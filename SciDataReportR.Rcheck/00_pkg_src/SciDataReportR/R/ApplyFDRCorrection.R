#' Apply multiple-comparison correction across a p-value matrix
#'
#' Central helper used by the SciDataReportR heatmap/matrix family to apply
#' multiple-comparison correction (FDR by default) with a selectable scope.
#' All non-finite p-values (`NA`, `NaN`) are left untouched and excluded from
#' the correction, matching the long-standing behavior of the plotting
#' functions that now delegate to this helper.
#'
#' @param pmat A numeric matrix or data frame of p-values, or a plain numeric
#'   vector. Matrices and data frames keep their dimensions and dimnames.
#' @param fdr_scope Either `"matrix"` (default) or `"per_outcome"`.
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
#' @export
ApplyFDRCorrection <- function(pmat,
                               fdr_scope = c("matrix", "per_outcome"),
                               outcome_margin = 2,
                               method = "fdr",
                               outcome_ids = NULL) {
  fdr_scope <- match.arg(fdr_scope)
  if (!outcome_margin %in% c(1, 2)) {
    stop("outcome_margin must be 1 (outcomes are rows) or 2 (outcomes are columns).")
  }

  adjust_vec <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method)
    out
  }

  # ---- vector input --------------------------------------------------------
  if (is.null(dim(pmat))) {
    p <- as.numeric(pmat)
    if (fdr_scope == "matrix") {
      res <- adjust_vec(p)
    } else {
      if (is.null(outcome_ids)) {
        stop("outcome_ids is required when pmat is a vector and fdr_scope = 'per_outcome'.")
      }
      if (length(outcome_ids) != length(p)) {
        stop("outcome_ids must have the same length as pmat.")
      }
      res <- rep(NA_real_, length(p))
      for (idx in split(seq_along(p), as.character(outcome_ids))) {
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
  res <- m
  if (fdr_scope == "matrix") {
    res[] <- adjust_vec(as.vector(m))
  } else if (outcome_margin == 2) {
    for (j in seq_len(ncol(m))) res[, j] <- adjust_vec(m[, j])
  } else {
    for (i in seq_len(nrow(m))) res[i, ] <- adjust_vec(m[i, ])
  }
  if (is_df) {
    res <- as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
  }
  res
}
