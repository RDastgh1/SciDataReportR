# Shared wording for adjusted significance in matrix plots.
ScidrAdjustedSignificanceSubtitle <- function(scope, method = "BH", family = NULL) {
  scope_label <- switch(scope,
    matrix = "across the whole matrix",
    per_outcome = "within each outcome",
    per_predictor = "within each predictor",
    per_group = "within each group-vs-referent contrast",
    per_variable = "within each variable",
    stop("Unknown significance subtitle scope: ", scope, call. = FALSE)
  )
  method_label <- if (method %in% c("BH", "fdr")) {
    "FDR-adjusted p-values (Benjamini-Hochberg)"
  } else {
    paste0(method, "-adjusted p-values")
  }
  paste0("Significance: ", method_label, "\nCorrection: ", scope_label,
    if (!is.null(family)) paste0("; ", family))
}
