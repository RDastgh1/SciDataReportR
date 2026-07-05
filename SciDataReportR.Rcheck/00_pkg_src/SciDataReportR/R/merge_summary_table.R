#' Combine safe_merge logs into a single summary table
#'
#' Bind the one-row `$log` tibbles produced by [safe_merge()] into a single
#' table, one row per merge, optionally filtered to merges that did not pass
#' cleanly. Useful as an end-of-pipeline rollup after a sequence of merges.
#'
#' @param merge_log A list of `$log` tibbles from [safe_merge()] results.
#'   Full [safe_merge()] result objects are also accepted (the `$log` element
#'   is extracted automatically), as is a single log tibble.
#' @param flagged_only Logical. If `TRUE`, only rows with
#'   `Status != "PASS"` are returned. Default is `FALSE`.
#'
#' @return A tibble with one row per merge and the same columns as a
#'   [safe_merge()] `$log` tibble.
#'
#' @seealso [safe_merge()]
#'
#' @examples
#' baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
#' labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
#' vitals <- data.frame(id = 1:4, sbp = c(120, 135, 118, 141))
#' m1 <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
#' m2 <- safe_merge(m1$data, vitals, by = "id", name = "+ vitals")
#' merge_summary_table(list(m1$log, m2$log))
#'
#' @export
merge_summary_table <- function(merge_log, flagged_only = FALSE) {

  if (is.data.frame(merge_log)) {
    merge_log <- list(merge_log)
  }

  if (!is.list(merge_log) || length(merge_log) == 0) {
    stop("merge_log must be a non-empty list of safe_merge() $log tibbles.")
  }

  logs <- lapply(merge_log, function(x) {
    if (is.data.frame(x)) {
      x
    } else if (is.list(x) && is.data.frame(x$log)) {
      x$log
    } else {
      stop(
        "Each element of merge_log must be a safe_merge() $log tibble ",
        "(or a full safe_merge() result)."
      )
    }
  })

  out <- dplyr::bind_rows(logs)

  if (isTRUE(flagged_only)) {
    out <- dplyr::filter(out, Status != "PASS")
  }

  out
}
