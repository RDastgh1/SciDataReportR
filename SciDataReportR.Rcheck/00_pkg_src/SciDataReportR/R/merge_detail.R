#' Print a plain-text detail report for one safe_merge result
#'
#' Print static `knitr::kable()` tables (no plots, no htmlwidgets) for the key
#' diagnostic sections of a [safe_merge()] result: the validation checks,
#' unmatched key combinations on each side, overlapping non-key variables, and
#' suspicious duplicate-variable conflicts. Sections with nothing to show are
#' skipped entirely.
#'
#' In R Markdown or Quarto, use a chunk with `results = "asis"` so the kable
#' markdown renders as tables.
#'
#' @param m A list returned by [safe_merge()] (must contain a `validation`
#'   element from [ValidateMerge()] and a `log` tibble).
#' @param TopN Integer. Maximum number of rows shown for the left-only and
#'   right-only unmatched-key previews. Default is `10`.
#'
#' @return `m`, invisibly. Called for its printed output.
#'
#' @seealso [safe_merge()], [merge_summary_table()], [ExploreMergeValidation()]
#'
#' @examples
#' baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
#' labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
#' m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
#' merge_detail(m)
#'
#' @export
merge_detail <- function(m, TopN = 10) {

  if (!is.list(m) || !all(c("validation", "log") %in% names(m))) {
    stop("m must be a list returned by safe_merge(), with validation and log elements.")
  }

  if (!is.numeric(TopN) || length(TopN) != 1 || TopN < 1) {
    stop("TopN must be a single positive number.")
  }

  TopN <- as.integer(TopN)
  v <- m$validation
  merge_name <- as.character(m$log$Merge[1])

  print_section <- function(df, caption) {
    print(knitr::kable(df, caption = caption))
    cat("\n")
  }

  print_section(
    v$Checks,
    paste0(merge_name, ": validation checks")
  )

  if (nrow(v$IDCoverage$LeftOnly) > 0) {
    print_section(
      utils::head(v$IDCoverage$LeftOnly, TopN),
      paste0(
        merge_name, ": keys only in left data (showing up to ", TopN,
        " of ", nrow(v$IDCoverage$LeftOnly), ")"
      )
    )
  }

  if (nrow(v$IDCoverage$RightOnly) > 0) {
    print_section(
      utils::head(v$IDCoverage$RightOnly, TopN),
      paste0(
        merge_name, ": keys only in right data (showing up to ", TopN,
        " of ", nrow(v$IDCoverage$RightOnly), ")"
      )
    )
  }

  if (nrow(v$OverlappingVariables) > 0) {
    print_section(
      v$OverlappingVariables,
      paste0(merge_name, ": overlapping non-key variables")
    )
  }

  if (nrow(v$SuspiciousConflicts) > 0) {
    print_section(
      v$SuspiciousConflicts,
      paste0(merge_name, ": suspicious duplicate-variable conflicts")
    )
  }

  invisible(m)
}
