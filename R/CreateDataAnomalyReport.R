#' Create shareable data anomaly review files in one call
#'
#' Run anomaly detection, build the review queue, and write the two artifacts a
#' collaborator needs: a standalone interactive HTML explorer they can open in a
#' browser, and an Excel workbook they can work through issue by issue and send
#' back. Source data is never modified.
#' @param data A data frame with uniquely named columns.
#' @param variables Model input columns, as in [DetectDataAnomalies()]. NULL
#'   (default) selects all supported columns except ID and strata.
#' @param id_var Optional participant ID column name, retained for review.
#' @param stratify_by Optional character vector defining independent fits.
#' @param file Destination `.html` path for the explorer, in an existing
#'   directory. NULL (default) writes `DataAnomalyReview_<date>.html` into
#'   `output_dir`. The workbook replaces the extension with `.xlsx`.
#' @param output_dir Directory used when `file` is NULL. Default the working
#'   directory.
#' @param longitudinal NULL (default) for cross-sectional screening only, or a
#'   named list of [DetectLongitudinalAnomalies()] arguments, at least
#'   `time_var` and `variables`, run on the same `data` and added to the same
#'   queue as a second analysis.
#' @param Title Dashboard title, also the report heading.
#' @param workbook Write the Excel review workbook? Default TRUE. Requires
#'   optional openxlsx and jsonlite.
#' @param explorer Write the standalone HTML explorer? Default TRUE. Requires
#'   optional DT and htmltools.
#' @param codebook Optional presentation codebook as in [DetectDataAnomalies()].
#' @param top_n Positive integer limit per method and stratum for unclassified
#'   score-only candidates added to the queue. Default 10.
#' @param layout Workbook layout. `"compact"` (the default) writes one visible
#'   sheet of one-line rows carrying the columns a reviewer scans and the
#'   editable decision fields; the supporting sheets and the remaining columns
#'   are hidden rather than dropped, so decisions still import. `"full"` shows
#'   every sheet and column with wrapped multi-line cells.
#' @param review_variables Optional character vector of screened column names
#'   limiting which findings the HTML explorer lists; whole-record findings are
#'   always listed. Detection and the Excel workbook are unaffected, so the
#'   workbook remains the complete record of findings.
#' @param overwrite Allow replacement of existing files? Default FALSE.
#' @param ... Passed to [DetectDataAnomalies()], for example `methods`,
#'   `thresholds`, `flag_percentile`, `method_args`, or the temporal arguments.
#'   The same `...` is not forwarded to the longitudinal detector; supply those
#'   settings inside `longitudinal`.
#' @details Missing optional packages are reported once before any modelling
#'   runs, rather than after a long detection. Detection, queue construction and
#'   file writing reuse [DetectDataAnomalies()],
#'   [DetectLongitudinalAnomalies()], [CreateAnomalyReviewQueue()],
#'   [SaveDataAnomalyExplorer()] and [ExportAnomalyReviewQueue()]; this function
#'   adds no separate detection behaviour. Flags in the written files are review
#'   candidates, not confirmed errors.
#' @return Invisibly, a list with `Queue` (the `AnomalyReviewQueue`), `Results`
#'   (named list of detection results) and `Files` (named character vector of
#'   the paths written). Each written file is also reported with a message.
#' @family data anomalies
#' @seealso [DetectDataAnomalies()], [ExploreDataAnomalies()],
#'   [SummarizeAnomalyDecisions()]
#' @examples
#' set.seed(24)
#' df_Example <- data.frame(ID = 1:60, Visit = as.Date("2024-01-01") + c(0:58, 900))
#' if (requireNamespace("DT", quietly = TRUE)) {
#'   review <- CreateDataAnomalyReport(df_Example, "Visit", id_var = "ID",
#'     methods = "dates", workbook = FALSE, output_dir = tempdir())
#'   review$Files
#' }
#' @export
CreateDataAnomalyReport <- function(data, variables = NULL, id_var = NULL,
    stratify_by = NULL, file = NULL, output_dir = ".", longitudinal = NULL,
    Title = "Data anomaly review", workbook = TRUE, explorer = TRUE,
    codebook = NULL, top_n = 10, overwrite = FALSE, review_variables = NULL,
    layout = c("compact", "full"), ...) {
  layout <- match.arg(layout)
  for (arg in c("workbook", "explorer", "overwrite")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(arg, " must be TRUE or FALSE.", call. = FALSE)
  }
  if (!explorer && !workbook) stop("Enable at least one of explorer or workbook.", call. = FALSE)
  if (!is.character(Title) || length(Title) != 1L || is.na(Title)) stop("Title must be one string.", call. = FALSE)
  if (!is.null(longitudinal) && (!is.list(longitudinal) || !length(longitudinal) ||
      is.null(names(longitudinal)) || any(!nzchar(names(longitudinal))) || anyDuplicated(names(longitudinal))))
    stop("longitudinal must be NULL or a uniquely named list of DetectLongitudinalAnomalies() arguments.", call. = FALSE)
  if (!is.null(longitudinal) && any(c("data", "codebook") %in% names(longitudinal)))
    stop("longitudinal must not set data or codebook; those are shared with the cross-sectional screen.", call. = FALSE)
  if (is.null(file)) {
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir)) stop("output_dir must be one path.", call. = FALSE)
    if (!dir.exists(output_dir)) stop("Directory does not exist: ", output_dir, call. = FALSE)
    file <- file.path(output_dir, paste0("DataAnomalyReview_", format(Sys.Date(), "%Y_%m_%d"), ".html"))
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) stop("file must be one path.", call. = FALSE)
  if (!dir.exists(dirname(file))) stop("Directory does not exist: ", dirname(file), call. = FALSE)
  fn_Workbook <- paste0(tools::file_path_sans_ext(file), ".xlsx")
  for (destination in c(if (explorer) file, if (workbook) fn_Workbook))
    if (file.exists(destination) && !overwrite) stop("File exists; use overwrite = TRUE: ", destination, call. = FALSE)

  settings <- list(...)
  methods <- if ("methods" %in% names(settings)) settings$methods else
    eval(formals(DetectDataAnomalies)$methods)
  if (!is.null(longitudinal)) methods <- c(methods, if ("methods" %in% names(longitudinal))
    longitudinal$methods else eval(formals(DetectLongitudinalAnomalies)$methods))
  explain <- identical(tryCatch(settings$explanation_args$method, error = function(e) NULL), "shap")
  ScidrAnomalyRequireDependencies(ScidrAnomalyDependencies(methods = unique(methods),
    explorer = explorer, workbook = workbook, plots = explorer, explain = explain))

  results <- list()
  results$`Cross-sectional` <- do.call(DetectDataAnomalies, c(list(data = data, variables = variables,
    id_var = id_var, stratify_by = stratify_by, codebook = codebook), settings))
  if (!is.null(longitudinal)) {
    arguments <- longitudinal
    if (is.null(arguments$id_var)) arguments$id_var <- id_var
    if (is.null(arguments$stratify_by)) arguments$stratify_by <- stratify_by
    results$Longitudinal <- do.call(DetectLongitudinalAnomalies,
      c(list(data = data, codebook = codebook), arguments))
  }
  queue <- CreateAnomalyReviewQueue(results, top_n = top_n, codebook = codebook)

  files <- character()
  if (explorer) files["Explorer"] <- SaveDataAnomalyExplorer(queue, file = file,
    overwrite = overwrite, Title = Title, variables = review_variables)
  if (workbook) {
    ExportAnomalyReviewQueue(queue, fn_Workbook, overwrite = overwrite, layout = layout)
    files["Workbook"] <- fn_Workbook
    message("Wrote workbook: ", fn_Workbook)
  }
  message(nrow(queue$Issues), " issue(s) for review. Flags are review candidates, not confirmed errors.")
  invisible(list(Queue = queue, Results = results, Files = files))
}
