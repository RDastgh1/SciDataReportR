#' Read a scientific data file with optional inspection
#'
#' `ReadSciData()` imports common scientific data file formats while preserving
#' original column names and labels. For Excel files, it can call `InspectFile()`
#' before import to flag common issues such as multiple sheets, metadata rows,
#' and color-coded formatting.
#'
#' @param path Path to the file.
#' @param sheet Sheet name or index for Excel files.
#' @param header_row Row containing column names for Excel files. If `NULL`,
#'   `ReadSciData()` may use the probable header row detected by `InspectFile()`.
#' @param col_names Passed to file readers where applicable. For Excel files,
#'   use `TRUE`, `FALSE`, or a character vector.
#' @param range Optional Excel range passed to `readxl::read_excel()`.
#' @param inspect Logical. If `TRUE`, inspect the file before importing.
#' @param print_inspection Logical. If `TRUE`, print inspection results when
#'   issues are detected.
#' @param strict Logical. If `TRUE`, stop when inspection detects potential issues.
#' @param guess_max Maximum rows used for type guessing where supported.
#' @param delim Delimiter for `.txt` files. Defaults to tab.
#' @param ... Additional arguments passed to the underlying reader.
#'
#' @return Imported data object. Inspection metadata is attached as the
#'   `scidata_inspection` attribute when `inspect = TRUE`.
#' @export
ReadSciData <- function(
  path,
  sheet = NULL,
  header_row = NULL,
  col_names = TRUE,
  range = NULL,
  inspect = TRUE,
  print_inspection = TRUE,
  strict = FALSE,
  guess_max = 10000,
  delim = NULL,
  ...
) {
  if (length(path) != 1) {
    stop("`path` must be a single file path.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  inspection <- NULL

  if (isTRUE(inspect)) {
    inspection <- InspectFile(
      path = path,
      sheet = sheet,
      quiet = TRUE
    )

    if (isTRUE(print_inspection) && length(inspection$issues) > 0) {
      InspectFile(
        path = path,
        sheet = sheet,
        quiet = FALSE
      )
    }

    if (isTRUE(strict) && length(inspection$issues) > 0) {
      stop(
        "InspectFile() found potential import issues. ",
        "Review `InspectFile(path)` or set `strict = FALSE`.",
        call. = FALSE
      )
    }
  }

  if (ext %in% c("xlsx", "xls")) {
    excel_sheet <- sheet

    if (is.null(excel_sheet) && !is.null(inspection$selected_sheet)) {
      excel_sheet <- inspection$selected_sheet
    }

    if (is.null(excel_sheet)) {
      excel_sheet <- 1
    }

    excel_skip <- 0

    if (!is.null(header_row)) {
      excel_skip <- header_row - 1
    } else if (
      isTRUE(inspect) &&
        !is.null(inspection$probable_header_row) &&
        inspection$probable_header_row > 1
    ) {
      excel_skip <- inspection$probable_header_row - 1
    }

    if (excel_skip < 0) {
      stop("`header_row` must be 1 or greater.", call. = FALSE)
    }
  }

  out <- switch(
    ext,

    csv = readr::read_csv(
      file = path,
      col_names = col_names,
      guess_max = guess_max,
      name_repair = "minimal",
      show_col_types = FALSE,
      ...
    ),

    tsv = readr::read_tsv(
      file = path,
      col_names = col_names,
      guess_max = guess_max,
      name_repair = "minimal",
      show_col_types = FALSE,
      ...
    ),

    txt = readr::read_delim(
      file = path,
      delim = if (is.null(delim)) "\t" else delim,
      col_names = col_names,
      guess_max = guess_max,
      name_repair = "minimal",
      show_col_types = FALSE,
      ...
    ),

    xlsx = readxl::read_excel(
      path = path,
      sheet = excel_sheet,
      skip = excel_skip,
      col_names = col_names,
      range = range,
      guess_max = guess_max,
      .name_repair = "minimal",
      ...
    ),

    xls = readxl::read_excel(
      path = path,
      sheet = excel_sheet,
      skip = excel_skip,
      col_names = col_names,
      range = range,
      guess_max = guess_max,
      .name_repair = "minimal",
      ...
    ),

    rds = readRDS(path),

    rda = {
      env <- new.env(parent = emptyenv())
      object_names <- load(path, envir = env)

      if (length(object_names) != 1) {
        stop(
          "The .rda file contains multiple objects: ",
          paste(object_names, collapse = ", "),
          ". ReadSciData() requires exactly one object for automatic import.",
          call. = FALSE
        )
      }

      env[[object_names]]
    },

    rdata = {
      env <- new.env(parent = emptyenv())
      object_names <- load(path, envir = env)

      if (length(object_names) != 1) {
        stop(
          "The .RData file contains multiple objects: ",
          paste(object_names, collapse = ", "),
          ". ReadSciData() requires exactly one object for automatic import.",
          call. = FALSE
        )
      }

      env[[object_names]]
    },

    sav = haven::read_sav(
      file = path,
      ...
    ),

    dta = haven::read_dta(
      file = path,
      ...
    ),

    sas7bdat = haven::read_sas(
      data_file = path,
      ...
    ),

    xpt = haven::read_xpt(
      file = path,
      ...
    ),

    parquet = arrow::read_parquet(
      file = path,
      as_data_frame = TRUE,
      ...
    ),

    feather = arrow::read_feather(
      file = path,
      as_data_frame = TRUE,
      ...
    ),

    json = jsonlite::read_json(
      path = path,
      simplifyVector = TRUE,
      ...
    ),

    stop("Unsupported file type: .", ext, " for file ", path, call. = FALSE)
  )

  attr(out, "scidata_source") <- normalizePath(path, winslash = "/", mustWork = FALSE)
  attr(out, "scidata_inspection") <- inspection

  out
}
