#' Read a scientific data file with optional inspection
#'
#' `ReadSciData()` imports common scientific data file formats while preserving
#' original column names and labels as much as possible. It optionally calls
#' `InspectFile()` before import to flag common issues such as multiple sheets,
#' metadata rows, unnamed columns, duplicate column names, and formatting.
#'
#' After import, `ReadSciData()` repairs only names that make the data frame
#' difficult or impossible to use in tidyverse workflows. Blank or `NA` names are
#' renamed to `...unnamed_POSITION`, and duplicate names are made unique only when
#' needed. Existing valid names are preserved.
#'
#' @param path Path to the file.
#' @param sheet Sheet name or index for Excel files.
#' @param header_row Row containing column names for Excel files. If `NULL`,
#'   `ReadSciData()` may use the probable header row detected by `InspectFile()`.
#' @param col_names Passed to file readers where applicable. For Excel files,
#'   use `TRUE`, `FALSE`, or a character vector.
#' @param range Optional Excel range passed to `readxl::read_excel()`.
#' @param inspect Logical. If `TRUE`, inspect the file before importing.
#' @param print_inspection Logical. If `TRUE`, print compact inspection results
#'   when issues are detected. Defaults to `interactive()`.
#' @param strict Logical. If `TRUE`, stop when inspection detects potential
#'   issues.
#' @param guess_max Maximum rows used for type guessing where supported.
#' @param delim Delimiter for `.txt` files. Defaults to tab.
#' @param repair_names Logical. If `TRUE`, repair blank and duplicate column
#'   names after import.
#' @param inspect_styles Logical. If `TRUE`, inspect Excel workbook formatting.
#'   This can be slower and noisier for large workbooks.
#' @param fast_delimited Logical. If `TRUE`, use `data.table::fread()` for
#'   delimited text files when available and when no extra reader arguments are
#'   supplied through `...`. This is usually much faster than `readr` for large
#'   `.csv`, `.tsv`, and `.txt` files.
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
  print_inspection = interactive(),
  strict = FALSE,
  guess_max = 10000,
  delim = NULL,
  repair_names = TRUE,
  inspect_styles = FALSE,
  fast_delimited = TRUE,
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
  reader_args <- list(...)

  if (!is.logical(fast_delimited) || length(fast_delimited) != 1 || is.na(fast_delimited)) {
    stop("`fast_delimited` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(inspect)) {
    inspection <- InspectFile(
      path = path,
      sheet = sheet,
      check_styles = inspect_styles,
      quiet = TRUE
    )

    if (isTRUE(print_inspection) && length(inspection$issues) > 0) {
      .PrintSciDataInspection(inspection)
    }

    if (isTRUE(strict) && length(inspection$issues) > 0) {
      stop(
        "InspectFile() found potential import issues. ",
        "Review `InspectFile(path)` or set `strict = FALSE`.",
        call. = FALSE
      )
    }
  }

  excel_sheet <- NULL
  excel_skip <- 0

  if (ext %in% c("xlsx", "xls")) {
    excel_sheet <- sheet

    if (is.null(excel_sheet) && !is.null(inspection$selected_sheet)) {
      excel_sheet <- inspection$selected_sheet
    }

    if (is.null(excel_sheet)) {
      excel_sheet <- 1
    }

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

    csv = .ReadSciDelimited(
      path = path,
      delim = ",",
      col_names = col_names,
      guess_max = guess_max,
      fast_delimited = fast_delimited,
      reader_args = reader_args
    ),

    tsv = .ReadSciDelimited(
      path = path,
      delim = "\t",
      col_names = col_names,
      guess_max = guess_max,
      fast_delimited = fast_delimited,
      reader_args = reader_args
    ),

    txt = .ReadSciDelimited(
      path = path,
      delim = if (is.null(delim)) "\t" else delim,
      col_names = col_names,
      guess_max = guess_max,
      fast_delimited = fast_delimited,
      reader_args = reader_args
    ),

    xlsx = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop(
          "Package `readxl` is required to read Excel files. ",
          "Install it with install.packages('readxl').",
          call. = FALSE
        )
      }

      readxl::read_excel(
        path = path,
        sheet = excel_sheet,
        skip = excel_skip,
        col_names = col_names,
        range = range,
        guess_max = guess_max,
        .name_repair = "minimal",
        ...
      )
    },

    xls = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop(
          "Package `readxl` is required to read Excel files. ",
          "Install it with install.packages('readxl').",
          call. = FALSE
        )
      }

      readxl::read_excel(
        path = path,
        sheet = excel_sheet,
        skip = excel_skip,
        col_names = col_names,
        range = range,
        guess_max = guess_max,
        .name_repair = "minimal",
        ...
      )
    },

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

    sav = haven::read_sav(file = path, ...),

    dta = haven::read_dta(file = path, ...),

    sas7bdat = haven::read_sas(data_file = path, ...),

    xpt = haven::read_xpt(file = path, ...),

    parquet = {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
          "Package `arrow` is required to read Parquet files. ",
          "Install it with install.packages('arrow').",
          call. = FALSE
        )
      }

      arrow::read_parquet(file = path, as_data_frame = TRUE, ...)
    },

    feather = {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
          "Package `arrow` is required to read Feather files. ",
          "Install it with install.packages('arrow').",
          call. = FALSE
        )
      }

      arrow::read_feather(file = path, as_data_frame = TRUE, ...)
    },

    json = {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop(
          "Package `jsonlite` is required to read JSON files. ",
          "Install it with install.packages('jsonlite').",
          call. = FALSE
        )
      }

      jsonlite::read_json(path = path, simplifyVector = TRUE, ...)
    },

    stop("Unsupported file type: .", ext, " for file ", path, call. = FALSE)
  )

  if (is.data.frame(out) && isTRUE(repair_names)) {
    name_repair <- .RepairSciDataNames(names(out))
    names(out) <- name_repair$repaired_names

    if (!is.null(inspection)) {
      inspection$import_name_repair <- name_repair
    }

    if (length(name_repair$blank_positions) > 0) {
      message(
        "ReadSciData(): ",
        length(name_repair$blank_positions),
        " unnamed column",
        ifelse(length(name_repair$blank_positions) == 1, "", "s"),
        " renamed to ",
        paste(name_repair$repaired_names[name_repair$blank_positions], collapse = ", "),
        "."
      )
    }

    if (length(name_repair$duplicate_changed_positions) > 0) {
      message(
        "ReadSciData(): duplicate column names repaired at column(s) ",
        paste(name_repair$duplicate_changed_positions, collapse = ", "),
        "."
      )
    }
  }

  attr(out, "scidata_source") <- normalizePath(path, winslash = "/", mustWork = FALSE)
  attr(out, "scidata_inspection") <- inspection

  out
}

.RepairSciDataNames <- function(x) {
  original_names <- x
  repaired_names <- x

  blank_positions <- which(is.na(repaired_names) | repaired_names == "")

  if (length(blank_positions) > 0) {
    repaired_names[blank_positions] <- paste0("...unnamed_", blank_positions)
  }

  before_unique <- repaired_names

  if (any(duplicated(repaired_names))) {
    repaired_names <- make.unique(repaired_names, sep = "_")
  }

  list(
    original_names = original_names,
    repaired_names = repaired_names,
    blank_positions = blank_positions,
    duplicate_changed_positions = which(before_unique != repaired_names),
    changed_positions = which(original_names != repaired_names)
  )
}

.ReadSciDelimited <- function(
  path,
  delim,
  col_names,
  guess_max,
  fast_delimited,
  reader_args
) {
  can_use_fread <- isTRUE(fast_delimited) &&
    requireNamespace("data.table", quietly = TRUE) &&
    length(reader_args) == 0 &&
    (isTRUE(col_names) || is.character(col_names))

  if (can_use_fread) {
    header <- isTRUE(col_names)
    out <- data.table::fread(
      file = path,
      sep = delim,
      header = header,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    )

    if (is.character(col_names)) {
      names(out) <- col_names
    }

    return(tibble::as_tibble(out, .name_repair = "minimal"))
  }

  read_fun <- if (identical(delim, ",")) {
    readr::read_csv
  } else if (identical(delim, "\t")) {
    readr::read_tsv
  } else {
    readr::read_delim
  }

  args <- c(
    list(
      file = path,
      col_names = col_names,
      guess_max = guess_max,
      name_repair = "minimal",
      show_col_types = FALSE
    ),
    if (!identical(delim, ",") && !identical(delim, "\t")) list(delim = delim),
    reader_args
  )

  do.call(read_fun, args)
}
