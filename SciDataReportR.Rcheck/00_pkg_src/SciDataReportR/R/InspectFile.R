#' Inspect a scientific data file before import
#'
#' `InspectFile()` checks common import risks before a file is read into an
#' analysis workflow. It is especially useful for Excel workbooks, where
#' multiple sheets, metadata rows, unnamed columns, duplicate column names, and
#' workbook formatting can change how the data should be interpreted.
#'
#' This function does not modify the file or the imported data. It only reports
#' what it finds.
#'
#' @param path Path to the file.
#' @param sheet Sheet name or index for Excel files. If `NULL`, the first sheet
#'   is inspected.
#' @param preview_rows Number of rows to preview when checking header and column
#'   name issues.
#' @param check_styles Logical. For `.xlsx` files, check whether workbook styles
#'   or formatting exist. This can be slower for large workbooks.
#' @param check_sheets Logical. For Excel files, check whether the workbook has
#'   multiple sheets.
#' @param check_header Logical. For Excel files, attempt to detect whether the
#'   header row is not row 1.
#' @param quiet Logical. If `FALSE`, print a compact inspection summary.
#'
#' @return Invisibly returns a list containing file inspection metadata.
#' @export
InspectFile <- function(
  path,
  sheet = NULL,
  preview_rows = 20,
  check_styles = TRUE,
  check_sheets = TRUE,
  check_header = TRUE,
  quiet = FALSE
) {
  if (length(path) != 1) {
    stop("`path` must be a single file path.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  file_info <- file.info(path)

  supported_extensions <- c(
    "csv", "tsv", "txt",
    "xlsx", "xls",
    "rds", "rda", "rdata",
    "sav", "dta", "sas7bdat", "xpt",
    "parquet", "feather",
    "json"
  )

  issues <- character()
  notes <- character()

  sheets <- NULL
  selected_sheet <- NULL
  probable_header_row <- NULL
  first_nonempty_row <- NULL
  preview <- NULL
  preview_names <- NULL
  blank_columns <- integer()
  duplicate_column_names <- character()
  styles_detected <- NA

  if (!ext %in% supported_extensions) {
    issues <- c(issues, paste0("Unsupported or unrecognized file extension: .", ext))
  }

  if (ext %in% c("xlsx", "xls")) {
    sheets <- tryCatch(
      readxl::excel_sheets(path),
      error = function(e) {
        issues <<- c(issues, paste0("Could not read Excel sheet names: ", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(sheets) && length(sheets) > 0) {
      if (is.null(sheet)) {
        selected_sheet <- sheets[[1]]
      } else if (is.numeric(sheet)) {
        if (length(sheet) != 1 || is.na(sheet)) {
          issues <- c(issues, "`sheet` must be a single sheet name or sheet index.")
          selected_sheet <- sheets[[1]]
        } else if (sheet < 1 || sheet > length(sheets)) {
          issues <- c(
            issues,
            paste0("`sheet = ", sheet, "` is outside the available sheet range 1-", length(sheets), ".")
          )
          selected_sheet <- sheets[[1]]
        } else {
          selected_sheet <- sheets[[sheet]]
        }
      } else {
        selected_sheet <- as.character(sheet)
        if (!selected_sheet %in% sheets) {
          issues <- c(issues, paste0("Requested sheet `", selected_sheet, "` was not found in the workbook."))
        }
      }

      if (isTRUE(check_sheets) && length(sheets) > 1 && is.null(sheet)) {
        issues <- c(
          issues,
          paste0(
            "Workbook has multiple sheets and no `sheet =` was supplied: ",
            paste(sheets, collapse = ", "),
            "."
          )
        )
      } else if (isTRUE(check_sheets) && length(sheets) > 1 && !is.null(sheet)) {
        notes <- c(
          notes,
          paste0(
            "Workbook has multiple sheets, and `sheet = ",
            selected_sheet,
            "` was selected."
          )
        )
      }
    }

    if (!is.null(selected_sheet) && !is.null(sheets) && selected_sheet %in% sheets) {
      preview <- tryCatch(
        readxl::read_excel(
          path = path,
          sheet = selected_sheet,
          col_names = FALSE,
          n_max = preview_rows,
          col_types = "text",
          .name_repair = "minimal"
        ),
        error = function(e) {
          issues <<- c(issues, paste0("Could not preview Excel sheet: ", conditionMessage(e)))
          NULL
        }
      )
    }

    if (isTRUE(check_header) && !is.null(preview) && nrow(preview) > 0) {
      row_scores <- rep(0, nrow(preview))

      for (i in seq_len(nrow(preview))) {
        vals <- unlist(preview[i, , drop = FALSE], use.names = FALSE)
        vals <- trimws(as.character(vals))
        vals <- vals[!is.na(vals) & vals != ""]

        if (length(vals) == 0) {
          row_scores[[i]] <- 0
        } else {
          unique_ratio <- length(unique(vals)) / length(vals)
          numeric_like <- grepl("^[-+]?[0-9,.]+%?$", vals)
          text_ratio <- mean(!numeric_like)
          non_empty_score <- log1p(length(vals))

          row_scores[[i]] <- unique_ratio + text_ratio + non_empty_score
        }
      }

      if (any(row_scores > 0)) {
        probable_header_row <- which.max(row_scores)
        first_nonempty_row <- which(row_scores > 0)[[1]]

        if (probable_header_row > 1) {
          issues <- c(
            issues,
            paste0(
              "Probable header row appears to be row ",
              probable_header_row,
              ", not row 1."
            )
          )
        }

        if (first_nonempty_row > 1) {
          issues <- c(
            issues,
            paste0("First non-empty row appears to be row ", first_nonempty_row, ".")
          )
        }

        preview_names <- unlist(preview[probable_header_row, , drop = FALSE], use.names = FALSE)
        preview_names <- trimws(as.character(preview_names))

        blank_columns <- which(is.na(preview_names) | preview_names == "")

        nonblank_names <- preview_names[!is.na(preview_names) & preview_names != ""]
        duplicate_column_names <- unique(nonblank_names[duplicated(nonblank_names)])

        if (length(blank_columns) > 0) {
          issues <- c(
            issues,
            paste0("Unnamed columns detected in probable header row: ", paste(blank_columns, collapse = ", "), ".")
          )
        }

        if (length(duplicate_column_names) > 0) {
          issues <- c(
            issues,
            paste0("Duplicate column names detected in probable header row: ", paste(duplicate_column_names, collapse = ", "), ".")
          )
        }
      }
    }

    if (isTRUE(check_styles) && ext == "xlsx") {
      styles_detected <- tryCatch(
        {
          wb <- openxlsx::loadWorkbook(path)
          length(openxlsx::getStyles(wb)) > 0
        },
        error = function(e) {
          notes <<- c(notes, paste0("Could not inspect workbook styles: ", conditionMessage(e)))
          NA
        }
      )

      if (isTRUE(styles_detected)) {
        notes <- c(
          notes,
          "Workbook contains custom formatting. If colors or styles encode meaning, verify manually."
        )
      }
    }
  }

  if (ext %in% c("csv", "tsv", "txt")) {
    delim <- switch(ext, csv = ",", tsv = "\t", txt = "\t")

    preview <- tryCatch(
      readr::read_delim(
        file = path,
        delim = delim,
        n_max = preview_rows,
        show_col_types = FALSE,
        name_repair = "minimal"
      ),
      error = function(e) {
        issues <<- c(issues, paste0("Could not preview delimited file: ", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(preview)) {
      preview_names <- names(preview)
      blank_columns <- which(is.na(preview_names) | preview_names == "")

      nonblank_names <- preview_names[!is.na(preview_names) & preview_names != ""]
      duplicate_column_names <- unique(nonblank_names[duplicated(nonblank_names)])

      if (length(blank_columns) > 0) {
        issues <- c(issues, paste0("Unnamed columns detected: ", paste(blank_columns, collapse = ", "), "."))
      }

      if (length(duplicate_column_names) > 0) {
        issues <- c(
          issues,
          paste0("Duplicate column names detected: ", paste(duplicate_column_names, collapse = ", "), ".")
        )
      }
    }
  }

  if (ext %in% c("sav", "dta", "sas7bdat", "xpt")) {
    notes <- c(
      notes,
      "Labelled data format detected. Variable and value labels will be preserved by ReadSciData()."
    )
  }

  out <- list(
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    file_name = basename(path),
    extension = ext,
    file_size_bytes = unname(file_info$size),
    file_size_mb = round(unname(file_info$size) / 1024^2, 3),
    modified_time = unname(file_info$mtime),
    sheets = sheets,
    selected_sheet = selected_sheet,
    probable_header_row = probable_header_row,
    first_nonempty_row = first_nonempty_row,
    blank_columns = blank_columns,
    duplicate_column_names = duplicate_column_names,
    styles_detected = styles_detected,
    preview = preview,
    preview_names = preview_names,
    issues = issues,
    notes = notes
  )

  if (!quiet) {
    .PrintSciDataInspection(out)
  }

  invisible(out)
}

.PrintSciDataInspection <- function(x) {
  cat("File inspection: ", x$file_name, "\n", sep = "")
  cat("Type: .", x$extension, " | Size: ", x$file_size_mb, " MB\n", sep = "")

  if (!is.null(x$sheets)) {
    cat("Sheets: ", paste(x$sheets, collapse = ", "), "\n", sep = "")
  }

  if (!is.null(x$selected_sheet)) {
    cat("Selected sheet: ", x$selected_sheet, "\n", sep = "")
  }

  if (!is.null(x$probable_header_row)) {
    cat("Header row: ", x$probable_header_row, "\n", sep = "")
  }

  if (length(x$blank_columns) > 0) {
    cat("Unnamed columns: ", paste(x$blank_columns, collapse = ", "), "\n", sep = "")
  }

  if (length(x$duplicate_column_names) > 0) {
    cat("Duplicate names: ", paste(x$duplicate_column_names, collapse = ", "), "\n", sep = "")
  }

  if (length(x$issues) == 0 && length(x$notes) == 0) {
    cat("No obvious import issues detected.\n")
    return(invisible(x))
  }

  if (length(x$issues) > 0) {
    cat("\nWarnings:\n")
    for (i in seq_along(x$issues)) {
      cat("  - ", x$issues[[i]], "\n", sep = "")
    }
  }

  if (length(x$notes) > 0) {
    cat("\nNotes:\n")
    for (i in seq_along(x$notes)) {
      cat("  - ", x$notes[[i]], "\n", sep = "")
    }
  }

  invisible(x)
}
