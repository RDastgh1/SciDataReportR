#' Inspect a scientific data file before import
#'
#' `InspectFile()` checks for common file import problems, especially in Excel
#' workbooks where metadata rows, multiple sheets, and cell formatting can affect
#' interpretation.
#'
#' @param path Path to the file.
#' @param sheet Sheet name or index for Excel files. If `NULL`, the first sheet is inspected.
#' @param preview_rows Number of rows to preview when checking for metadata/header issues.
#' @param check_styles Logical. For Excel files, check whether custom styles/formatting exist.
#' @param check_sheets Logical. For Excel files, check for multiple sheets.
#' @param check_header Logical. For Excel files, attempt to detect whether the header is not row 1.
#' @param quiet Logical. If `FALSE`, print the inspection result.
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

  issues <- character()
  notes <- character()
  sheets <- NULL
  selected_sheet <- NULL
  probable_header_row <- NULL
  preview <- NULL
  styles_detected <- NA
  recommendation <- NULL

  supported_extensions <- c(
    "csv", "tsv", "txt",
    "xlsx", "xls",
    "rds", "rda", "rdata",
    "sav", "dta", "sas7bdat", "xpt",
    "parquet", "feather",
    "json"
  )

  if (!ext %in% supported_extensions) {
    issues <- c(
      issues,
      paste0("Unsupported or unrecognized file extension: .", ext)
    )
  }

  if (ext %in% c("xlsx", "xls")) {
    sheets <- tryCatch(
      readxl::excel_sheets(path),
      error = function(e) {
        issues <<- c(
          issues,
          paste0("Could not read Excel sheet names: ", conditionMessage(e))
        )
        NULL
      }
    )

    if (!is.null(sheets) && length(sheets) > 0) {
      if (is.null(sheet)) {
        selected_sheet <- sheets[[1]]
      } else if (is.numeric(sheet)) {
        if (sheet < 1 || sheet > length(sheets)) {
          issues <- c(
            issues,
            paste0(
              "`sheet = ", sheet,
              "` is outside the available sheet range 1-",
              length(sheets), "."
            )
          )
          selected_sheet <- sheets[[1]]
        } else {
          selected_sheet <- sheets[[sheet]]
        }
      } else {
        selected_sheet <- as.character(sheet)

        if (!selected_sheet %in% sheets) {
          issues <- c(
            issues,
            paste0(
              "Requested sheet `", selected_sheet,
              "` was not found in the workbook."
            )
          )
        }
      }

      if (check_sheets && length(sheets) > 1) {
        issues <- c(
          issues,
          paste0(
            "Workbook has multiple sheets: ",
            paste(sheets, collapse = ", "),
            ". Specify `sheet =` to avoid importing the wrong sheet."
          )
        )
      }
    }

    if (!is.null(selected_sheet) && selected_sheet %in% sheets) {
      preview <- tryCatch(
        readxl::read_excel(
          path,
          sheet = selected_sheet,
          col_names = FALSE,
          n_max = preview_rows,
          col_types = "text",
          .name_repair = "minimal"
        ),
        error = function(e) {
          issues <<- c(
            issues,
            paste0("Could not preview Excel sheet: ", conditionMessage(e))
          )
          NULL
        }
      )
    }

    if (check_header && !is.null(preview) && nrow(preview) > 0) {
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
              ", not row 1. There may be metadata above the column names."
            )
          )
        }

        if (first_nonempty_row > 1) {
          issues <- c(
            issues,
            paste0(
              "The first non-empty row appears to be row ",
              first_nonempty_row,
              ". The sheet may contain leading blank rows."
            )
          )
        }
      }
    }

    if (check_styles && ext == "xlsx") {
      styles_detected <- tryCatch(
        {
          wb <- openxlsx::loadWorkbook(path)
          length(openxlsx::getStyles(wb)) > 0
        },
        error = function(e) {
          notes <<- c(
            notes,
            paste0("Could not inspect workbook styles: ", conditionMessage(e))
          )
          NA
        }
      )

      if (isTRUE(styles_detected)) {
        issues <- c(
          issues,
          "Workbook contains custom styles or formatting. Color coding may contain information that will not be captured by a standard data import."
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

  if (length(issues) == 0) {
    recommendation <- "No obvious import issues detected."
  } else {
    recommendation <- "Review the issues before trusting an automated import."
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
    styles_detected = styles_detected,
    preview = preview,
    issues = issues,
    notes = notes,
    recommendation = recommendation
  )

  if (!quiet) {
    cat("File inspection\n")
    cat("---------------\n")
    cat("File: ", out$file_name, "\n", sep = "")
    cat("Type: .", out$extension, "\n", sep = "")
    cat("Size: ", out$file_size_mb, " MB\n", sep = "")

    if (!is.null(out$sheets)) {
      cat("Sheets: ", paste(out$sheets, collapse = ", "), "\n", sep = "")
    }

    if (!is.null(out$selected_sheet)) {
      cat("Selected sheet: ", out$selected_sheet, "\n", sep = "")
    }

    if (!is.null(out$probable_header_row)) {
      cat("Probable header row: ", out$probable_header_row, "\n", sep = "")
    }

    if (!is.na(out$styles_detected)) {
      cat(
        "Styles detected: ",
        ifelse(isTRUE(out$styles_detected), "yes", "no"),
        "\n",
        sep = ""
      )
    }

    cat("\n")

    if (length(out$issues) == 0) {
      cat("No obvious import issues detected.\n")
    } else {
      cat("Potential issues:\n")

      for (i in seq_along(out$issues)) {
        cat(i, ". ", out$issues[[i]], "\n", sep = "")
      }
    }

    if (length(out$notes) > 0) {
      cat("\nNotes:\n")

      for (i in seq_along(out$notes)) {
        cat(i, ". ", out$notes[[i]], "\n", sep = "")
      }
    }

    cat("\nRecommendation: ", out$recommendation, "\n", sep = "")
  }

  invisible(out)
}
