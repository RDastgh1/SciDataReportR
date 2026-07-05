#' Harmonize merge key columns before a safe merge
#'
#' Internal helper used by safe_merge().
#'
#' @return A list with elements `left` and `right` (the input data frames with
#'   harmonized key columns) and `report` (a data frame describing the
#'   harmonization applied to each key).
#' @keywords internal
HarmonizeMergeKeys <- function(
  df_before,
  df_add,
  by,
  key_parser = NULL,
  stop_on_failed_numeric = TRUE,
  n_examples = 5
) {
  if (!is.data.frame(df_before)) {
    stop("df_before must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(df_add)) {
    stop("df_add must be a data.frame.", call. = FALSE)
  }

  if (missing(by) || !is.character(by) || length(by) == 0) {
    stop("by must be a character vector of one or more merge keys.", call. = FALSE)
  }

  missing_before <- setdiff(by, names(df_before))
  missing_add <- setdiff(by, names(df_add))

  if (length(missing_before) > 0) {
    stop(
      "The following merge key(s) are missing from df_before: ",
      paste(missing_before, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_add) > 0) {
    stop(
      "The following merge key(s) are missing from df_add: ",
      paste(missing_add, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(key_parser) && is.null(names(key_parser))) {
    stop("key_parser must be a named list of functions.", call. = FALSE)
  }

  can_be_numeric <- function(x) {
    x_chr <- trimws(as.character(x))
    x_chr <- x_chr[!is.na(x_chr) & x_chr != ""]

    if (length(x_chr) == 0) {
      return(TRUE)
    }

    converted <- suppressWarnings(as.numeric(x_chr))
    all(!is.na(converted))
  }

  failed_numeric_values <- function(x, n = 5) {
    x_chr <- trimws(as.character(x))
    x_chr <- unique(x_chr[!is.na(x_chr) & x_chr != ""])
    converted <- suppressWarnings(as.numeric(x_chr))
    bad <- x_chr[is.na(converted)]

    head(bad, n)
  }

  apply_parser <- function(x, key) {
    if (!is.null(key_parser) && key %in% names(key_parser)) {
      parser_fun <- key_parser[[key]]

      if (!is.function(parser_fun)) {
        stop(
          "key_parser[['",
          key,
          "']] must be a function.",
          call. = FALSE
        )
      }

      return(parser_fun(x))
    }

    x
  }

  harmonized_left <- df_before
  harmonized_right <- df_add

  report <- tibble::tibble(
    Key = character(),
    LeftTypeBefore = character(),
    RightTypeBefore = character(),
    LeftTypeAfter = character(),
    RightTypeAfter = character(),
    ParserUsed = logical(),
    HarmonizedTo = character(),
    Action = character(),
    Status = character(),
    FailedLeftExamples = character(),
    FailedRightExamples = character()
  )

  for (key in by) {
    left_original <- harmonized_left[[key]]
    right_original <- harmonized_right[[key]]

    left_type_before <- paste(class(left_original), collapse = "/")
    right_type_before <- paste(class(right_original), collapse = "/")

    parser_used <- !is.null(key_parser) && key %in% names(key_parser)

    left_parsed <- apply_parser(left_original, key)
    right_parsed <- apply_parser(right_original, key)

    left_numeric <- can_be_numeric(left_parsed)
    right_numeric <- can_be_numeric(right_parsed)

    left_failed <- failed_numeric_values(left_parsed, n_examples)
    right_failed <- failed_numeric_values(right_parsed, n_examples)

    if (left_numeric && right_numeric) {
      left_harmonized <- suppressWarnings(as.numeric(trimws(as.character(left_parsed))))
      right_harmonized <- suppressWarnings(as.numeric(trimws(as.character(right_parsed))))

      harmonized_to <- "numeric"
      action <- dplyr::case_when(
        parser_used ~ "Applied key parser and converted both key columns to numeric.",
        left_type_before != right_type_before ~ "Converted both key columns to numeric.",
        TRUE ~ "Confirmed both key columns are numeric-compatible."
      )
      status <- "PASS"
    } else {
      if (isTRUE(stop_on_failed_numeric) &&
          (grepl("numeric|integer|double", left_type_before) ||
             grepl("numeric|integer|double", right_type_before))) {
        stop(
          "Merge key harmonization failed for key `",
          key,
          "`.\n\n",
          "One side is numeric-like, but the other side contains values that cannot be safely converted to numeric.\n\n",
          "Left type: ",
          left_type_before,
          "\n",
          "Right type: ",
          right_type_before,
          "\n\n",
          "Examples from df_before that cannot be numeric: ",
          ifelse(length(left_failed) == 0, "None", paste(left_failed, collapse = ", ")),
          "\n",
          "Examples from df_add that cannot be numeric: ",
          ifelse(length(right_failed) == 0, "None", paste(right_failed, collapse = ", ")),
          "\n\n",
          "Clean the key column before merging, or pass a key_parser to safe_merge().\n",
          "Example:\n",
          "safe_merge(..., key_parser = list(",
          key,
          " = function(x) stringr::str_remove(as.character(x), '^BRAIN_')))",
          call. = FALSE
        )
      }

      left_harmonized <- as.character(left_parsed)
      right_harmonized <- as.character(right_parsed)

      harmonized_to <- "character"
      action <- dplyr::case_when(
        parser_used ~ "Applied key parser and converted both key columns to character.",
        TRUE ~ "Converted both key columns to character."
      )
      status <- "WARNING"
    }

    harmonized_left[[key]] <- left_harmonized
    harmonized_right[[key]] <- right_harmonized

    report <- dplyr::bind_rows(
      report,
      tibble::tibble(
        Key = key,
        LeftTypeBefore = left_type_before,
        RightTypeBefore = right_type_before,
        LeftTypeAfter = paste(class(left_harmonized), collapse = "/"),
        RightTypeAfter = paste(class(right_harmonized), collapse = "/"),
        ParserUsed = parser_used,
        HarmonizedTo = harmonized_to,
        Action = action,
        Status = status,
        FailedLeftExamples = ifelse(length(left_failed) == 0, "", paste(left_failed, collapse = ", ")),
        FailedRightExamples = ifelse(length(right_failed) == 0, "", paste(right_failed, collapse = ", "))
      )
    )
  }

  list(
    left = harmonized_left,
    right = harmonized_right,
    report = report
  )
}
