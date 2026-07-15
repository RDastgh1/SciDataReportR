#' Prepare a data frame for SPSS export
#'
#' `PrepSPSS()` converts column names to SPSS-safe variable names while
#' preserving the original column names as variable labels. It also returns and
#' optionally writes a name map so users can track original and exported names.
#'
#' SPSS caps value labels (factor levels and `haven::labelled()` value labels)
#' at 120 bytes and variable labels at 256 bytes; `haven::write_sav()` refuses
#' to write files that exceed these limits. `PrepSPSS()` truncates over-long
#' labels to fit (marking them with `"..."`), de-duplicates factor levels that
#' collide after truncation, and records every truncation in a label map so no
#' information is silently lost.
#'
#' @param data A data frame or tibble.
#' @param path Optional file path ending in `.sav`. If supplied, the prepared
#'   data are written using `haven::write_sav()`.
#' @param name_map_path Optional file path for saving the name map as a CSV.
#' @param label_map_path Optional file path for saving the label map (all
#'   truncated labels alongside their original text) as a CSV.
#' @param return One of `"list"`, `"data"`, or `"map"`.
#' @param quiet Logical. If `FALSE`, prints a compact summary.
#' @param show_map Logical. If `TRUE`, prints the full original-to-SPSS name map.
#'   The default is `FALSE` because large scientific datasets can have thousands
#'   of renamed variables.
#' @param max_length Maximum SPSS variable-name length. Defaults to 64.
#' @param max_label_length Maximum value-label length in bytes. Defaults to 120,
#'   the SPSS limit enforced by `haven::write_sav()`. Factor levels and value
#'   labels longer than this are truncated and recorded in the label map.
#' @param compress Compression type passed to `haven::write_sav()`.
#' @param ... Additional arguments passed to `haven::write_sav()` if `path` is
#'   supplied.
#'
#' @return Depending on `return`, either a list (with elements `data`,
#'   `name_map`, and `label_map`), the prepared data frame, or the name map.
#' @export
PrepSPSS <- function(
  data,
  path = NULL,
  name_map_path = NULL,
  label_map_path = NULL,
  return = c("list", "data", "map"),
  quiet = FALSE,
  show_map = FALSE,
  max_length = 64,
  max_label_length = 120,
  compress = "byte",
  ...
) {
  return <- match.arg(return)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (ncol(data) == 0) {
    stop("`data` must have at least one column.", call. = FALSE)
  }

  if (!is.numeric(max_length) || length(max_length) != 1 || max_length < 8) {
    stop("`max_length` must be a single number greater than or equal to 8.", call. = FALSE)
  }

  if (!is.numeric(max_label_length) || length(max_label_length) != 1 || max_label_length < 8) {
    stop("`max_label_length` must be a single number greater than or equal to 8.", call. = FALSE)
  }

  max_length <- as.integer(max_length)
  max_label_length <- as.integer(max_label_length)

  original_names <- names(data)

  if (is.null(original_names)) {
    stop("`data` must have column names.", call. = FALSE)
  }

  if (!requireNamespace("janitor", quietly = TRUE)) {
    stop(
      "Package `janitor` is required to clean variable names for SPSS. ",
      "Install it with install.packages('janitor').",
      call. = FALSE
    )
  }

  spss_names <- janitor::make_clean_names(
    original_names,
    case = "snake",
    ascii = TRUE,
    use_make_names = FALSE,
    allow_dupes = FALSE
  )

  spss_names <- gsub("[^a-z0-9_]", "_", spss_names)
  spss_names <- gsub("_+", "_", spss_names)
  spss_names <- gsub("^_|_$", "", spss_names)
  spss_names <- ifelse(spss_names == "", "var", spss_names)
  spss_names <- ifelse(grepl("^[0-9]", spss_names), paste0("v", spss_names), spss_names)
  spss_names <- .MakeUniqueWithLimit(spss_names, max_length = max_length)

  name_map <- tibble::tibble(
    original_name = original_names,
    spss_name = spss_names,
    changed = original_names != spss_names
  )

  data_spss <- data
  names(data_spss) <- spss_names

  label_records <- list()
  max_variable_label_bytes <- 256L

  for (i in seq_along(data_spss)) {
    # Variable label: the original column name, capped at the SPSS 256-byte
    # limit for variable labels.
    variable_label <- original_names[[i]]

    if (nchar(variable_label, type = "bytes") > max_variable_label_bytes) {
      truncated_label <- .TruncateLabel(variable_label, max_variable_label_bytes)

      label_records[[length(label_records) + 1]] <- tibble::tibble(
        spss_name = spss_names[[i]],
        label_type = "variable_label",
        value = NA_character_,
        original_label = variable_label,
        spss_label = truncated_label
      )

      variable_label <- truncated_label
    }

    attr(data_spss[[i]], "label") <- variable_label

    # Value labels from factor levels, capped at the SPSS 120-byte limit.
    if (is.factor(data_spss[[i]])) {
      old_levels <- levels(data_spss[[i]])
      too_long <- nchar(old_levels, type = "bytes") > max_label_length

      if (any(too_long)) {
        new_levels <- old_levels
        new_levels[too_long] <- .TruncateLabel(old_levels[too_long], max_label_length)
        new_levels <- .DedupeLabels(new_levels, changed = too_long, max_bytes = max_label_length)

        label_records[[length(label_records) + 1]] <- tibble::tibble(
          spss_name = spss_names[[i]],
          label_type = "value_label",
          value = as.character(which(too_long)),
          original_label = old_levels[too_long],
          spss_label = new_levels[too_long]
        )

        levels(data_spss[[i]]) <- new_levels
      }
    }

    # Value labels from haven-style labelled vectors, same 120-byte limit.
    value_labels <- attr(data_spss[[i]], "labels", exact = TRUE)

    if (!is.null(value_labels) && !is.null(names(value_labels))) {
      old_labels <- names(value_labels)
      too_long <- nchar(old_labels, type = "bytes") > max_label_length

      if (any(too_long)) {
        new_labels <- old_labels
        new_labels[too_long] <- .TruncateLabel(old_labels[too_long], max_label_length)
        new_labels <- .DedupeLabels(new_labels, changed = too_long, max_bytes = max_label_length)

        label_records[[length(label_records) + 1]] <- tibble::tibble(
          spss_name = spss_names[[i]],
          label_type = "value_label",
          value = as.character(unname(value_labels[too_long])),
          original_label = old_labels[too_long],
          spss_label = new_labels[too_long]
        )

        names(value_labels) <- new_labels
        attr(data_spss[[i]], "labels") <- value_labels
      }
    }
  }

  if (length(label_records) > 0) {
    label_map <- do.call(rbind, label_records)
  } else {
    label_map <- tibble::tibble(
      spss_name = character(),
      label_type = character(),
      value = character(),
      original_label = character(),
      spss_label = character()
    )
  }

  if (!is.null(name_map_path)) {
    readr::write_csv(name_map, name_map_path)
  }

  if (!is.null(label_map_path)) {
    readr::write_csv(label_map, label_map_path)
  }

  if (!is.null(path)) {
    haven::write_sav(
      data = data_spss,
      path = path,
      compress = compress,
      ...
    )
  }

  if (!quiet) {
    n_changed <- sum(name_map$changed)
    n_truncated <- nrow(label_map)

    cat("SPSS preparation\n")
    cat("----------------\n")
    cat("Columns checked: ", nrow(name_map), "\n", sep = "")
    cat("Columns renamed: ", n_changed, "\n", sep = "")
    cat("Original names saved as variable labels.\n")

    if (n_truncated > 0) {
      cat(
        "Labels truncated to fit SPSS limits: ", n_truncated,
        " (full text kept in the label map).\n",
        sep = ""
      )
    }

    if (!is.null(path)) {
      cat("SPSS file written: ", path, "\n", sep = "")
    }

    if (!is.null(name_map_path)) {
      cat("Name map written: ", name_map_path, "\n", sep = "")
    }

    if (!is.null(label_map_path)) {
      cat("Label map written: ", label_map_path, "\n", sep = "")
    }

    if (isTRUE(show_map) && n_changed > 0) {
      cat("\nRenamed columns:\n")
      changed_rows <- name_map[name_map$changed, , drop = FALSE]

      for (i in seq_len(nrow(changed_rows))) {
        cat(
          i,
          ". ",
          changed_rows$original_name[[i]],
          " -> ",
          changed_rows$spss_name[[i]],
          "\n",
          sep = ""
        )
      }
    }
  }

  out <- list(
    data = data_spss,
    name_map = name_map,
    label_map = label_map
  )

  if (return == "data") {
    return(data_spss)
  }

  if (return == "map") {
    return(name_map)
  }

  out
}

.MakeUniqueWithLimit <- function(x, max_length = 64) {
  out <- character(length(x))
  used <- character()

  for (i in seq_along(x)) {
    base <- substr(x[[i]], 1, max_length)

    if (!base %in% used) {
      out[[i]] <- base
      used <- c(used, base)
    } else {
      counter <- 1

      repeat {
        suffix <- paste0("_", counter)
        base_length <- max_length - nchar(suffix)
        candidate <- paste0(substr(base, 1, base_length), suffix)

        if (!candidate %in% used) {
          out[[i]] <- candidate
          used <- c(used, candidate)
          break
        }

        counter <- counter + 1
      }
    }
  }

  out
}

# Cut a single string so its UTF-8 representation fits within `max_bytes`,
# without splitting a multi-byte character.
.CutBytes <- function(x, max_bytes) {
  n <- nchar(x)

  while (n > 0 && nchar(substr(x, 1, n), type = "bytes") > max_bytes) {
    n <- n - 1
  }

  substr(x, 1, n)
}

# Truncate labels to at most `max_bytes` bytes, marking truncated labels with
# `marker`. Labels already within the limit are returned unchanged.
.TruncateLabel <- function(x, max_bytes, marker = "...") {
  marker_bytes <- nchar(marker, type = "bytes")

  vapply(
    x,
    function(label) {
      if (is.na(label) || nchar(label, type = "bytes") <= max_bytes) {
        return(label)
      }

      paste0(.CutBytes(label, max_bytes - marker_bytes), marker)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# Resolve collisions created by truncation. Only entries flagged in `changed`
# may be renamed further; untouched labels always keep their original text.
.DedupeLabels <- function(labels, changed, max_bytes) {
  used <- labels[!changed]

  for (i in which(changed)) {
    candidate <- labels[[i]]
    counter <- 1

    while (candidate %in% used) {
      suffix <- paste0("_", counter)
      candidate <- paste0(
        .CutBytes(labels[[i]], max_bytes - nchar(suffix, type = "bytes")),
        suffix
      )
      counter <- counter + 1
    }

    labels[[i]] <- candidate
    used <- c(used, candidate)
  }

  labels
}
