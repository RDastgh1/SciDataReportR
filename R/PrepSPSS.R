#' Prepare a data frame for SPSS export
#'
#' `PrepSPSS()` converts column names to SPSS-safe variable names while
#' preserving the original column names as variable labels. It also returns a
#' name map so users can track the original and exported names.
#'
#' This is useful before writing data to `.sav` using `haven::write_sav()`.
#'
#' @param data A data frame or tibble.
#' @param path Optional file path ending in `.sav`. If supplied, the prepared
#'   data are written using `haven::write_sav()`.
#' @param name_map_path Optional file path for saving the name map as a CSV.
#' @param return One of `"list"`, `"data"`, or `"map"`.
#' @param quiet Logical. If `FALSE`, prints a short summary of renamed columns.
#' @param max_length Maximum variable-name length. Defaults to 64.
#' @param compress Compression type passed to `haven::write_sav()`.
#' @param ... Additional arguments passed to `haven::write_sav()` if `path` is supplied.
#'
#' @return Depending on `return`, either a list, data frame, or name map.
#' @export
PrepSPSS <- function(
  data,
  path = NULL,
  name_map_path = NULL,
  return = c("list", "data", "map"),
  quiet = FALSE,
  max_length = 64,
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

  max_length <- as.integer(max_length)

  original_names <- names(data)

  if (is.null(original_names)) {
    stop("`data` must have column names.", call. = FALSE)
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

  make_unique_with_limit <- function(x, max_length = 64) {
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

  spss_names <- make_unique_with_limit(spss_names, max_length = max_length)

  name_map <- tibble::tibble(
    original_name = original_names,
    spss_name = spss_names,
    changed = original_names != spss_names
  )

  data_spss <- data
  names(data_spss) <- spss_names

  for (i in seq_along(data_spss)) {
    attr(data_spss[[i]], "label") <- original_names[[i]]
  }

  if (!quiet) {
    n_changed <- sum(name_map$changed)

    cat("SPSS preparation\n")
    cat("----------------\n")
    cat("Columns checked: ", nrow(name_map), "\n", sep = "")
    cat("Columns renamed: ", n_changed, "\n", sep = "")

    if (n_changed > 0) {
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

    cat("\nOriginal names were saved as variable labels.\n")
  }

  if (!is.null(name_map_path)) {
    readr::write_csv(name_map, name_map_path)
  }

  if (!is.null(path)) {
    haven::write_sav(
      data = data_spss,
      path = path,
      compress = compress,
      ...
    )
  }

  out <- list(
    data = data_spss,
    name_map = name_map
  )

  if (return == "data") {
    return(data_spss)
  }

  if (return == "map") {
    return(name_map)
  }

  out
}
