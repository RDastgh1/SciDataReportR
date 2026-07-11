#' Derive Freesurfer bilateral totals and ICV-adjusted measures
#'
#' Automatically derives bilateral Freesurfer volume measures from ASEG and DKT
#' outputs, then creates intracranial-volume-adjusted ratios.
#'
#' This function is designed for Freesurfer data frames that contain ASEG,
#' DKT, and global Freesurfer volume variables. It supports both native
#' Freesurfer-style ASEG names using hyphens and cleaned names using
#' underscores.
#'
#' Supported ASEG-style bilateral pairs include names such as:
#'
#' * `Left-Hippocampus`
#' * `Right-Hippocampus`
#' * `Left_Caudate`
#' * `Right_Caudate`
#'
#' Supported DKT-style bilateral cortical pairs include names such as:
#'
#' * `lh_fusiform_volume`
#' * `rh_fusiform_volume`
#'
#' Supported intracranial volume columns are:
#'
#' * `EstimatedTotalIntraCranialVol`
#' * `eTIV`
#'
#' If both `EstimatedTotalIntraCranialVol` and `eTIV` are present, the function
#' checks that they are equivalent before deriving ICV-adjusted variables. If
#' the two columns differ, the function stops and asks the user to resolve which
#' intracranial volume variable should be used.
#'
#' Bilateral totals are computed as:
#'
#' `left + right`
#'
#' ICV-adjusted ratios are computed as:
#'
#' `volume / intracranial volume`
#'
#' ICV-adjusted ratios are produced for each individual left/right (or
#' `lh`/`rh`) measure, for the derived bilateral totals, and for the global
#' Freesurfer volumes. For example, `Left-Hippocampus` and `Right-Hippocampus`
#' yield `Left_Hippocampus_icv`, `Right_Hippocampus_icv`, and
#' `Hippocampus_total_icv`.
#'
#' The function returns only the newly derived variables, not the original data.
#' This makes it convenient to append the derived columns with `cbind()` or
#' `dplyr::bind_cols()`.
#'
#' @param data A data frame containing Freesurfer ASEG and/or DKT variables.
#' @param verbose Logical. If `TRUE`, prints a short summary of derived variables.
#'   Default is `TRUE`.
#'
#' @return A data frame containing only newly derived variables, with the same
#'   number of rows as `data`. A derivation log is stored in the attribute
#'   `"Freesurfer_derivation_log"`.
#'
#' @examples
#' \dontrun{
#' fs_derived <- DeriveFreesurferVolumes(df_freesurfer)
#'
#' df_freesurfer <- dplyr::bind_cols(
#'   df_freesurfer,
#'   DeriveFreesurferVolumes(df_freesurfer)
#' )
#'
#' attr(fs_derived, "Freesurfer_derivation_log")
#' }
#'
#' @export
DeriveFreesurferVolumes <- function(
  data,
  verbose = TRUE
) {

  ############################################################
  ## Input checks
  ############################################################

  if (!is.data.frame(data)) {
    stop(
      "`data` must be a data frame.",
      call. = FALSE
    )
  }

  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop(
      "`verbose` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  ############################################################
  ## Internal helper functions
  ############################################################

  clean_output_name <- function(x) {
    x <- gsub("-", "_", x)
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  get_aseg_side <- function(x) {
    ifelse(
      grepl("^Left[-_]", x),
      "Left",
      ifelse(
        grepl("^Right[-_]", x),
        "Right",
        NA_character_
      )
    )
  }

  get_aseg_structure <- function(x) {
    structure <- sub("^Left[-_]", "", x)
    structure <- sub("^Right[-_]", "", structure)
    structure
  }

  ############################################################
  ## Determine intracranial volume variable
  ############################################################

  has_estimated_icv <- "EstimatedTotalIntraCranialVol" %in% names(data)
  has_etiv <- "eTIV" %in% names(data)

  if (!has_estimated_icv && !has_etiv) {
    stop(
      paste(
        "Neither `EstimatedTotalIntraCranialVol` nor `eTIV` was found.",
        "One of these columns is required for ICV-adjusted variables."
      ),
      call. = FALSE
    )
  }

  if (has_estimated_icv && !is.numeric(data$EstimatedTotalIntraCranialVol)) {
    stop(
      "`EstimatedTotalIntraCranialVol` must be numeric.",
      call. = FALSE
    )
  }

  if (has_etiv && !is.numeric(data$eTIV)) {
    stop(
      "`eTIV` must be numeric.",
      call. = FALSE
    )
  }

  if (has_estimated_icv && has_etiv) {

    same_icv <- isTRUE(
      all.equal(
        data$EstimatedTotalIntraCranialVol,
        data$eTIV,
        check.attributes = FALSE
      )
    )

    if (!same_icv) {
      stop(
        paste(
          "`EstimatedTotalIntraCranialVol` and `eTIV` differ.",
          "These should represent the same Freesurfer estimated intracranial volume.",
          "Please inspect both columns and decide which one should be used."
        ),
        call. = FALSE
      )
    }

    icv <- data$EstimatedTotalIntraCranialVol
    icv_source <- "EstimatedTotalIntraCranialVol"

  } else if (has_estimated_icv) {

    icv <- data$EstimatedTotalIntraCranialVol
    icv_source <- "EstimatedTotalIntraCranialVol"

  } else {

    icv <- data$eTIV
    icv_source <- "eTIV"

  }

  if (any(is.na(icv))) {
    warning(
      "The selected ICV variable contains missing values. Derived `_icv` variables will be missing for those rows.",
      call. = FALSE
    )
  }

  if (any(icv == 0, na.rm = TRUE)) {
    stop(
      "The selected ICV variable contains zero values. Cannot divide by zero.",
      call. = FALSE
    )
  }

  ############################################################
  ## Initialize output and tracking objects
  ############################################################

  out <- data.frame(row.names = seq_len(nrow(data)))

  aseg_total_cols <- character()
  dkt_total_cols <- character()
  global_icv_cols <- character()
  bilateral_icv_cols <- character()
  lateralized_icv_cols <- character()

  unmatched_aseg_left <- character()
  unmatched_aseg_right <- character()
  unmatched_dkt_left <- character()
  unmatched_dkt_right <- character()

  skipped_non_numeric_aseg <- character()
  skipped_non_numeric_dkt <- character()
  skipped_non_numeric_global <- character()

  ############################################################
  ## ASEG bilateral totals
  ##
  ## Supports:
  ## Left-Hippocampus + Right-Hippocampus = Hippocampus_total
  ## Left_Hippocampus + Right_Hippocampus = Hippocampus_total
  ############################################################

  aseg_cols <- names(data)[grepl("^(Left|Right)[-_]", names(data))]

  if (length(aseg_cols) > 0) {

    aseg_map <- data.frame(
      original_col = aseg_cols,
      side = get_aseg_side(aseg_cols),
      structure_raw = get_aseg_structure(aseg_cols),
      stringsAsFactors = FALSE
    )

    aseg_map$structure_key <- clean_output_name(aseg_map$structure_raw)

    left_map <- aseg_map[aseg_map$side == "Left", , drop = FALSE]
    right_map <- aseg_map[aseg_map$side == "Right", , drop = FALSE]

    for (i in seq_len(nrow(left_map))) {

      left_col <- left_map$original_col[i]
      structure_key <- left_map$structure_key[i]

      matching_right <- right_map$original_col[
        right_map$structure_key == structure_key
      ]

      if (length(matching_right) == 0) {
        unmatched_aseg_left <- c(
          unmatched_aseg_left,
          left_col
        )
        next
      }

      if (length(matching_right) > 1) {
        stop(
          paste(
            "Multiple right-side ASEG matches were found for",
            paste0("`", left_col, "`."),
            "Please inspect duplicated or inconsistently named Freesurfer columns."
          ),
          call. = FALSE
        )
      }

      right_col <- matching_right[1]

      if (!is.numeric(data[[left_col]]) || !is.numeric(data[[right_col]])) {
        skipped_non_numeric_aseg <- c(
          skipped_non_numeric_aseg,
          left_col,
          right_col
        )
        next
      }

      new_col <- paste0(
        structure_key,
        "_total"
      )

      out[[new_col]] <- data[[left_col]] + data[[right_col]]

      aseg_total_cols <- c(
        aseg_total_cols,
        new_col
      )
    }

    for (i in seq_len(nrow(right_map))) {

      right_col <- right_map$original_col[i]
      structure_key <- right_map$structure_key[i]

      matching_left <- left_map$original_col[
        left_map$structure_key == structure_key
      ]

      if (length(matching_left) == 0) {
        unmatched_aseg_right <- c(
          unmatched_aseg_right,
          right_col
        )
      }
    }
  }

  ############################################################
  ## DKT bilateral totals
  ##
  ## Example:
  ## lh_fusiform_volume + rh_fusiform_volume = fusiform_volume_total
  ############################################################

  lh_cols <- grep(
    "^lh_",
    names(data),
    value = TRUE
  )

  rh_cols <- grep(
    "^rh_",
    names(data),
    value = TRUE
  )

  for (lh_col in lh_cols) {

    region <- sub(
      "^lh_",
      "",
      lh_col
    )

    rh_col <- paste0(
      "rh_",
      region
    )

    if (!rh_col %in% names(data)) {
      unmatched_dkt_left <- c(
        unmatched_dkt_left,
        lh_col
      )
      next
    }

    if (!is.numeric(data[[lh_col]]) || !is.numeric(data[[rh_col]])) {
      skipped_non_numeric_dkt <- c(
        skipped_non_numeric_dkt,
        lh_col,
        rh_col
      )
      next
    }

    region_name <- sub(
      "_volume$",
      "",
      region
    )

    region_name <- clean_output_name(region_name)

    new_col <- paste0(
      region_name,
      "_volume_total"
    )

    out[[new_col]] <- data[[lh_col]] + data[[rh_col]]

    dkt_total_cols <- c(
      dkt_total_cols,
      new_col
    )
  }

  for (rh_col in rh_cols) {

    region <- sub(
      "^rh_",
      "",
      rh_col
    )

    lh_col <- paste0(
      "lh_",
      region
    )

    if (!lh_col %in% names(data)) {
      unmatched_dkt_right <- c(
        unmatched_dkt_right,
        rh_col
      )
    }
  }

  ############################################################
  ## Global Freesurfer ICV-adjusted variables
  ##
  ## These are not left/right totals, but are commonly adjusted by ICV.
  ############################################################

  global_volumes <- c(
    "BrainSegVol",
    "BrainSegVolNotVent_ASEG",
    "BrainSegVolNotVent_DKT",
    "CortexVol",
    "CerebralWhiteMatterVol",
    "SubCortGrayVol",
    "TotalGrayVol",
    "SupraTentorialVol",
    "SupraTentorialVolNotVent",
    "MaskVol"
  )

  global_volumes <- intersect(
    global_volumes,
    names(data)
  )

  for (global_col in global_volumes) {

    if (!is.numeric(data[[global_col]])) {
      skipped_non_numeric_global <- c(
        skipped_non_numeric_global,
        global_col
      )
      next
    }

    new_col <- paste0(
      global_col,
      "_icv"
    )

    out[[new_col]] <- data[[global_col]] / icv

    global_icv_cols <- c(
      global_icv_cols,
      new_col
    )
  }

  ############################################################
  ## ICV-adjusted lateralized (left/right) variables
  ##
  ## Each individual left/right ASEG column and lh/rh DKT column
  ## is divided by ICV, so ratios exist for the sides, not only
  ## for the bilateral totals.
  ##
  ## Formula:
  ## side_volume / ICV
  ############################################################

  lateralized_side_cols <- c(aseg_cols, lh_cols, rh_cols)

  for (side_col in lateralized_side_cols) {

    if (!is.numeric(data[[side_col]])) {
      next
    }

    new_col <- paste0(
      clean_output_name(side_col),
      "_icv"
    )

    out[[new_col]] <- data[[side_col]] / icv

    lateralized_icv_cols <- c(
      lateralized_icv_cols,
      new_col
    )
  }

  ############################################################
  ## ICV-adjusted bilateral total variables
  ##
  ## Formula:
  ## derived_total / ICV
  ############################################################

  total_cols <- c(
    aseg_total_cols,
    dkt_total_cols
  )

  for (total_col in total_cols) {

    if (!is.numeric(out[[total_col]])) {
      next
    }

    new_col <- paste0(
      total_col,
      "_icv"
    )

    out[[new_col]] <- out[[total_col]] / icv

    bilateral_icv_cols <- c(
      bilateral_icv_cols,
      new_col
    )
  }

  ############################################################
  ## Derivation log
  ############################################################

  derivation_log <- list(
    icv_source = icv_source,
    n_rows_input = nrow(data),
    n_rows_output = nrow(out),
    n_aseg_total_cols = length(aseg_total_cols),
    n_dkt_total_cols = length(dkt_total_cols),
    n_global_icv_cols = length(global_icv_cols),
    n_lateralized_icv_cols = length(lateralized_icv_cols),
    n_bilateral_icv_cols = length(bilateral_icv_cols),
    n_cols_output = ncol(out),
    aseg_total_cols = aseg_total_cols,
    dkt_total_cols = dkt_total_cols,
    global_icv_cols = global_icv_cols,
    lateralized_icv_cols = lateralized_icv_cols,
    bilateral_icv_cols = bilateral_icv_cols,
    unmatched_aseg_left = unmatched_aseg_left,
    unmatched_aseg_right = unmatched_aseg_right,
    unmatched_dkt_left = unmatched_dkt_left,
    unmatched_dkt_right = unmatched_dkt_right,
    skipped_non_numeric_aseg = unique(skipped_non_numeric_aseg),
    skipped_non_numeric_dkt = unique(skipped_non_numeric_dkt),
    skipped_non_numeric_global = unique(skipped_non_numeric_global)
  )

  attr(out, "Freesurfer_derivation_log") <- derivation_log

  ############################################################
  ## Optional summary message
  ############################################################

  if (verbose) {

    message("Derived Freesurfer measures:")
    message("  ICV source: ", icv_source)
    message("  ASEG bilateral totals: ", length(aseg_total_cols))
    message("  DKT bilateral totals: ", length(dkt_total_cols))
    message("  Global ICV-adjusted measures: ", length(global_icv_cols))
    message("  Lateralized (left/right) ICV-adjusted measures: ", length(lateralized_icv_cols))
    message("  Bilateral ICV-adjusted measures: ", length(bilateral_icv_cols))
    message("  Returned derived variables: ", ncol(out))

    n_unmatched <- length(unmatched_aseg_left) +
      length(unmatched_aseg_right) +
      length(unmatched_dkt_left) +
      length(unmatched_dkt_right)

    n_skipped_non_numeric <- length(unique(skipped_non_numeric_aseg)) +
      length(unique(skipped_non_numeric_dkt)) +
      length(unique(skipped_non_numeric_global))

    if (n_unmatched > 0) {
      message("  Unmatched left/right variables: ", n_unmatched)
      message("  Inspect `attr(output, 'Freesurfer_derivation_log')` for details.")
    }

    if (n_skipped_non_numeric > 0) {
      message("  Skipped non-numeric variables: ", n_skipped_non_numeric)
      message("  Inspect `attr(output, 'Freesurfer_derivation_log')` for details.")
    }
  }

  ############################################################
  ## Return only derived variables
  ############################################################

  out
}
