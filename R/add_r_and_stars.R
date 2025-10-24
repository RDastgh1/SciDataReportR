#' Add r-values and significance stars to a correlations heatmap (pass `res`)
#' Works with the list returned by PlotCorrelationsHeatmap().
#'
#' @param res The list returned by PlotCorrelationsHeatmap()
#' @param star_from One of:
#'   "existing" (use whatever the chosen plot already mapped as stars),
#'   "raw" (use `stars` or compute from `P`),
#'   "fdr" (use `stars_FDR` or compute from `P_adj`),
#'   "P","P_adj" (compute from those columns),
#'   "column" (use `star_col`).
#' @param star_col Column name to use when star_from = "column"
#' @param r_var Column name for correlation values (default "R")
#' @param r_digits Decimal places for r labels
#' @param r_size,star_size Text sizes for r and stars
#' @param r_color,star_color Colors for r and stars
#' @param r_nudge_y,star_nudge_y Vertical nudges (r down, stars up)
#' @param remove_existing_stars If TRUE, remove any pre-existing star text layers
#' @param p_breaks,p_labels Cutpoints/labels for computing stars from p
#' @return A ggplot with r-values and stars added
#' @export
add_r_and_stars <- function(res,
                            star_from = c("existing","raw","fdr","P","P_adj","column"),
                            star_col = NULL,
                            r_var = "R",
                            r_digits = 2,
                            r_size = 3,
                            r_color = "black",
                            r_nudge_y = -0.28,
                            star_size = 6,
                            star_color = "black",
                            star_nudge_y = 0.28,
                            remove_existing_stars = TRUE,
                            p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            p_labels = c("***","**","*","")) {

  star_from <- match.arg(star_from)

  has_fdr <- is.list(res$FDRCorrected) && inherits(res$FDRCorrected$plot, "ggplot")
  has_raw <- is.list(res$Unadjusted)   && inherits(res$Unadjusted$plot,   "ggplot")

  # helpers ------------------------------------------------------------
  safe_as_label <- function(x) if (is.null(x)) NA_character_ else paste(deparse(x), collapse = "")
  is_text_layer <- function(layer) inherits(layer$geom, "GeomText") || inherits(layer$geom, "GeomLabel")

  plot_has_star_mapping <- function(p) {
    if (!inherits(p, "ggplot")) return(FALSE)
    any(vapply(p$layers, function(ly) {
      if (!is_text_layer(ly)) return(FALSE)
      lab <- safe_as_label(ly$mapping$label)
      # More inclusive check - look for "star" anywhere in the mapping
      !is.na(lab) && grepl("star", lab, ignore.case = TRUE)
    }, logical(1)))
  }

  pick_plot <- function() {
    if (star_from %in% c("fdr","P_adj")) {
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from %in% c("raw","P")) {
      if (has_raw) return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
    } else if (star_from == "existing") {
      if (has_fdr && plot_has_star_mapping(res$FDRCorrected$plot)) return(res$FDRCorrected$plot)
      if (has_raw && plot_has_star_mapping(res$Unadjusted$plot))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from == "column") {
      if (is.null(star_col)) stop("star_from='column' requires a `star_col`.")
      if (has_fdr && star_col %in% names(res$FDRCorrected$plot$data)) return(res$FDRCorrected$plot)
      if (has_raw && star_col %in% names(res$Unadjusted$plot$data))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    }
    stop("Couldn't find a valid ggplot in `res`.")
  }

  p <- pick_plot()
  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # detect original stars mapping (if any) -----------------------------------
  existing_star_col <- {
    labs <- vapply(p$layers, function(ly)
      if (is_text_layer(ly)) safe_as_label(ly$mapping$label) else NA_character_, character(1))
    labs <- labs[!is.na(labs)]
    # More inclusive pattern - look for anything with "star" in it
    labs <- labs[grepl("star", labs, ignore.case = TRUE)]
    labs[labs %in% names(d)][1]
  }

  stars_from_p <- function(pvec) as.character(cut(pvec, breaks = p_breaks, labels = p_labels))

  pick_stars <- function() {
    cols <- names(d)
    if (star_from == "existing") {
      if (!is.na(existing_star_col)) return(as.character(d[[existing_star_col]]))
      if ("stars_FDR" %in% cols) return(as.character(d$stars_FDR))
      if ("P_adj"    %in% cols)  return(stars_from_p(d$P_adj))
      if ("stars"    %in% cols)  return(as.character(d$stars))
      if ("P"        %in% cols)  return(stars_from_p(d$P))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "raw") {
      if ("stars" %in% cols) return(as.character(d$stars))
      if ("P"     %in% cols) return(stars_from_p(d$P))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "fdr") {
      if ("stars_FDR" %in% cols) return(as.character(d$stars_FDR))
      if ("P_adj"    %in% cols)  return(stars_from_p(d$P_adj))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "P")     return(stars_from_p(d$P))
    if (star_from == "P_adj") return(stars_from_p(d$P_adj))
    if (star_from == "column") {
      if (is.null(star_col) || !star_col %in% cols) {
        stop("star_from='column' requires a valid `star_col` in the plot data.")
      }
      return(as.character(d[[star_col]]))
    }
    rep(NA_character_, nrow(d))
  }

  # labels -------------------------------------------------------------------
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- pick_stars()

  # remove any pre-existing star text layers - MORE AGGRESSIVE REMOVAL ----
  if (remove_existing_stars && length(p$layers)) {
    keep <- vapply(p$layers, function(ly) {
      if (!is_text_layer(ly)) return(TRUE)

      # Get the label mapping
      lab <- safe_as_label(ly$mapping$label)
      if (is.na(lab)) return(TRUE)

      # Remove if it contains "star" anywhere OR if it's exactly "stars"
      contains_star <- grepl("star", lab, ignore.case = TRUE)

      # Also check if the layer has static text that might be stars
      # This catches cases where stars are added as static text
      static_label <- ly$aes_params$label
      if (!is.null(static_label) && length(static_label) == 1) {
        # Check if it's a star pattern
        is_star_pattern <- grepl("^\\*+$", static_label)
        if (is_star_pattern) return(FALSE)
      }

      # Remove if it maps to anything with "star" in the name
      return(!contains_star)
    }, logical(1))
    p$layers <- p$layers[keep]
  }

  # compose ------------------------------------------------------------------
  p +
    ggplot2::geom_text(
      data = d,
      mapping = ggplot2::aes(label = label_r),
      inherit.aes = TRUE,
      size = r_size, color = r_color,
      position = ggplot2::position_nudge(y = r_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = d,
      mapping = ggplot2::aes(label = label_star),
      inherit.aes = TRUE,
      size = star_size, color = star_color,
      position = ggplot2::position_nudge(y = star_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    )
}
