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
                            star_nudge_y = 0.15,
                            remove_existing_stars = TRUE,
                            p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            p_labels = c("***","**","*","")) {

  star_from <- match.arg(star_from)

  # Check structure of res
  has_fdr <- !is.null(res$FDRCorrected) && !is.null(res$FDRCorrected$plot) && inherits(res$FDRCorrected$plot, "ggplot")
  has_raw <- !is.null(res$Unadjusted) && !is.null(res$Unadjusted$plot) && inherits(res$Unadjusted$plot, "ggplot")

  # Helper functions
  is_text_layer <- function(layer) {
    inherits(layer$geom, "GeomText") || inherits(layer$geom, "GeomLabel")
  }

  # Function to detect if a plot has star mapping
  plot_has_star_mapping <- function(p) {
    if (!inherits(p, "ggplot")) return(FALSE)

    for (layer in p$layers) {
      if (is_text_layer(layer)) {
        # Check if the layer has a label aesthetic
        if (!is.null(layer$mapping$label)) {
          label_expr <- deparse(layer$mapping$label)
          if (grepl("star", label_expr, ignore.case = TRUE)) {
            return(TRUE)
          }
        }
      }
    }
    return(FALSE)
  }

  # Select the appropriate plot based on star_from
  pick_plot <- function() {
    if (star_from %in% c("fdr", "P_adj")) {
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) {
        warning("Requested FDR stars but only unadjusted plot available")
        return(res$Unadjusted$plot)
      }
    } else if (star_from %in% c("raw", "P")) {
      if (has_raw) return(res$Unadjusted$plot)
      if (has_fdr) {
        warning("Requested raw stars but only FDR-corrected plot available")
        return(res$FDRCorrected$plot)
      }
    } else if (star_from == "existing") {
      # Try to use whichever plot has stars mapped
      if (has_fdr && plot_has_star_mapping(res$FDRCorrected$plot)) return(res$FDRCorrected$plot)
      if (has_raw && plot_has_star_mapping(res$Unadjusted$plot)) return(res$Unadjusted$plot)
      # Default to FDR if available
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from == "column") {
      if (is.null(star_col)) stop("star_from='column' requires a `star_col`.")
      if (has_fdr && star_col %in% names(res$FDRCorrected$plot$data)) return(res$FDRCorrected$plot)
      if (has_raw && star_col %in% names(res$Unadjusted$plot$data)) return(res$Unadjusted$plot)
      stop(sprintf("Column '%s' not found in either plot's data.", star_col))
    }

    stop("No valid ggplot found in res. Check that res is output from PlotCorrelationsHeatmap().")
  }

  p <- pick_plot()
  d <- p$data

  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # Function to compute stars from p-values
  stars_from_p <- function(pvec) {
    as.character(cut(pvec, breaks = p_breaks, labels = p_labels))
  }

  # Detect which star column was originally mapped (if any)
  existing_star_col <- NA_character_
  for (layer in p$layers) {
    if (is_text_layer(layer) && !is.null(layer$mapping$label)) {
      label_expr <- deparse(layer$mapping$label)
      # Check if it's mapping to stars or stars_FDR
      if (grepl("^stars_FDR$", label_expr)) {
        existing_star_col <- "stars_FDR"
        break
      } else if (grepl("^stars$", label_expr)) {
        existing_star_col <- "stars"
        break
      }
    }
  }

  # Select stars based on star_from
  pick_stars <- function() {
    cols <- names(d)

    if (star_from == "existing") {
      # Use whatever was originally mapped, or default to appropriate stars
      if (!is.na(existing_star_col) && existing_star_col %in% cols) {
        return(as.character(d[[existing_star_col]]))
      }
      # If FDR plot, default to stars_FDR
      if (has_fdr && identical(p, res$FDRCorrected$plot) && "stars_FDR" %in% cols) {
        return(as.character(d$stars_FDR))
      }
      # Otherwise use regular stars
      if ("stars" %in% cols) return(as.character(d$stars))
      # Compute from p-values if no star columns
      if ("P_adj" %in% cols) return(stars_from_p(d$P_adj))
      if ("P" %in% cols) return(stars_from_p(d$P))
      return(rep("", nrow(d)))
    }

    if (star_from == "raw") {
      if ("stars" %in% cols) return(as.character(d$stars))
      if ("P" %in% cols) return(stars_from_p(d$P))
      warning("No raw p-values or stars found")
      return(rep("", nrow(d)))
    }

    if (star_from == "fdr") {
      if ("stars_FDR" %in% cols) return(as.character(d$stars_FDR))
      if ("P_adj" %in% cols) return(stars_from_p(d$P_adj))
      warning("No FDR-adjusted p-values or stars found")
      return(rep("", nrow(d)))
    }

    if (star_from == "P") {
      if ("P" %in% cols) return(stars_from_p(d$P))
      stop("Column 'P' not found in plot data")
    }

    if (star_from == "P_adj") {
      if ("P_adj" %in% cols) return(stars_from_p(d$P_adj))
      stop("Column 'P_adj' not found in plot data")
    }

    if (star_from == "column") {
      if (is.null(star_col) || !star_col %in% cols) {
        stop(sprintf("star_col '%s' not found in plot data", star_col))
      }
      return(as.character(d[[star_col]]))
    }

    rep("", nrow(d))
  }

  # Create labels
  d$label_r <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- pick_stars()

  # Remove existing text layers if requested
  if (remove_existing_stars && length(p$layers) > 0) {
    keep <- vapply(p$layers, function(layer) {
      if (!is_text_layer(layer)) return(TRUE)

      # Check if this layer is mapping to stars
      if (!is.null(layer$mapping$label)) {
        label_expr <- deparse(layer$mapping$label)
        if (grepl("star", label_expr, ignore.case = TRUE)) {
          return(FALSE)
        }
      }

      # Check for static star labels
      if (!is.null(layer$aes_params$label)) {
        if (all(grepl("^\\*+$|^$|^ns$", layer$aes_params$label))) {
          return(FALSE)
        }
      }

      return(TRUE)
    }, logical(1))

    p$layers <- p$layers[keep]
  }

  # Add new layers
  p +
    ggplot2::geom_text(
      data = d,
      mapping = ggplot2::aes(label = label_r),
      inherit.aes = TRUE,
      size = r_size,
      color = r_color,
      position = ggplot2::position_nudge(y = r_nudge_y),
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = d,
      mapping = ggplot2::aes(label = label_star),
      inherit.aes = TRUE,
      size = star_size,
      color = star_color,
      position = ggplot2::position_nudge(y = star_nudge_y),
      na.rm = TRUE
    )
}
