#' Add r-values and (optionally) a fresh set of significance stars to a heatmap
#' built by PlotCorrelationsHeatmap().
#'
#' @param res  List returned by PlotCorrelationsHeatmap()
#' @param star_from One of:
#'   "existing"  (keep the plot's current stars; only add r),
#'   "raw"       (use `stars` or compute from `P`),
#'   "fdr"       (use `stars_FDR` or compute from `P_adj`),
#'   "P","P_adj" (compute from those columns),
#'   "column"    (use `star_col`).
#' @param star_col Column name when star_from = "column"
#' @param r_var    Column for r-values in plot data (default "R")
#' @param r_digits Decimal places for r labels
#' @param r_size   Text size for r labels
#' @param r_color  Color for r labels
#' @param r_nudge_y Vertical nudge for r (negative = below center)
#' @param star_size Text size for NEW star labels (ignored if star_from="existing")
#' @param star_color Color for NEW star labels
#' @param star_nudge_y Vertical nudge for NEW star labels (ignored if "existing")
#' @param p_breaks Cutpoints for p→stars
#' @param p_labels Labels for p→stars
#' @return ggplot
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
                            p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            p_labels = c("***","**","*","")) {

  star_from <- match.arg(star_from)

  # helpers --------------------------------------------------------------
  .is_text <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")
  .label_name <- function(expr) {
    if (is.null(expr)) return(NA_character_)
    s <- paste(deparse(expr), collapse = "")
    s <- tolower(s)
    s <- sub("^after_stat\\((.*)\\)$", "\\1", s)  # after_stat(stars) -> stars
    s <- sub("^\\.?data\\$", "", s)               # .data$stars -> stars
    s <- gsub("`", "", s)                         # `stars` -> stars
    gsub("[^a-z0-9_]", "", s)                     # drop punctuation/spaces
  }
  .stars_from_p <- function(pvec, n) {
    if (!length(pvec)) return(rep("", n))
    as.character(cut(pvec, breaks = p_breaks, labels = p_labels))
  }
  .has_plot <- function(x) is.list(x) && inherits(x$plot, "ggplot")

  # choose base plot (intuitive defaulting) ------------------------------
  if (star_from %in% c("raw","P")) {
    if (!.has_plot(res$Unadjusted)) stop("Unadjusted plot not found in `res`.")
    p <- res$Unadjusted$plot
  } else if (star_from %in% c("fdr","P_adj")) {
    if (!.has_plot(res$FDRCorrected)) stop("FDRCorrected plot not found in `res`.")
    p <- res$FDRCorrected$plot
  } else if (star_from == "column") {
    if (is.null(star_col)) stop("star_from='column' requires `star_col`.")
    if (.has_plot(res$Unadjusted) && star_col %in% names(res$Unadjusted$plot$data)) {
      p <- res$Unadjusted$plot
    } else if (.has_plot(res$FDRCorrected) && star_col %in% names(res$FDRCorrected$plot$data)) {
      p <- res$FDRCorrected$plot
    } else stop("`star_col` not found in either plot's data.")
  } else { # "existing" → prefer the plot that actually mapped stars
    if (.has_plot(res$Unadjusted)) {
      pU <- res$Unadjusted$plot
      usedU <- any(vapply(pU$layers, \(ly) .is_text(ly) && .label_name(ly$mapping$label) %in% c("stars","starsfdr"), logical(1)))
      if (usedU) p <- pU
    }
    if (!exists("p") && .has_plot(res$FDRCorrected)) {
      pF <- res$FDRCorrected$plot
      usedF <- any(vapply(pF$layers, \(ly) .is_text(ly) && .label_name(ly$mapping$label) %in% c("stars","starsfdr"), logical(1)))
      if (usedF) p <- pF
    }
    if (!exists("p")) {
      # fall back to Unadjusted then FDR if neither mapped stars
      if (.has_plot(res$Unadjusted)) p <- res$Unadjusted$plot
      else if (.has_plot(res$FDRCorrected)) p <- res$FDRCorrected$plot
      else stop("Couldn't find a valid ggplot in `res`.")
    }
  }

  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # format r labels
  d$label_r <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])

  if (star_from == "existing") {
    # Keep the plot's original stars exactly as-is; only add r
    return(
      p +
        ggplot2::geom_text(
          data = d, ggplot2::aes(label = label_r),
          inherit.aes = TRUE, size = r_size, color = r_color,
          position = ggplot2::position_nudge(y = r_nudge_y),
          check_overlap = TRUE, na.rm = TRUE
        )
    )
  }

  # Otherwise: compute the star labels you asked for and replace any old star layers
  label_star <- switch(
    star_from,
    raw = {
      if      ("stars"    %in% names(d)) as.character(d$stars)
      else if ("P"        %in% names(d)) .stars_from_p(d$P, nrow(d))
      else stop("For RAW stars, need 'stars' or 'P' in plot data.")
    },
    fdr = {
      if      ("stars_FDR" %in% names(d)) as.character(d$stars_FDR)
      else if ("P_adj"     %in% names(d))  .stars_from_p(d$P_adj, nrow(d))
      else stop("For FDR stars, need 'stars_FDR' or 'P_adj' in plot data.")
    },
    P = {
      if (!"P" %in% names(d)) stop("Column 'P' not in plot data.")
      .stars_from_p(d$P, nrow(d))
    },
    P_adj = {
      if (!"P_adj" %in% names(d)) stop("Column 'P_adj' not in plot data.")
      .stars_from_p(d$P_adj, nrow(d))
    },
    column = {
      if (!star_col %in% names(d)) stop("`star_col` not found in plot data.")
      as.character(d[[star_col]])
    }
  )
  d$label_star <- label_star

  # Remove any prior star text layers (but leave non-text/title/etc.)
  if (length(p$layers)) {
    keep <- vapply(p$layers, function(ly) {
      if (!.is_text(ly)) return(TRUE)
      nm <- .label_name(ly$mapping$label)
      is.na(nm) || !(nm %in% c("stars","starsfdr","labelstar"))
    }, logical(1))
    p$layers <- p$layers[keep]
  }

  # Add r (below) and your chosen stars (above)
  p +
    ggplot2::geom_text(
      data = d, ggplot2::aes(label = label_r),
      inherit.aes = TRUE, size = r_size, color = r_color,
      position = ggplot2::position_nudge(y = r_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = d, ggplot2::aes(label = label_star),
      inherit.aes = TRUE, size = star_size, color = star_color,
      position = ggplot2::position_nudge(y = star_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    )
}
