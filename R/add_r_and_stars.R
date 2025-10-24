#' Add r-values and significance stars to a correlations heatmap (pass `res`)
#'
#' Works with the list returned by PlotCorrelationsHeatmap(). Chooses the
#' appropriate ggplot inside `res` based on `star_from`, draws r-values slightly
#' below center, and draws significance stars slightly above center.
#'
#' @param res The list returned by PlotCorrelationsHeatmap()
#' @param star_from One of:
#'   - "existing" : use whatever the chosen plot already mapped as stars
#'   - "raw"      : use `stars` if present, else compute from `P`
#'   - "fdr"      : use `stars_FDR` if present, else compute from `P_adj`
#'   - "P"        : compute stars from `P` cutpoints
#'   - "P_adj"    : compute stars from `P_adj` cutpoints
#'   - "column"   : use a custom column name given by `star_col`
#'   Notes:
#'     * "raw" and "P" prefer res$Unadjusted$plot; "fdr" and "P_adj" prefer res$FDRCorrected$plot.
#'     * "existing" picks the plot that already mapped a stars column (FDR preferred).
#' @param star_col  Column name to use when star_from = "column"
#' @param r_var     Column name for correlation values (default "R")
#' @param r_digits  Decimal places for r labels
#' @param r_size    Text size for r labels
#' @param r_color   Color for r labels
#' @param r_nudge_y Vertical nudge for r (negative = below center)
#' @param star_size Text size for star labels
#' @param star_color Color for star labels
#' @param star_nudge_y Vertical nudge for stars (positive = above center)
#' @param remove_existing_stars If TRUE, remove any pre-existing star text layers
#' @param p_breaks Cutpoints for p→stars (default c(-Inf, .001, .01, .05, Inf))
#' @param p_labels Labels for those cutpoints (default c("***","**","*",""))
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

  # helpers -------------------------------------------------------------
  is_text <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")
  label_to_name <- function(expr) {
    if (is.null(expr)) return(NA_character_)
    s <- tryCatch({
      if (requireNamespace("rlang", quietly = TRUE)) rlang::as_label(expr) else paste(deparse(expr), collapse = "")
    }, error = function(e) paste(deparse(expr), collapse = ""))
    # normalize: drop backticks, after_stat(), .data$, punctuation, case
    s <- tolower(s)
    s <- gsub("^after_stat\\(|\\)$", "", s)
    s <- gsub("^\\.data\\$", "", s)
    s <- gsub("[^a-z0-9_]", "", s)
    s
  }

  plot_has_star_map <- function(p) {
    if (!inherits(p, "ggplot")) return(FALSE)
    any(vapply(p$layers, function(ly) {
      if (!is_text(ly)) return(FALSE)
      nm <- label_to_name(ly$mapping$label)
      !is.na(nm) && nm %in% c("stars", "starsfdr")
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
      if (has_fdr && plot_has_star_map(res$FDRCorrected$plot)) return(res$FDRCorrected$plot)
      if (has_raw && plot_has_star_map(res$Unadjusted$plot))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from == "column") {
      if (is.null(star_col)) stop("star_from='column' requires `star_col`.")
      if (has_fdr && star_col %in% names(res$FDRCorrected$plot$data)) return(res$FDRCorrected$plot)
      if (has_raw && star_col %in% names(res$Unadjusted$plot$data))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    }
    stop("Couldn't find a valid ggplot inside `res`.")
  }

  p <- pick_plot()
  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # existing stars column used by this plot (if any)
  existing_star_col <- {
    nms <- vapply(p$layers, function(ly) if (is_text(ly)) label_to_name(ly$mapping$label) else NA_character_, character(1))
    nms <- nms[!is.na(nms)]
    cand <- nms[nms %in% c("stars","starsfdr")]
    if (length(cand)) ifelse(cand[1] == "stars", "stars", "stars_FDR") else NA_character_
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
      if (!star_col %in% names(d)) stop("`star_col` not found in plot data.")
      return(as.character(d[[star_col]]))
    }
    rep(NA_character_, nrow(d))
  }

  # labels ---------------------------------------------------------------
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- pick_stars()

  # strong, idempotent cleanup (catches "stars", "`stars`", ".data$stars", etc.)
  if (remove_existing_stars && length(p$layers)) {
    drop_names <- c("stars","starsfdr","labelstar","labelr","labelcomb")
    keep <- vapply(p$layers, function(ly) {
      if (!is_text(ly)) return(TRUE)
      nm <- label_to_name(ly$mapping$label)
      if (is.na(nm)) return(TRUE)
      !(nm %in% drop_names)
    }, logical(1))
    p$layers <- p$layers[keep]
  }

  # compose --------------------------------------------------------------
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
